/*
  Copyright (c) 2015 Genome Research Ltd.
  Copyright (c) 2014 Jouni Siren

  Author: Jouni Siren <jouni.siren@iki.fi>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#ifndef _RELATIVE_FM_RLZ_H
#define _RELATIVE_FM_RLZ_H

#include <vector>

#include <sdsl/suffix_arrays.hpp>

#include "support.h"

using namespace sdsl;

namespace relative
{

//------------------------------------------------------------------------------

/*
  Find a LZ77 parsing of 'text' relative to 'reference'. The parsing is written into
  the last three parameters; each phrase is reference[start, start + length - 2], followed
  by mismatch. Use bv_fmi constructor to build index for the reverse reference.
*/
void relativeLZSuccinct(const bit_vector& text, const bit_vector& reference,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, bit_vector& mismatches);

struct bv_fmi;

// Use this to reuse an existing index for the reverse reference.
void relativeLZSuccinct(const bit_vector& text, const bv_fmi& reference,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, bit_vector& mismatches);

/*
  This version uses temporary files on disk. Alternatively, use reverseIndex(reference, csa)
  to build index for the reverse reference, and then call the template version. The sequence
  and the reference may contain character value 0 or character value 1, but not both.
*/
void relativeLZ(const int_vector<8>& text, const int_vector<8>& reference,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, int_vector<8>& mismatches);

//------------------------------------------------------------------------------

template<class IntVector, bool differential>
struct CharAt
{
  inline static uint64_t at(const IntVector& seq, uint64_t i)
  {
    if(i >= seq.size()) { return 0; }
    return seq[i];
  }
};

template<class IntVector>
struct CharAt<IntVector, true>
{
  inline static uint64_t at(const IntVector& seq, uint64_t i)
  {
    if(i >= seq.size()) { return 0; }
    uint64_t prev = (i > 0 ? seq[i - 1] : 0), curr = seq[i];
    return DiffEncoderNZ::encode(curr, prev);
  }
};

// FIXME Handle this in a more general way.
const uint64_t MAX_RLZ_PHRASE_LENGTH = 1024;

/*
  This version is intended for integer sequences. It parses the reference using a prebuilt
  suffix array. The text must not contain character value 0, while the reference has it only
  as an implicit endmarker.
  If mismatches == 0, the algorithm does not output them.
  If differential == true, the parsing is based on differential values instead of absolute
  values.
*/
template<class IntVector, bool differential>
void
relativeLZ(const IntVector& text, const IntVector& reference, const int_vector<0>& sa,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, std::vector<uint64_t>* mismatches)
{
  if(text.size() == 0) { return; }

  starts.clear(); lengths.clear();
  if(mismatches != 0) { mismatches->clear(); }
  uint64_t text_pos = 0;
  while(text_pos < text.size())
  {
    uint64_t sp = 0, ep = sa.size() - 1, matched = 0;
    uint64_t limit = std::min(text.size() - text_pos, MAX_RLZ_PHRASE_LENGTH);
    while(matched < limit)
    {
      uint64_t low = sp, high = ep, next = CharAt<IntVector, differential>::at(text, text_pos + matched);
      while(low < high) // Find the first suffix that matches the next character.
      {
        uint64_t mid = low + (high - low) / 2;
        uint64_t val = CharAt<IntVector, differential>::at(reference, sa[mid] + matched);
        if(val < next) { low = mid + 1; }
        else if(val > next) { high = mid - 1; }
        else { high = mid; }
      }
      if(CharAt<IntVector, differential>::at(reference, sa[low] + matched) != next) { break; }
      sp = low;

      high = ep;
      while(low < high) // Find the last suffix that matches the next character.
      {
        uint64_t mid = low + (high + 1 - low) / 2;
        uint64_t val = CharAt<IntVector, differential>::at(reference, sa[mid] + matched);
        if(val > next) { high = mid - 1; }
        else { low = mid; }
      }
      ep = high; matched++;
    }

    starts.push_back(sa[sp]); // FIXME: Find the match closest to text_pos?
    if(text_pos + matched < text.size()) { matched++; } // Add the mismatch.
    lengths.push_back(matched);
    if(mismatches != 0) { mismatches->push_back(text[text_pos + matched - 1]); }
    text_pos += matched;
  }
}

//------------------------------------------------------------------------------

/*
  The sequence and the reference may contain character value 0 or character value 1,
  but not both.
*/
template<class IntVector, class CSA>
void
relativeLZ(const IntVector& text, const CSA& reference,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, IntVector& mismatches)
{
  if(text.size() == 0) { return; }

  // Use backward searching on the BWT of the reverse reference to parse the text.
  uint64_t text_pos = 0;
  starts.clear(); lengths.clear();
  std::vector<typename IntVector::value_type> char_buffer;
  while(text_pos < text.size())
  {
    uint64_t sp = 0, ep = reference.size() - 1;
    uint64_t len = 0; // We have matched len characters in the reverse reference.
    while(text_pos + len < text.size())
    {
      // Does the current character exist in the reference?
      if(reference.char2comp[text[text_pos + len]] == 0) { break; }
      uint64_t new_sp = 0, new_ep = 0;
      backward_search(reference, sp, ep, text[text_pos + len], new_sp, new_ep);
      if(new_sp > new_ep) { break; }
      else { sp = new_sp; ep = new_ep; len++; }
    }
    uint64_t reverse_pos = reference[sp]; // Position of the reverse pattern in reverse reference.
    starts.push_back(reference.size() - reverse_pos - len - 1); // Convert into actual pattern position.
    if(text_pos + len < text.size()) { len++; } // Add the mismatching character.
    lengths.push_back(len);
    char_buffer.push_back(text[text_pos + len - 1]);
    text_pos += len;
  }
  mismatches.resize(char_buffer.size());
  for(uint64_t i = 0; i < char_buffer.size(); i++) { mismatches[i] = char_buffer[i]; }
}

//------------------------------------------------------------------------------

/*
  If the sequence contains character value 0, it is converted into 1.
  Therefore the sequence should not contain both 0s and 1s.
*/

template<class IntVector, class CSA>
void
reverseIndex(const IntVector& sequence, CSA& csa)
{
  cache_config config;
  std::string filename = tmp_file(config, "text");
  std::ofstream out(filename.c_str(), std::ios_base::binary);
  if(!out)
  {
    std::cerr << "reverseIndex(): Cannot open temporary file " << filename << " for writing" << std::endl;
    return;
  }
  IntVector reverse(sequence.size());
  for(uint64_t i = 0; i < sequence.size(); i++)
  {
    reverse[i] = std::max(sequence[sequence.size() - 1 - i], (typename IntVector::value_type)1);
  }
  reverse.serialize(out); util::clear(reverse);
  out.close();

  construct(csa, filename, 0);
  sdsl::remove(filename);
}

//------------------------------------------------------------------------------

/*
  An FM-index for the reverse of a bitvector.
*/
struct bv_fmi
{
public:
  typedef std::pair<uint64_t, uint64_t> range_type;

  const static uint64_t DEFAULT_BLOCK_SIZE  = 64 * MEGABYTE;
  const static uint64_t DEFAULT_SAMPLE_RATE = 127;  // Should be prime with SA order sampling.

  explicit bv_fmi(const bit_vector& source,
    uint64_t block_size = DEFAULT_BLOCK_SIZE, uint64_t _sample_rate = DEFAULT_SAMPLE_RATE);
  explicit bv_fmi(std::istream& in);

  uint64_t serialize(std::ostream& out);

  bit_vector              bwt;
  bit_vector::rank_1_type rank;
  uint64_t                zeros;      // Number of 0-bits, excluding the endmarker.
  uint64_t                endmarker;  // Position of the endmarker, encoded by a 0-bit.

  int_vector<0>           sa_samples; // Will contain SA[i * sample_rate + 1] for all i.
  uint64_t                sample_rate;

  inline uint64_t LF1(uint64_t pos) const
  {
    return this->zeros + 1 + this->rank(pos);
  }

  inline uint64_t LF0(uint64_t pos) const
  {
    return 1 + pos - this->rank(pos) - (pos > this->endmarker ? 1 : 0);
  }

  inline range_type bitRange(bool bit) const
  {
    return (bit ? range_type(this->zeros + 1, this->bwt.size() - 1) : range_type(1, this->zeros));
  }

  // pos % sample_rate is assumed to be 1.
  inline uint64_t sampleAt(uint64_t pos) const
  {
    return this->sa_samples[pos / this->sample_rate];
  }

private:
  void incrementalBWT(uint64_t block_size);
  void lastBWTBlock(uint64_t offset);
  void prevBWTBlock(uint64_t offset, uint64_t block_size);
  void sampleSA();
};

//------------------------------------------------------------------------------

/*
  This struct encodes the starting positions of phrases in the reference relative to
  the starting positions in text. The init() function is compatible with the RLZ
  parser, taking phrase lengths instead of text positions as input.
*/
struct relative_encoder
{
  typedef uint64_t size_type;

  int_vector<0>           values; // Phrase starts encoded as (ref_pos - text_pos).
  bit_vector              rle;    // rle[i] is set, if phrase i is stored.
  bit_vector::rank_1_type rank;

  relative_encoder();
  relative_encoder(const relative_encoder& r);
  relative_encoder(relative_encoder&& r);
  relative_encoder& operator=(const relative_encoder& r);
  relative_encoder& operator=(relative_encoder&& r);

  template<class Container>
  void init(const Container& ref_pos, const Container& lengths)
  {
    util::clear(this->values); util::clear(this->rle); util::clear(this->rank);

    // Check whether run-length encoding can help.
    size_type text_pos = 0, prev = ~(size_type)0, rle_max = 0, direct_max = 0, runs = 0;
    for(size_type i = 0; i < ref_pos.size(); i++)
    {
      size_type temp = encode(ref_pos[i], text_pos);
      if(i == 0 || (temp != prev && lengths[i] > 1))
      {
        rle_max = std::max(rle_max, temp);
        prev = temp; runs++;
      }
      direct_max = std::max(direct_max, temp);
      text_pos += lengths[i];
    }
    double rle_bits = runs * bitlength(rle_max) + ref_pos.size() * 1.25;
    double direct_bits = ref_pos.size() * bitlength(direct_max);

    if(rle_bits >= direct_bits) // Use direct encoding.
    {
#ifdef VERBOSE_STATUS_INFO
      std::cout << "Using direct encoding for starting positions." << std::endl;
#endif
      int_vector<0> buffer(ref_pos.size(), 0, bitlength(direct_max));
      text_pos = 0;
      for(size_type i = 0; i < ref_pos.size(); i++)
      {
        buffer[i] = encode(ref_pos[i], text_pos);
        text_pos += lengths[i];
      }
      util::assign(this->values, buffer);
    }
    else  // Use run-length encoding.
    {
#ifdef VERBOSE_STATUS_INFO
      std::cout << "Using run-length encoding for starting positions." << std::endl;
#endif
      int_vector<0> buffer(runs, 0, bitlength(rle_max));
      util::assign(this->rle, bit_vector(ref_pos.size(), 0));
      text_pos = 0; prev = 0; runs = 0;
      for(size_type i = 0; i < ref_pos.size(); i++)
      {
        size_type temp = encode(ref_pos[i], text_pos);
        if(i == 0 || (temp != prev && lengths[i] > 1))
        {
          buffer[runs] = temp;
          prev = temp; runs++;
          this->rle[i] = 1;
        }
        text_pos += lengths[i];
      }
      util::assign(this->values, buffer);
      util::init_support(this->rank, &(this->rle));
    }
  }

  size_type reportSize() const;
  void load(std::istream& input);
  size_type serialize(std::ostream& output, structure_tree_node* v = nullptr, std::string name = "") const;

  // Returns a reference position for given text position relative to given phrase.
  inline size_type decode(size_type phrase, size_type text_pos) const
  {
    size_type temp = (this->rle.size() > 0 ? this->values[this->rank(phrase + 1) - 1] : this->values[phrase]);
    if(temp & 1) { return text_pos - (temp + 1) / 2; }
    return text_pos + temp / 2;
  }

  // Encodes (ref_pos - text_pos) as unsigned integer.
  inline static size_type encode(size_type ref_pos, size_type text_pos)
  {
    if(ref_pos >= text_pos) { return 2 * (ref_pos - text_pos); }
    return 2 * (text_pos - ref_pos) - 1;
  }

private:
  void copy(const relative_encoder& r);
};

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_RLZ_H
