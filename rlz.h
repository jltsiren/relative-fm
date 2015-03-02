#ifndef _RELATIVE_FM_RLZ_H
#define _RELATIVE_FM_RLZ_H

#include <vector>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/vectors.hpp>

#include "utils.h"

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
  These versions use temporary files on disk. Use reverseIndex(reference, csa)
  to build index for the reverse reference. The sequence and the reference may contain
  character value 0 or character value 1, but not both.

  For the int_vector<0> variant, it is advisable to use util::bit_compress(mismatches)
  afterwards.
*/
void relativeLZ(const int_vector<8>& text, const int_vector<8>& reference,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, int_vector<8>& mismatches);

void relativeLZ(const int_vector<0>& text, const int_vector<0>& reference,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, int_vector<0>& mismatches);

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

#ifdef VERBOSE_STATUS_INFO
  std::cout << "Parsed the text as " << starts.size() << " phrases." << std::endl;
#endif
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
    uint64_t text_pos = 0, prev = ~(uint64_t)0, rle_max = 0, direct_max = 0, runs = 0;
    for(uint64_t i = 0; i < ref_pos.size(); i++)
    {
      uint64_t temp = encode(ref_pos[i], text_pos);
      if(i == 0 || (temp != prev && lengths[i] > 1))
      {
        rle_max = std::max(rle_max, temp);
        prev = temp; runs++;
      }
      direct_max = std::max(direct_max, temp);
      text_pos += lengths[i];
    }
    double rle_bits = runs * (bits::hi(rle_max) + 1) + ref_pos.size() * 1.25;
    double direct_bits = ref_pos.size() * (bits::hi(direct_max) + 1);

    if(rle_bits >= direct_bits) // Use direct encoding.
    {
#ifdef VERBOSE_STATUS_INFO
      std::cout << "Using direct encoding for starting positions." << std::endl;
#endif
      this->values.width(bits::hi(direct_max) + 1); this->values.resize(ref_pos.size());
      text_pos = 0;
      for(uint64_t i = 0; i < ref_pos.size(); i++)
      {
        this->values[i] = encode(ref_pos[i], text_pos);
        text_pos += lengths[i];
      }
    }
    else  // Use run-length encoding.
    {
#ifdef VERBOSE_STATUS_INFO
      std::cout << "Using run-length encoding for starting positions." << std::endl;
#endif
      this->values.width(bits::hi(rle_max) + 1); this->values.resize(runs);
      this->rle.resize(ref_pos.size()); util::set_to_value(this->rle, 0);
      text_pos = 0; prev = 0; runs = 0;
      for(uint64_t i = 0; i < ref_pos.size(); i++)
      {
        uint64_t temp = encode(ref_pos[i], text_pos);
        if(i == 0 || (temp != prev && lengths[i] > 1))
        {
          this->values[runs] = temp;
          prev = temp; runs++;
          this->rle[i] = 1;
        }
        text_pos += lengths[i];
      }
      util::init_support(this->rank, &(this->rle));
    }
  }

  uint64_t reportSize() const;
  void load(std::istream& input);
  uint64_t serialize(std::ostream& output) const;

  // Returns a reference position for given text position relative to given phrase.
  inline uint64_t decode(uint64_t phrase, uint64_t text_pos) const
  {
    uint64_t temp = (this->rle.size() > 0 ? this->values[this->rank(phrase + 1) - 1] : this->values[phrase]);
    if(temp & 1) { return text_pos - (temp + 1) / 2; }
    return text_pos + temp / 2;
  }

  // Encodes (ref_pos - text_pos) as unsigned integer.
  inline static uint64_t encode(uint64_t ref_pos, uint64_t text_pos)
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
