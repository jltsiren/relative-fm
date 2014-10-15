#ifndef _RELATIVE_FM_RLZ_H
#define _RELATIVE_FM_RLZ_H


#include <vector>

#include <sdsl/suffix_arrays.hpp>
#include <sdsl/vectors.hpp>


namespace sdsl
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
    std::cerr << "reverseIndex(): Cannot open temporary file " << filename << " for writing!" << std::endl;
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

  const static uint64_t DEFAULT_BLOCK_SIZE  = 64 * 1048576;
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
  This struct encodes an ordered set of items stored in an ordered set of blocks.
  Both items and blocks are 0-based, and some of the blocks may be empty.
*/
struct rlz_helper
{
  sd_vector<>                 v;
  sd_vector<>::rank_1_type    v_rank;
  sd_vector<>::select_1_type  v_select;

  bit_vector                  nonzero;
  bit_vector::rank_1_type     nz_rank;
  bit_vector::select_1_type   nz_select;

  rlz_helper();
  rlz_helper(const rlz_helper& r);
  rlz_helper(rlz_helper&& r);
  rlz_helper& operator=(const rlz_helper& r);
  rlz_helper& operator=(rlz_helper&& r);

  void init(const std::vector<uint64_t>& values);

  inline uint64_t blockFor(uint64_t item) const
  {
    uint64_t block = this->v_rank(item);
    if(this->nonzero.size() > 0) { block = this->nz_select(block + 1); }
    return block;
  }

  inline uint64_t itemsAfter(uint64_t block) const
  {
    block++;  // Convert to 1-based blocks.
    if(this->nonzero.size() > 0)
    {
      block = this->nz_rank(block);
      if(block == 0) { return 0; }
    }
    return this->v_select(block) + 1;
  }

  inline bool isLast(uint64_t item) const
  {
    return v[item];
  }

  // This assumes that all blocks have non-zero size.
  inline uint64_t blockSize(uint64_t block) const
  {
    uint64_t res = this->v_select(block + 1) + 1;
    if(block > 0) { res -= this->v_select(block) + 1; }
    return res;
  }

  uint64_t reportSize() const;
  void load(std::istream& input);
  uint64_t serialize(std::ostream& output) const;

private:
  void copy(const rlz_helper& r);
};

//------------------------------------------------------------------------------

} // namespace sdsl

#endif // _RELATIVE_FM_RLZ_H
