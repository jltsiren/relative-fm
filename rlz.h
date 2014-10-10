#ifndef _RELATIVE_FM_RLZ_H
#define _RELATIVE_FM_RLZ_H


#include <vector>

#include <sdsl/vectors.hpp>


namespace sdsl
{

//------------------------------------------------------------------------------

/*
  Find a LZ77 parsing of 'text' relative to 'reference'. The parsing is written into
  the last three parameters; each phrase is reference[start, start + length - 2], followed
  by mismatch.
*/
void relativeLZ(const bit_vector& text, const bit_vector& reference,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, bit_vector& mismatches);

struct bv_fmi;

// Use this if BWT and SA samples have already been built.
void relativeLZ(const bit_vector& text, const bv_fmi& reference,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, bit_vector& mismatches);

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

} // namespace sdsl

#endif // _RELATIVE_FM_RLZ_H
