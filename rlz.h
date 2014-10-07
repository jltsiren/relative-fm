#ifndef _RELATIVE_FM_RLZ_H
#define _RELATIVE_FM_RLZ_H


#include <vector>

#include <sdsl/vectors.hpp>


using namespace sdsl;

//------------------------------------------------------------------------------

const uint64_t RLZ_SA_SAMPLE_RATE = 128;

/*
  Find a LZ77 parsing of 'text' relative to 'reference'. The parsing is written into
  the last three parameters; each phrase is reference[start, start + length - 2], followed
  by mismatch.
*/
void relativeLZ(const bit_vector& text, const bit_vector& reference,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, bit_vector& mismatches);

// This version is several times slower but requires much less memory.
void relativeLZSuccinct(const bit_vector& text, const bit_vector& reference,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, bit_vector& mismatches);

//------------------------------------------------------------------------------

const uint64_t INCBWT_BLOCK_SIZE = 16 * 1048576;

/*
  Input:    Bit sequence bwt[0, n-1], with bwt[n] = 0 as an endmarker.
  Output:   BWT(bwt + $) stored in the same bit sequence. The endmarker is encoded with a 0-bit.
  Returns:  Position in bwt that contains the endmarker.
*/
uint64_t incrementalBWT(bit_vector& bwt, uint64_t block_size = INCBWT_BLOCK_SIZE);

/*
  Input:    BWT, the position of the endmarker, and the desired sample rate.
            The endmarker is assumed to be encoded with a 0-bit.
  Output:   sa_samples contains SA[i * sample_rate + 1] for all i.
*/
void sampleSA(bit_vector& bwt, bit_vector::rank_1_type& bwt_rank, uint64_t endmarker,
  int_vector<0>& sa_samples, uint64_t sample_rate);

inline uint64_t
_LF1(bit_vector::rank_1_type& rank, uint64_t pos, uint64_t zeros)
{
  return zeros + 1 + rank(pos);
}

// The endmarker is encoded as a 0-bit, which should be ignored in rank_0.
inline uint64_t
_LF0(bit_vector::rank_1_type& rank, uint64_t pos, uint64_t endmarker)
{
  return 1 + pos - rank(pos) - (pos > endmarker ? 1 : 0);
}

//------------------------------------------------------------------------------

#endif // _RELATIVE_FM_RLZ_H
