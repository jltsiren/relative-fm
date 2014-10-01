#ifndef _RELATIVE_FM_RLZ_VECTOR_H
#define _RELATIVE_FM_RLZ_VECTOR_H


#include <fstream>
#include <iostream>

#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>

#include "utils.h"

//------------------------------------------------------------------------------

class RLZVector
{
public:
  RLZVector(const bit_vector& text, const bit_vector& _reference,
    const bit_vector::rank_1_type& _ref_rank,
    const bit_vector::select_1_type& _ref_select);
  RLZVector(std::ifstream& input, const bit_vector& _reference,
    const bit_vector::rank_1_type& _ref_rank,
    const bit_vector::select_1_type& _ref_select);
  ~RLZVector();

  uint64_t reportSize() const;
  void writeTo(std::ofstream& output) const;

  uint64_t size() const { return this->lengths.size() - 1; }
  uint64_t items() const { return this->ones.size() - 1; }

  /*
    These follow SDSL conventions.
  */
  uint64_t rank(uint64_t i) const;
  uint64_t select(uint64_t i) const;
  bool operator[](uint64_t i) const;

private:
  const bit_vector&                reference;
  const bit_vector::rank_1_type&   ref_rank;
  const bit_vector::select_1_type& ref_select;

  int_vector<0> phrases;

  sd_vector<> lengths;
  sd_vector<>::rank_1_type length_rank;
  sd_vector<>::select_1_type length_select;

  sd_vector<> ones;
  sd_vector<>::rank_1_type one_rank;
  sd_vector<>::select_1_type one_select;

  void buildRankSelect();

  // Counts the number of 1-bits in reference[ref_pos, ref_pos + phrase_length - 1].
  inline uint64_t oneBits(uint64_t ref_pos, uint64_t phrase_length) const
  {
    if(phrase_length == 0) { return 0; }
    return this->ref_rank(ref_pos + phrase_length) - this->ref_rank(ref_pos);
  }

  // Returns the offset of the i'th 1-bit in reference[ref_pos, ...] (i is 1-based).
  inline uint64_t findBit(uint64_t ref_pos, uint64_t i) const
  {
    return this->ref_select(this->ref_rank(ref_pos) + i) - ref_pos;
  }

  inline uint64_t bitAt(uint64_t ref_pos, uint64_t offset) const
  {
    return this->reference[ref_pos + offset];
  }

  RLZVector();
  RLZVector(const RLZVector&);
  RLZVector& operator=(const RLZVector&);
};

//------------------------------------------------------------------------------

#endif // _RELATIVE_FM_RLZ_VECTOR_H
