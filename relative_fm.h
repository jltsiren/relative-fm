#ifndef _RELATIVE_FM_RELATIVE_FM_H
#define _RELATIVE_FM_RELATIVE_FM_H


#include <fstream>
#include <iostream>

#include <sdsl/csa_alphabet_strategy.hpp>
#include <sdsl/wavelet_trees.hpp>


#include "utils.h"

/*
  There is probably some double counting in the reportSize() functions. SDSL function
  size_in_bytes() apparently returns the size of a serialized object, not the size of
  the object in memory. It is also unclear how to best handle the case, where objects
  contain other objects instead of pointers or references.
*/

//------------------------------------------------------------------------------

class SimpleFM
{
public:
  const static std::string EXTENSION; // .bwt

  explicit SimpleFM(const std::string& base_name);
  explicit SimpleFM(std::ifstream& input);
  ~SimpleFM();

  uint64_t reportSize(bool print = false) const;
  void writeTo(const std::string& base_name) const;
  void writeTo(std::ofstream& output) const;

  template<class Iter> range_type find(Iter begin, Iter end) const;

  bwt_type      bwt;
  alphabet_type alpha;

private:
  void loadFrom(std::ifstream& input);

  SimpleFM();
  SimpleFM(const SimpleFM&);
  SimpleFM& operator=(const SimpleFM&);
};

template<class Iter>
range_type
SimpleFM::find(Iter begin, Iter end) const
{
  range_type res(0, this->bwt.size() - 1);
  while(begin != end)
  {
    if(!hasChar(this->alpha, *begin)) { return range_type(1, 0); }
    res = LF(this->bwt, this->alpha, res, *begin);
    if(isEmpty(res)) { return range_type(1, 0); }
    ++begin;
  }
  return res;
}

//------------------------------------------------------------------------------

class RelativeFM
{
public:
  const static std::string EXTENSION;          // .rfm
  const static uint64_t    BLOCK_SIZE = 1024;  // Split the BWTs into blocks of this size or less.
  const static uint        MAX_DEPTH  = 32;    // Maximum length of a pattern used to split the BWTs.

  /*
    Maximum diagonal in LCS computation. If further diagonals would be needed, only the most frequent
    character in the ranges will be matched. Maximal memory usage will be around
    4 * MAX_D * MAX_D bytes (+ SimpleFM for the reference and the target sequence).
  */
  const static int MAX_D = 50000;

  RelativeFM(const SimpleFM& ref, const SimpleFM& seq, bool print = false);
  RelativeFM(const SimpleFM& ref, const std::string& base_name);
  RelativeFM(const SimpleFM& ref, std::ifstream& input);
  ~RelativeFM();

  uint64_t reportSize(bool print = false) const;
  void writeTo(const std::string& base_name) const;
  void writeTo(std::ofstream& output) const;

  template<class Iter> range_type find(Iter begin, Iter end) const;

  /*
    The hybrid bitvector of Juha Kärkkäinen, Dominik Kempa, and Simon J. Puglisi would probably
    be good here, but it only supports vectors shorter than 2^31 bits.
  */
//  typedef hybrid_vector<> vector_type;
  typedef rrr_vector<63> vector_type;

  const SimpleFM& reference;
  bwt_type        ref_minus_lcs, seq_minus_lcs;
  vector_type     ref_lcs, seq_lcs;
  alphabet_type   alpha;
  uint64_t        size;

  vector_type::rank_1_type*   ref_rank;
  vector_type::select_1_type* ref_select;
  vector_type::rank_1_type*   seq_rank;
  vector_type::select_1_type* seq_select;

private:
  void loadFrom(std::ifstream& input);
  void buildRankSelect();

  /*
    The straightforward implementation counts the number of c's up to position i.
    To get the proper result, we need to check whether position i is in LCS or its complement,
    and ignore the last position when doing rank in that sequence.
  */
  inline uint64_t rank(uint64_t i, uint8_t c) const
  {
    uint64_t res = 0;
    uint64_t lcs_bits = this->seq_rank->rank(i + 1); // Number of LCS bits up to i in seq.
    bool check_lcs = (this->seq_lcs[i] == 1); // Is position i in LCS.
    if(lcs_bits < i + 1) // There are i + 1 - lcs_bits non-LCS bits in seq.
    {
      res += this->seq_minus_lcs.rank(i + check_lcs - lcs_bits, c);
    }
    if(lcs_bits > 0)
    {
      uint64_t ref_pos = this->ref_select->select(lcs_bits);  // Select is 1-based.
      res += this->reference.bwt.rank(ref_pos + 1 - check_lcs, c);
      if(lcs_bits < ref_pos + 1) // At least one non-LCS bit in ref.
      {
        res -= this->ref_minus_lcs.rank(ref_pos + 1 - lcs_bits, c);
      }
    }
    return res;
  }

  RelativeFM();
  RelativeFM(const RelativeFM&);
  RelativeFM& operator=(const RelativeFM&);
};

template<class Iter>
range_type
RelativeFM::find(Iter begin, Iter end) const
{
  range_type res(0, this->size - 1);
  while(begin != end)
  {
    if(!hasChar(this->alpha, *begin)) { return range_type(1, 0); }
    uint64_t pos = cumulative(this->alpha, *begin);
    res.first = pos + this->rank(res.first, *begin);
    res.second = pos + this->rank(res.second + 1, *begin) - 1;
    if(isEmpty(res)) { return range_type(1, 0); }
    ++begin;
  }
  return res;
}

std::vector<std::pair<int, int> >
greedyLCS(const bwt_type& ref, const bwt_type& seq, range_type ref_range, range_type seq_range, bool onlyNs);

std::pair<bit_vector, bit_vector>
alignBWTs(const SimpleFM& ref, const SimpleFM& seq, uint64_t block_size, uint max_depth, uint64_t& lcs, bool print);

//------------------------------------------------------------------------------

#endif // _RELATIVE_FM_RELATIVE_FM_H
