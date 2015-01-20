#ifndef _RELATIVE_FM_RELATIVE_FM_H
#define _RELATIVE_FM_RELATIVE_FM_H


#include <fstream>
#include <iostream>

#include <sdsl/csa_alphabet_strategy.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "utils.h"

//------------------------------------------------------------------------------

/*
  A simple FM-index for byte alphabet.

  The storage format is an int_vector<8> containing the BWT (with a contiguous alphabet)
  in base_name.bwt, and a byte_alphabet mapping the original alphabet to the contiguous
  alphabet in base_name.alpha.

  RankStructure must support the following operations:
    - Queries: [i], rank(i, c)
    - Constructor: (), (int_vector_buffer<8>, uint64_t)
    - Basic operations: swap(), size()
    - With the structure as a parameter: size_in_bytes(bwt), serialize(bwt, .., ..)
*/
template<class RankStructure = bwt_type>
class SimpleFM
{
public:
  explicit SimpleFM(const std::string& base_name)
  {
    {
      int_vector_buffer<8> buffer(base_name + BWT_EXTENSION);
      RankStructure temp(buffer, buffer.size());
      this->bwt.swap(temp);
    }
    {
      std::string filename = base_name + ALPHA_EXTENSION;
      std::ifstream in(filename.c_str(), std::ios_base::binary);
      if(!in)
      {
        std::cerr << "SimpleFM()::SimpleFM(): Cannot open alphabet file " << filename << std::endl;
        return;
      }
      this->alpha.load(in);
      in.close();
    }
  }

  ~SimpleFM()
  {
  }

  uint64_t size() const { return this->bwt.size(); }
  uint64_t sequences() const { return this->alpha.C[1]; }

  uint64_t reportSize(bool print = false) const
  {
    uint64_t bytes = size_in_bytes(this->bwt) + size_in_bytes(this->alpha);

    if(print)
    {
      printSize("Simple FM", bytes, this->size());
      std::cout << std::endl;
    }

    return bytes;
  }

  void writeTo(const std::string& base_name) const
  {
    {
      std::string filename = base_name + ALPHA_EXTENSION;
      std::ofstream output(filename.c_str(), std::ios_base::binary);
      if(!output)
      {
        std::cerr << "SimpleFM::writeTo(): Cannot open alphabet file " << filename << std::endl;
        return;
      }
      this->alpha.serialize(output); output.close();
    }

    {
      std::string filename = base_name + BWT_EXTENSION;
      int_vector_buffer<8> buffer(filename, std::ios::out);
      for(uint64_t i = 0; i < this->size(); i++)
      {
        buffer[i] = this->alpha.comp2char[this->bwt[i]];
      }
      buffer.close();
    }
  }

  template<class Iter> range_type find(Iter begin, Iter end) const
  {
    range_type res(0, this->size() - 1);
    while(begin != end)
    {
      --end;
      res = LF(this->bwt, this->alpha, res, *end);
      if(length(res) == 0) { return range_type(1, 0); }
    }
    return res;
  }

  template<class ByteVector>
  void extractBWT(range_type range, ByteVector& buffer) const
  {
    if(isEmpty(range) || range.second >= this->size()) { return; }
    buffer.resize(length(range));
    for(uint64_t i = 0; i < buffer.size(); i++)
    {
      buffer[i] = this->alpha.comp2char[this->bwt[range.first + i]];
    }
  }

  template<class ByteVector>
  void extractBWT(ByteVector& buffer) const
  {
    this->extractBWT(range_type(0, this->size() - 1), buffer);
  }

  RankStructure bwt;
  alphabet_type alpha;

private:
  SimpleFM();
  SimpleFM(const SimpleFM&);
  SimpleFM(SimpleFM&&);
  SimpleFM& operator=(const SimpleFM&);
  SimpleFM& operator=(SimpleFM&&);
};

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

  RelativeFM(const SimpleFM<>& ref, const SimpleFM<>& seq, bool print = false);
  RelativeFM(const SimpleFM<>& ref, const std::string& base_name);
  RelativeFM(const SimpleFM<>& ref, std::istream& input);
  ~RelativeFM();

  uint64_t size() const { return this->m_size; }
  uint64_t sequences() const { return this->alpha.C[1]; }

  uint64_t reportSize(bool print = false) const;
  void writeTo(const std::string& base_name) const;
  void writeTo(std::ostream& output) const;

  template<class Iter> range_type find(Iter begin, Iter end) const;

#ifdef USE_SPARSE_BITVECTORS
  typedef sd_vector<> vector_type;
#else
  typedef rrr_vector<63> vector_type;
#endif

  const SimpleFM<>&           reference;
  bwt_type                    ref_minus_lcs, seq_minus_lcs;
  vector_type                 ref_lcs, seq_lcs;
  alphabet_type               alpha;
  uint64_t                    m_size;

  vector_type::rank_1_type    seq_rank;
#ifdef USE_SPARSE_BITVECTORS
  vector_type::select_0_type  ref_select;
#else
  vector_type::select_1_type  ref_select;
#endif

private:
  void loadFrom(std::istream& input);
  void buildRankSelect();

  /*
    The straightforward implementation counts the number of c's up to position i.
    To get the proper result, we need to check whether position i is in LCS or its complement,
    and ignore the last position when doing rank in that sequence.
    Note that c is a real character.
  */
  inline uint64_t rank(uint64_t i, uint8_t c) const
  {
    uint64_t res = 0;
    uint8_t ref_c = this->reference.alpha.char2comp[c], seq_c = this->alpha.char2comp[c];

#ifdef USE_SPARSE_BITVECTORS  // LCS is marked with 0-bits.
    uint64_t lcs_bits = i + 1 - this->seq_rank.rank(i + 1); // Number of LCS bits up to i in seq.
    bool check_lcs = (this->seq_lcs[i] == 0); // Is position i in LCS.
#else
    uint64_t lcs_bits = this->seq_rank.rank(i + 1); // Number of LCS bits up to i in seq.
    bool check_lcs = (this->seq_lcs[i] == 1); // Is position i in LCS.
#endif
    if(lcs_bits < i + 1) // There are i + 1 - lcs_bits non-LCS bits in seq.
    {
      res += this->seq_minus_lcs.rank(i + check_lcs - lcs_bits, seq_c);
    }
    if(lcs_bits > 0)
    {
      uint64_t ref_pos = this->ref_select.select(lcs_bits);  // Select is 1-based.
      res += this->reference.bwt.rank(ref_pos + 1 - check_lcs, ref_c);
      if(lcs_bits < ref_pos + 1) // At least one non-LCS bit in ref.
      {
        res -= this->ref_minus_lcs.rank(ref_pos + 1 - lcs_bits, ref_c);
      }
    }

    return res;
  }

  RelativeFM();
  RelativeFM(const RelativeFM&);
  RelativeFM(RelativeFM&&);
  RelativeFM& operator=(const RelativeFM&);
  RelativeFM& operator==(RelativeFM&);
};

template<class Iter>
range_type
RelativeFM::find(Iter begin, Iter end) const
{
  range_type res(0, this->size() - 1);
  while(begin != end)
  {
    --end;
    uint64_t begin = cumulative(this->alpha, *end);
    res.first = begin + this->rank(res.first, *end);
    res.second = begin + this->rank(res.second + 1, *end) - 1;
    if(length(res) == 0) { return range_type(1, 0); }
  }
  return res;
}

std::vector<std::pair<int, int> >
greedyLCS(const SimpleFM<>& ref, const SimpleFM<>& seq, range_type ref_range, range_type seq_range, bool onlyNs);

std::pair<bit_vector, bit_vector>
alignBWTs(const SimpleFM<>& ref, const SimpleFM<>& seq, uint64_t block_size, uint max_depth, uint64_t& lcs, bool print);

//------------------------------------------------------------------------------

#endif // _RELATIVE_FM_RELATIVE_FM_H
