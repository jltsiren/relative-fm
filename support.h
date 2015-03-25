#ifndef _RELATIVE_FM_SUPPORT_H
#define _RELATIVE_FM_SUPPORT_H

#include "utils.h"

using namespace sdsl;

namespace relative
{

//------------------------------------------------------------------------------

// This function has been specialized for RLSequence.
template<class ByteVector>
void
characterCounts(const ByteVector& sequence, int_vector<64>& counts)
{
  for(uint64_t c = 0; c < counts.size(); c++) { counts[c] = 0; }
  for(uint64_t i = 0; i < sequence.size(); i++) { counts[sequence[i]]++; }
}

/*
  This replaces the SDSL byte_alphabet. The main improvements are:
    - The alphabet can be built from an existing sequence.
    - The comp order does not need to be the same as character order, as long as \0 is the first character.
*/

class Alphabet
{
public:
  typedef uint64_t size_type;
  const static size_type MAX_SIGMA = 256;

  Alphabet();
  Alphabet(const Alphabet& s);
  Alphabet(Alphabet&& s);
  ~Alphabet();

  /*
    ByteVector only has to support operator[]. If there is a clearly faster way for sequential
    access, function characterCounts() should be specialized.
  */
  template<class ByteVector>
  explicit Alphabet(const ByteVector& sequence) :
    char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma)
  {
    util::assign(this->m_char2comp, int_vector<8>(MAX_SIGMA, 0));
    util::assign(this->m_comp2char, int_vector<8>(MAX_SIGMA, 0));
    util::assign(this->m_C, int_vector<64>(MAX_SIGMA + 1, 0));
    this->m_sigma = 0;
    if(sequence.size() == 0) { return; }

    // Step 1: Count the occurrences and temporarily store them in m_C.
    characterCounts(sequence, this->m_C);

    // Step 2: Determine the effective alphabet and compact the alphabet.
    for(size_type i = 0; i < MAX_SIGMA; i++)
    {
      if(this->m_C[i] > 0)
      {
        this->m_char2comp[i] = this->m_sigma;
        this->m_comp2char[this->m_sigma] = i;
        this->m_C[this->m_sigma] = this->m_C[i];
        this->m_sigma++;
      }
    }
    this->m_comp2char.resize(this->m_sigma);
    this->m_C.resize(this->m_sigma + 1);

    // Step 3: Determine the cumulative counts.
    for(size_type i = this->m_sigma; i > 0; i--) { this->m_C[i] = this->m_C[i - 1]; }
    this->m_C[0] = 0;
    for(size_type i = 1; i <= this->m_sigma; i++) { this->m_C[i] += this->m_C[i - 1]; }
  }

  void swap(Alphabet& v);
  Alphabet& operator=(const Alphabet& v);
  Alphabet& operator=(Alphabet&& v);

  size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  /*
    Change the alphabet to the one specified by the string. The length of the string
    must be sigma, and each character can occur at most once.
  */
  bool assign(const std::string& alphabet_string);

private:
  int_vector<8>  m_char2comp, m_comp2char;
  int_vector<64> m_C;
  size_type      m_sigma;

  void copy(const Alphabet& v);

public:
  const int_vector<8>&  char2comp;
  const int_vector<8>&  comp2char;
  const int_vector<64>& C;
  const size_type&      sigma;
};  // class Alphabet

//------------------------------------------------------------------------------

/*
  This class uses an sd_vector to encode the cumulative sum of an array of integers.
  The array contains sum() items in size() elements. The array uses 0-based indexes.
*/

class CumulativeArray
{
public:
  typedef uint64_t size_type;

  CumulativeArray();
  CumulativeArray(const CumulativeArray& s);
  CumulativeArray(CumulativeArray&& s);
  ~CumulativeArray();

  /*
    The IntVector has to support operator[] that returns a non-const reference.
    The input is the original array, which gets overwritten by a cumulative array.
    This can be reversed by calling CumulativeArray::cumulativeToOriginal().
  */
  template<class IntVector>
  explicit CumulativeArray(IntVector& sequence)
  {
    this->m_size = sequence.size();

    for(size_type i = 1; i < this->size(); i++) { sequence[i] += sequence[i - 1] + 1; }
    this->v = sd_vector<>(sequence.begin(), sequence.end());

    util::init_support(this->rank, &(this->v));
    util::init_support(this->select_1, &(this->v));
    util::init_support(this->select_0, &(this->v));
  }

  template<class IntVector>
  static void
  cumulativeToOriginal(IntVector& sequence, size_type size)
  {
    for(size_type i = size - 1; i > 0; i--) { sequence[i] -= sequence[i - 1] + 1; }
  }

  void swap(CumulativeArray& v);
  CumulativeArray& operator=(const CumulativeArray& v);
  CumulativeArray& operator=(CumulativeArray&& v);

  size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  inline size_type size() const { return this->m_size; }

  // The sum of all elements.
  inline size_type sum() const { return this->v.size() - this->size(); }

  // The sum of the first k elements.
  inline size_type sum(size_type k) const
  {
    if(k == 0) { return 0; }
    if(k > this->size()) { k = this->size(); }

    return this->select_1(k) - k + 1;
  }

  inline size_type operator[](size_type i) const { return this->sum(i + 1) - this->sum(i); }

  // The inverse of sum(). Returns the element for item i.
  inline size_type inverse(size_type i) const
  {
    if(i >= this->sum()) { return this->size(); }

    return this->select_0(i + 1) - i;
  }

  // Is item i the last item in its element.
  inline bool isLast(size_type i) const
  {
    if(i >= this->sum()) { return false; }
    return this->v[this->select_0(i + 1) + 1];
  }

private:
  sd_vector<>                v;
  sd_vector<>::rank_1_type   rank;
  sd_vector<>::select_1_type select_1;
  sd_vector<>::select_0_type select_0;
  size_type                  m_size;  // Size of the original array.

  void copy(const CumulativeArray& v);
};  // class CumulativeArray

//------------------------------------------------------------------------------

/*
  This class implements a mapping between the LCS positions in two sequences.
*/
class LCS
{
public:
  typedef uint64_t size_type;

#ifdef USE_HYBRID_BITVECTORS
  typedef hyb_vector<> vector_type;
#else
  typedef rrr_vector<63> vector_type;
#endif

  LCS();
  LCS(const bit_vector& a, const bit_vector& b, size_type _lcs_size);
  LCS(const LCS& l);
  LCS(LCS&& l);
  ~LCS();

  void swap(LCS& l);
  LCS& operator=(const LCS& l);
  LCS& operator=(LCS&& l);

  uint64_t serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  inline uint64_t size() const { return this->lcs_size; }
  inline uint64_t ref_size() const { return this->ref.size(); }
  inline uint64_t seq_size() const { return this->seq.size(); }

  /*
    To find how many LCS bits are in ref before the 0-based position i, use ref_rank(i).
    To find the 1-based LCS bit i in seq, use seq_select(i).
  */

  vector_type                ref;
  vector_type::rank_1_type   ref_rank;
#ifdef USE_HYBRID_BITVECTORS
  inline uint64_t ref_select(uint64_t i) const { return this->select(this->ref, this->ref_rank, i); }
#else
  vector_type::select_1_type  ref_select;
#endif

  vector_type                seq;
  vector_type::rank_1_type   seq_rank;
#ifdef USE_HYBRID_BITVECTORS
  inline uint64_t seq_select(uint64_t i) const { return this->select(this->seq, this->seq_rank, i); }
#else
  vector_type::select_1_type  seq_select;
#endif

  size_type lcs_size;

private:
  void copy(const LCS& l);
  void set_vectors();

#ifdef USE_HYBRID_BITVECTORS
  uint64_t select(const vector_type& vec, const vector_type::rank_1_type& rank, uint64_t i) const;
#endif
};  // class LCS

//------------------------------------------------------------------------------

/*
  This class stores an integer array in two parts. Small values less than 255 are stored in
  an int_vector<8>, while large values are marked with 255 and stored separately in an
  int_vector<0>.
*/
class SLArray
{
public:
  typedef uint64_t size_type;
  typedef uint64_t value_type;

  const static size_type  BLOCK_SIZE  = 64;
  const static value_type LARGE_VALUE = 255;

  SLArray();
  SLArray(const SLArray& s);
  SLArray(SLArray&& s);
  ~SLArray();

  template<class IntVector>
  explicit SLArray(const IntVector& source)
  {
    util::assign(this->small, int_vector<8>(source.size(), 0));

    size_type  large_values = 0;
    value_type max_large = 0;
    for(size_type i = 0; i < this->size(); i++)
    {
      if(source[i] < LARGE_VALUE) { this->small[i] = source[i]; }
      else
      {
        this->small[i] = LARGE_VALUE;
        large_values++; max_large = std::max(max_large, (value_type)(source[i]));
      }
    }

    if(large_values > 0)
    {
      size_type blocks = (this->size() + BLOCK_SIZE - 1) / BLOCK_SIZE;
      util::assign(this->large, int_vector<0>(large_values, 0, bitlength(max_large)));
      util::assign(this->samples, int_vector<0>(blocks, 0, bitlength(large_values)));
      for(size_type block = 0, pos = 0, large_pos = 0; block < blocks; block++)
      {
        this->samples[block] = large_pos;
        size_type limit = std::min(this->size(), (block + 1) * BLOCK_SIZE);
        while(pos < limit)
        {
          if(source[pos] >= LARGE_VALUE)
          {
            this->large[large_pos] = source[pos]; large_pos++;
          }
          pos++;
        }
      }
    }
  }

  void swap(SLArray& s);
  SLArray& operator=(const SLArray& s);
  SLArray& operator=(SLArray&& s);

  uint64_t serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  inline size_type size() const { return this->small.size(); }
  inline size_type largeValues() const { return this->large.size(); }

  inline value_type operator[] (size_type i) const
  {
    if(this->small[i] < LARGE_VALUE) { return this->small[i]; }
    else { return this->large[this->large_rank(i)]; }
  }

  // Semiopen interval [from, to). No sanity checking.
  inline int_vector<64> extract(size_type from, size_type to) const
  {
    int_vector<64> result(to - from, 0);
    if(this->largeValues() == 0)
    {
      for(size_type i = from; i < to; i++) { result[i - from] = this->small[i]; }
    }
    else
    {
      for(size_type i = from, j = this->large_rank(from); i < to; i++)
      {
        if(this->small[i] < LARGE_VALUE) { result[i - from] = this->small[i]; }
        else { result[i - from] = this->large[j]; j++; }
      }
    }
    return result;
  }

  inline size_type large_rank(size_type i) const
  {
    size_type res = this->samples[i / BLOCK_SIZE];
    for(size_type j = i - i % BLOCK_SIZE; j < i; j++)
    {
      if(this->small[j] >= LARGE_VALUE) { res++; }
    }
    return res;
  }

  int_vector<8> small;
  int_vector<0> large;
  int_vector<0> samples;

private:
  void copy(const SLArray& s);
};  // class SLArray

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_SUPPORT_H
