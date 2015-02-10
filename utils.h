#ifndef _RELATIVE_FM_UTILS_H
#define _RELATIVE_FM_UTILS_H


#include <fstream>
#include <iostream>

#include <sdsl/wavelet_trees.hpp>


using namespace sdsl;

//------------------------------------------------------------------------------

extern const std::string BWT_EXTENSION;               // .bwt
extern const std::string NATIVE_BWT_EXTENSION;        // .cbwt
extern const std::string ALPHA_EXTENSION;             // .alpha
extern const std::string SIMPLE_FM_DEFAULT_ALPHABET;  // \0ACGNT
extern const std::string ROPEBWT2_ALPHABET;           // \0ACGTN

//------------------------------------------------------------------------------

typedef std::pair<uint64_t, uint64_t> range_type;

enum LoadMode { mode_plain, mode_native, mode_ropebwt2 };

template<class A, class B>
std::ostream& operator<<(std::ostream& stream, const std::pair<A, B>& data)
{
  return stream << "(" << data.first << ", " << data.second << ")";
}

inline uint64_t
length(range_type range)
{
  return range.second + 1 - range.first;
}

// An empty range is usually encoded as (x, y), where x > y.
// Use length for checking the emptiness of ranges of type (0, (uint64_t)-1).
inline bool
isEmpty(range_type range)
{
  return (range.first > range.second);
}

inline uint64_t
bitlength(uint64_t val)
{
  return bits::hi(val) + 1;
}

inline uint mapToUint(int val)
{
  if(val >= 0) { return 2 * val; }
  else         { return 2 * (-val) - 1; }
}

//------------------------------------------------------------------------------

template<class A>
void
readInteger(std::ifstream& input, A a)
{
  input.read((char*)&a, sizeof(A));
}

template<class A>
void
writeInteger(std::ofstream& output, A a)
{
  output.write((char*)&a, sizeof(A));
}

//------------------------------------------------------------------------------

inline double
inMegabytes(uint64_t bytes)
{
  return bytes / 1048576.0;
}

inline double
inBPC(uint64_t bytes, uint64_t size)
{
  return (8.0 * bytes) / size;
}

void printSize(const std::string& header, uint64_t bytes, uint64_t data_size, uint indent = 18);
void printTime(const std::string& header, uint64_t found, uint64_t matches, uint64_t bytes, double seconds, uint indent = 18);

//------------------------------------------------------------------------------

double readTimer();
uint64_t memoryUsage(); // Peak memory usage in bytes.

// Returns the total length of the rows, excluding line ends.
uint64_t readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows);

//------------------------------------------------------------------------------

// This function has been specialized for RLSequence.
template<class ByteVector>
void
characterCounts(const ByteVector& sequence, uint64_t size, int_vector<64>& counts)
{
  for(uint64_t c = 0; c < counts.size(); c++) { counts[c] = 0; }
  for(uint64_t i = 0; i < size; i++) { counts[sequence[i]]++; }
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
  Alphabet(const ByteVector& sequence, size_type size) :
    char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma)
  {
    util::assign(this->m_char2comp, int_vector<8>(MAX_SIGMA, 0));
    util::assign(this->m_comp2char, int_vector<8>(MAX_SIGMA, 0));
    util::assign(this->m_C, int_vector<64>(MAX_SIGMA + 1, 0));
    this->m_sigma = 0;
    if(size == 0) { return; }

    // Step 1: Count the occurrences and temporarily store them in m_C.
    characterCounts(sequence, size, this->m_C);

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

typedef wt_huff<bit_vector, rank_support_v5<> > bwt_type;

/*
  All the following functions assume that c is a real character, not a comp value in a
  contiguous alphabet.
*/

inline uint64_t
cumulative(const Alphabet& alpha, uint8_t c)
{
  return alpha.C[alpha.char2comp[c]];
}

inline bool
hasChar(const Alphabet& alpha, uint8_t c)
{
  return (c == 0 || alpha.char2comp[c] > 0);
}

template<class RankStructure>
inline range_type
charRange(const RankStructure& bwt, const Alphabet& alpha, uint8_t c)
{
  uint64_t comp = alpha.char2comp[c];
  if(comp < alpha.sigma) { return range_type(alpha.C[comp], alpha.C[comp + 1] - 1); }
  else { return range_type(alpha.C[comp], bwt.size() - 1); }
}

template<class RankStructure>
inline range_type
LF(const RankStructure& bwt, const Alphabet& alpha, range_type rng, uint8_t c)
{
  c = alpha.char2comp[c];
  uint64_t begin = alpha.C[c];
  return range_type(begin + bwt.rank(rng.first, c), begin + bwt.rank(rng.second + 1, c) - 1);
}

//------------------------------------------------------------------------------

/*
  Some SDSL extensions.
*/

template<class element>
uint64_t
write_vector(const std::vector<element>& vec, std::ostream& out, structure_tree_node* v, std::string name)
{
  structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(vec));
  uint64_t written_bytes = 0;
  written_bytes += write_member(vec.size(), out, child, "size");
  out.write((char*)(vec.data()), vec.size() * sizeof(element));
  written_bytes += vec.size() * sizeof(element);
  structure_tree::add_size(v, written_bytes);
  return written_bytes;
}

template<class element>
void
read_vector(std::vector<element>& vec, std::istream& in)
{
  uint64_t size = 0;
  read_member(size, in);
  vec.clear();
  std::vector<element> temp(size);
  in.read((char*)(temp.data()), vec.size() * sizeof(element));
  vec.swap(temp);
}

/*
  Generic in-memory construction from int_vector_buffer<8> and size. Not very space-efficient, as it
  duplicates the data.
*/
template<class Type>
void
directConstruct(Type& structure, const int_vector<8>& data)
{
  std::string ramfile = ram_file_name(util::to_string(&structure));
  store_to_file(data, ramfile);
  {
    int_vector_buffer<8> buffer(ramfile); // Must remove the buffer before removing the ramfile.
    Type temp(buffer, data.size());
    structure.swap(temp);
  }
  ram_fs::remove(ramfile);
}

//------------------------------------------------------------------------------

/*
  Determines the number of runs and the gap0, gap1, and run measures in bits.
*/
template<class VectorType>
void
countRuns(const VectorType& vec, uint64_t& runs, uint64_t& gap0, uint64_t& gap1, uint64_t& run, uint64_t& delta)
{
  bool in_run = false;
  uint64_t prev0 = 0, prev1 = 0, prevr = 0, onebits = 0;
  for(uint64_t i = 0; i < vec.size(); i++)
  {
    if(vec[i] == 1)
    {
      gap1 += bitlength(i + 1 - prev1);
      if(!in_run)
      {
        run += bitlength(i + 1 - prevr);
        runs++; prevr = i + 1;
      }
      in_run = true; prev1 = i + 1;
      onebits++;
    }
    else
    {
      gap0 += bitlength(i + 1 - prev0);
      if(in_run)
      {
        run += bitlength(i + 1 - prevr);
        prevr = i + 1;
      }
      in_run = false; prev0 = i + 1;
    }
  }

  int_vector<64> buffer(onebits);
  for(uint64_t i = 0, j = 0; i < vec.size(); i++)
  {
    if(vec[i] == 1) { buffer[j] = i; j++; }
  }
  enc_vector<> delta_vec(buffer);
  delta = size_in_bytes(delta_vec);
}

//------------------------------------------------------------------------------

/*
  This class uses an sd_vector to encode the cumulative sum of an array of integers.
  The array contains sum() items in size() elements.
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
  CumulativeArray(IntVector& sequence, size_type _size)
  {
    this->m_size = _size;

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

    return this->rank(this->select_0(i + 1));
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

#endif // _RELATIVE_FM_UTILS_H
