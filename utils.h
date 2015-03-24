#ifndef _RELATIVE_FM_UTILS_H
#define _RELATIVE_FM_UTILS_H

#include <algorithm>
#include <fstream>
#include <iostream>

#include <sdsl/wavelet_trees.hpp>

/*
  Critical sections:

  alignbwts   Function alignBWTs.
  stderr      Writing to std::cerr within parallel sections.
  stdout      Writing to std::cout within parallel sections.
*/
#include <omp.h>

using namespace sdsl;

namespace relative
{

//------------------------------------------------------------------------------

extern const std::string BWT_EXTENSION;               // .bwt
extern const std::string NATIVE_BWT_EXTENSION;        // .cbwt
extern const std::string ALPHA_EXTENSION;             // .alpha
extern const std::string SAMPLE_EXTENSION;            // .samples
extern const std::string SIMPLE_FM_DEFAULT_ALPHABET;  // \0ACGNT
extern const std::string ROPEBWT2_ALPHABET;           // \0ACGTN

const uint64_t MEGABYTE = 1048576;
const double MEGABYTE_DOUBLE = 1048576.0;

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

const uint64_t FNV_OFFSET_BASIS = 0xcbf29ce484222325UL;
const uint64_t FNV_PRIME        = 0x100000001b3UL;

inline uint64_t fnv1a_hash(uint8_t c, uint64_t seed)
{
  return (seed ^ c) * FNV_PRIME;
}

inline uint64_t fnv1a_hash(uint64_t val, uint64_t seed)
{
  uint8_t* chars = (uint8_t*)&val;
  for(uint64_t i = 0; i < 8; i++) { seed = fnv1a_hash(chars[i], seed); }
  return seed;
}

//------------------------------------------------------------------------------

inline double
inMegabytes(uint64_t bytes)
{
  return bytes / MEGABYTE_DOUBLE;
}

inline double
inBPC(uint64_t bytes, uint64_t size)
{
  return (8.0 * bytes) / size;
}

void printSize(const std::string& header, uint64_t bytes, uint64_t data_size, uint64_t indent = 18);
void printTime(const std::string& header, uint64_t found, uint64_t matches, uint64_t bytes, double seconds, bool occs, uint64_t indent = 18);

//------------------------------------------------------------------------------

double readTimer();
uint64_t memoryUsage(); // Peak memory usage in bytes.

// Returns the total length of the rows, excluding line ends.
uint64_t readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows);

template<class Iterator, class Comparator>
void
parallelQuickSort(Iterator first, Iterator last, const Comparator& comp)
{
#ifdef _GLIBCXX_PARALLEL
  std::sort(first, last, comp, __gnu_parallel::balanced_quicksort_tag());
#else
  std::sort(first, last, comp);
#endif
}

template<class Iterator>
void
parallelQuickSort(Iterator first, Iterator last)
{
#ifdef _GLIBCXX_PARALLEL
  std::sort(first, last, __gnu_parallel::balanced_quicksort_tag());
#else
  std::sort(first, last);
#endif
}

template<class Iterator, class Comparator>
void
parallelMergeSort(Iterator first, Iterator last, const Comparator& comp)
{
#ifdef _GLIBCXX_PARALLEL
  std::sort(first, last, comp, __gnu_parallel::multiway_mergesort_tag());
#else
  std::sort(first, last, comp);
#endif
}

template<class Iterator>
void
parallelMergeSort(Iterator first, Iterator last)
{
#ifdef _GLIBCXX_PARALLEL
  std::sort(first, last, __gnu_parallel::multiway_mergesort_tag());
#else
  std::sort(first, last);
#endif
}

template<class Iterator, class Comparator>
void
sequentialSort(Iterator first, Iterator last, const Comparator& comp)
{
#ifdef _GLIBCXX_PARALLEL
  std::sort(first, last, comp, __gnu_parallel::sequential_tag());
#else
  std::sort(first, last, comp);
#endif
}

template<class Iterator>
void
sequentialSort(Iterator first, Iterator last)
{
#ifdef _GLIBCXX_PARALLEL
  std::sort(first, last, __gnu_parallel::sequential_tag());
#else
  std::sort(first, last);
#endif
}

//------------------------------------------------------------------------------

typedef wt_huff<bit_vector, rank_support_v5<> > bwt_type;

/*
  All the following functions assume that c is a real character, not a comp value in a
  contiguous alphabet.
*/

template<class AlphabetType>
inline uint64_t
cumulative(const AlphabetType& alpha, uint8_t c)
{
  return alpha.C[alpha.char2comp[c]];
}

template<class AlphabetType>
inline bool
hasChar(const AlphabetType& alpha, uint8_t c)
{
  return (c == 0 || alpha.char2comp[c] > 0);
}

template<class RankStructure, class AlphabetType>
inline range_type
charRange(const RankStructure& bwt, const AlphabetType& alpha, uint8_t c)
{
  uint64_t comp = alpha.char2comp[c];
  if(comp < alpha.sigma) { return range_type(alpha.C[comp], alpha.C[comp + 1] - 1); }
  else { return range_type(alpha.C[comp], bwt.size() - 1); }
}

template<class RankStructure, class AlphabetType>
inline range_type
LF(const RankStructure& bwt, const AlphabetType& alpha, range_type rng, uint8_t c)
{
  c = alpha.char2comp[c];
  uint64_t begin = alpha.C[c];
  return range_type(begin + bwt.rank(rng.first, c), begin + bwt.rank(rng.second + 1, c) - 1);
}

/*
  Returns SA[i]. Call index.supportsLocate() first.
*/
template<class Index>
uint64_t
locate(const Index& index, uint64_t i)
{
  if(i >= index.size()) { return 0; }

  uint64_t steps = 0;
  while(i % index.sa_sample_rate != 0)
  {
    range_type res = index.LF(i);
    if(res.second == 0) { return steps; }
    i = res.first; steps++;
  }

  return index.sa_samples[i / index.sa_sample_rate] + steps;
}

/*
  Extracts text[range]. Call index.supportsExtract() first.
*/
template<class Index>
std::string
extract(const Index& index, range_type range)
{
  if(range.second >= index.size()) { range.second = index.size() - 1; }
  if(isEmpty(range)) { return std::string(); }

  uint64_t bwt_pos = index.inverse(range.second), text_pos = range.second, c = 0;
  while(index.alpha.C[c + 1] <= bwt_pos) { c++; }

  // Extract the sequence.
  std::string result(length(range), 0);
  while(text_pos > range.first)
  {
    result[text_pos - range.first] = index.alpha.comp2char[c];
    range_type temp = index.LF(bwt_pos);
    bwt_pos = temp.first; c = temp.second; text_pos--;
  }
  result[0] = index.alpha.comp2char[c];

  return result;
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
  {
    std::vector<element> temp;
    vec.swap(temp);
  }
  {
    std::vector<element> temp(size);
    in.read((char*)(temp.data()), temp.size() * sizeof(element));
    vec.swap(temp);
  }
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

} // namespace relative

#endif // _RELATIVE_FM_UTILS_H
