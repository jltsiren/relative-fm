/*
  Copyright (c) 2015 Genome Research Ltd.
  Copyright (c) 2014 Jouni Siren

  Author: Jouni Siren <jouni.siren@iki.fi>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

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

namespace relative
{

using namespace sdsl;

//------------------------------------------------------------------------------

extern const std::string BWT_EXTENSION;               // .bwt
extern const std::string NATIVE_BWT_EXTENSION;        // .cbwt
extern const std::string ALPHA_EXTENSION;             // .alpha
extern const std::string SAMPLE_EXTENSION;            // .samples
extern const std::string LCP_EXTENSION;               // .lcp
extern const std::string DLCP_INDEX_EXTENSION;        // .dlcp_index
extern const std::string SIMPLE_FM_DEFAULT_ALPHABET;  // \0ACGNT
extern const std::string ROPEBWT_ALPHABET;            // \0ACGTN

const uint64_t MEGABYTE = 1048576;
const double MEGABYTE_DOUBLE = 1048576.0;

//------------------------------------------------------------------------------

typedef std::pair<uint64_t, uint64_t> range_type;

enum LoadMode { mode_plain, mode_native, mode_ropebwt, mode_ropebwt2 };

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
// The function works now for (0, -1) as well.
inline bool
isEmpty(range_type range)
{
  return (range.first + 1 > range.second + 1);
}

inline uint64_t
bitlength(uint64_t val)
{
  return bits::hi(val) + 1;
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

inline double
inMicroseconds(double seconds)
{
  return seconds * 1000000;
}

void printSize(const std::string& header, uint64_t bytes, uint64_t data_size, uint64_t indent = 18);
void printTime(const std::string& header, uint64_t found, uint64_t matches, uint64_t bytes, double seconds, bool occs, uint64_t indent = 18);
void printTime(const std::string& header, uint64_t queries, double seconds, uint64_t indent = 18);

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
  The comp value or the character value of T[SA[i]]. No sanity checking.
*/
template<class AlphabetType>
uint8_t
findComp(const AlphabetType& alpha, uint64_t bwt_pos)
{
  uint64_t comp = 0;
  while(alpha.C[comp + 1] <= bwt_pos) { comp++; }
  return comp;
}

template<class AlphabetType>
uint8_t
findChar(const AlphabetType& alpha, uint64_t bwt_pos)
{
  return alpha.comp2char[findComp(alpha, bwt_pos)];
}

/*
  Returns SA[i]. Call index.supportsLocate() first.
*/
template<class Index>
uint64_t
locate(const Index& index, uint64_t i)
{
  if(i >= index.size()) { return index.size(); }

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

  uint64_t bwt_pos = index.inverse(range.second), text_pos = range.second;
  uint64_t c = findComp(index.alpha, bwt_pos);

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

/*
  Iterates index.Psi() k times. Uses inverse() and locate() if 2 * k >= threshold.
  Return index.size() if i >= index.size() or SA[i]+k >= index.size().
*/
template<class Index>
uint64_t
Psi(const Index& index, uint64_t i, uint64_t k, uint64_t threshold)
{
  if(2 * k >= threshold) { return index.inverse(index.locate(i) + k); }
  else
  {
    for(uint64_t j = 0; j < k; j++) { i = index.Psi(i); }
    return i;
  }
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

/*
  Extracts the given range from source, overwriting target.
*/
template<class VectorType>
void
extractBits(const VectorType& source, range_type range, bit_vector& target)
{
  if(isEmpty(range) || range.second >= source.size()) { return; }

  util::assign(target, bit_vector(length(range), 0));
  for(uint64_t i = 0; i < target.size(); i += 64)
  {
    uint64_t len = std::min((uint64_t)64, target.size() - i);
    target.set_int(i, source.get_int(range.first + i, len), len);
  }
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
  This encodes (val - prev) as a non-negative integer.
*/
struct DiffEncoder
{
  inline static uint64_t encode(uint64_t val, uint64_t prev)
  {
    if(val >= prev) { return 2 * (val - prev); }
    else            { return 2 * (prev - val) - 1; }
  }

  inline static uint64_t decode(uint64_t code, uint64_t prev)
  {
    if(code & 1) { return prev - (code + 1) / 2; }
    else         { return code / 2 + prev; }
  }
};

/*
  This encodes (val - prev) as a positive integer. SA construction does not work with
  character value 0.
*/
struct DiffEncoderNZ
{
  inline static uint64_t encode(uint64_t val, uint64_t prev)
  {
    if(val >= prev) { return 2 * (val - prev) + 1; }
    else            { return 2 * (prev - val); }
  }

  inline static uint64_t decode(uint64_t code, uint64_t prev)
  {
    if(code & 1) { return (code - 1) / 2 + prev; }
    else         { return prev - code / 2; }
  }
};

/*
  Use endmarker == true if you want to add an endmarker with value 0. In that case,
  you should probably use DiffEncoderNZ instead of DiffEncoder. Set target width in
  advance to save memory.
*/
template<class IntVector, class Encoder = DiffEncoder>
void
differentialArray(const IntVector& source, int_vector<0>& target, bool endmarker)
{
  util::assign(target, int_vector<0>(source.size() + endmarker, 0, target.width()));
  for(uint64_t i = 0, prev = 0; i < source.size(); i ++)
  {
    uint64_t curr = source[i];
    target[i] = Encoder::encode(curr, prev);
    prev = curr;
  }
}

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_UTILS_H
