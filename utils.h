/*
  Copyright (c) 2015, 2016, 2017 Genome Research Ltd.
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

//------------------------------------------------------------------------------

typedef std::uint64_t size_type;
typedef std::uint8_t  char_type;
typedef std::uint8_t  comp_type;
typedef std::uint8_t  byte_type;

const size_type BYTE_BITS = 8;
const size_type WORD_BITS = 64;

const size_type KILOBYTE     = 1024;
const size_type MILLION      = 1000000;
const size_type MEGABYTE     = KILOBYTE * KILOBYTE;
const size_type GIGABYTE     = KILOBYTE * MEGABYTE;

const double KILOBYTE_DOUBLE = 1024.0;
const double MILLION_DOUBLE  = 1000000.0;
const double MEGABYTE_DOUBLE = KILOBYTE_DOUBLE * KILOBYTE_DOUBLE;
const double GIGABYTE_DOUBLE = KILOBYTE_DOUBLE * MEGABYTE_DOUBLE;

using sdsl::bit_vector;

//------------------------------------------------------------------------------

extern const std::string BWT_EXTENSION;               // .bwt
extern const std::string NATIVE_BWT_EXTENSION;        // .cbwt
extern const std::string ALPHA_EXTENSION;             // .alpha
extern const std::string SAMPLE_EXTENSION;            // .samples
extern const std::string LCP_EXTENSION;               // .lcp
extern const std::string DLCP_EXTENSION;              // .dlcp
extern const std::string DLCP_INDEX_EXTENSION;        // .dlcp_index
extern const std::string SIMPLE_FM_DEFAULT_ALPHABET;  // \0ACGNT

extern const std::string RLCP_EXTENSION;              // .rlcp

//------------------------------------------------------------------------------

/*
  range_type stores a closed range [first, second]. Empty ranges are indicated by
  first > second. The emptiness check uses +1 to handle the common special case
  [0, -1].
*/

typedef std::pair<size_type, size_type> range_type;

struct Range
{
  inline static size_type length(range_type range)
  {
    return range.second + 1 - range.first;
  }

  inline static bool empty(range_type range)
  {
    return (range.first + 1 > range.second + 1);
  }

  inline static bool empty(size_type sp, size_type ep)
  {
    return (sp + 1 > ep + 1);
  }

  inline static size_type bound(size_type value, range_type bounds)
  {
    return bound(value, bounds.first, bounds.second);
  }

  inline static size_type bound(size_type value, size_type low, size_type high)
  {
    return std::max(std::min(value, high), low);
  }

  inline static range_type empty_range()
  {
    return range_type(1, 0);
  }
};

template<class A, class B>
std::ostream& operator<<(std::ostream& stream, const std::pair<A, B>& data)
{
  return stream << "(" << data.first << ", " << data.second << ")";
}

//------------------------------------------------------------------------------

template<class IntegerType>
inline size_type
bit_length(IntegerType val)
{
  return sdsl::bits::hi(val) + 1;
}

//------------------------------------------------------------------------------

const size_type FNV_OFFSET_BASIS = 0xcbf29ce484222325UL;
const size_type FNV_PRIME        = 0x100000001b3UL;

inline size_type fnv1a_hash(byte_type c, size_type seed)
{
  return (seed ^ c) * FNV_PRIME;
}

inline size_type fnv1a_hash(size_type val, size_type seed)
{
  byte_type* chars = (byte_type*)&val;
  for(size_type i = 0; i < 8; i++) { seed = fnv1a_hash(chars[i], seed); }
  return seed;
}

//------------------------------------------------------------------------------

inline double
inMegabytes(size_type bytes)
{
  return bytes / MEGABYTE_DOUBLE;
}

inline double
inGigabytes(size_type bytes)
{
  return bytes / GIGABYTE_DOUBLE;
}

inline double
inBPC(size_type bytes, size_type size)
{
  return (8.0 * bytes) / size;
}

inline double
inMicroseconds(double seconds)
{
  return seconds * MILLION_DOUBLE;
}

const size_type DEFAULT_INDENT = 18;

void printHeader(const std::string& header, size_type indent = DEFAULT_INDENT);
void printSize(const std::string& header, size_type bytes, size_type data_size, size_type indent = DEFAULT_INDENT);
void printTime(const std::string& header, size_type found, size_type matches, size_type bytes, double seconds, bool occs, size_type indent = DEFAULT_INDENT);
void printTime(const std::string& header, size_type queries, double seconds, size_type indent = DEFAULT_INDENT);

//------------------------------------------------------------------------------

double readTimer();       // Seconds from an arbitrary time point.
size_type memoryUsage();  // Peak memory usage in bytes.

//------------------------------------------------------------------------------

// Returns the total length of the rows, excluding line ends.
size_type readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows);

//------------------------------------------------------------------------------

/*
  parallelQuickSort() uses less working space than parallelMergeSort(). Calling omp_set_nested(1)
  improves the speed of parallelQuickSort().
*/

template<class Iterator, class Comparator>
void
parallelQuickSort(Iterator first, Iterator last, const Comparator& comp)
{
#ifdef _GLIBCXX_PARALLEL
  int nested = omp_get_nested();
  omp_set_nested(1);
  std::sort(first, last, comp, __gnu_parallel::balanced_quicksort_tag());
  omp_set_nested(nested);
#else
  std::sort(first, last, comp);
#endif
}

template<class Iterator>
void
parallelQuickSort(Iterator first, Iterator last)
{
#ifdef _GLIBCXX_PARALLEL
  int nested = omp_get_nested();
  omp_set_nested(1);
  std::sort(first, last, __gnu_parallel::balanced_quicksort_tag());
  omp_set_nested(nested);
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

template<class Element>
void
removeDuplicates(std::vector<Element>& vec, bool parallel)
{
  if(parallel) { parallelQuickSort(vec.begin(), vec.end()); }
  else         { sequentialSort(vec.begin(), vec.end()); }
  vec.resize(std::unique(vec.begin(), vec.end()) - vec.begin());
}

//------------------------------------------------------------------------------

typedef sdsl::wt_huff<> bwt_type;

/*
  All the following functions assume that c is a real character, not a comp value in a
  contiguous alphabet.
*/

template<class AlphabetType>
inline size_type
cumulative(const AlphabetType& alpha, char_type c)
{
  return alpha.C[alpha.char2comp[c]];
}

template<class AlphabetType>
inline bool
hasChar(const AlphabetType& alpha, char_type c)
{
  return (c == 0 || alpha.char2comp[c] > 0);
}

template<class RankStructure, class AlphabetType>
inline range_type
charRange(const RankStructure& bwt, const AlphabetType& alpha, char_type c)
{
  size_type comp = alpha.char2comp[c];
  if(comp < alpha.sigma) { return range_type(alpha.C[comp], alpha.C[comp + 1] - 1); }
  else { return range_type(alpha.C[comp], bwt.size() - 1); }
}

template<class RankStructure, class AlphabetType>
inline range_type
LF(const RankStructure& bwt, const AlphabetType& alpha, range_type rng, char_type c)
{
  c = alpha.char2comp[c];
  size_type begin = alpha.C[c];
  return range_type(begin + bwt.rank(rng.first, c), begin + bwt.rank(rng.second + 1, c) - 1);
}

/*
  The comp value or the character value of T[SA[i]]. No sanity checking.
*/
template<class AlphabetType>
comp_type
findComp(const AlphabetType& alpha, size_type bwt_pos)
{
  size_type comp = 0;
  while(alpha.C[comp + 1] <= bwt_pos) { comp++; }
  return comp;
}

template<class AlphabetType>
char_type
findChar(const AlphabetType& alpha, size_type bwt_pos)
{
  return alpha.comp2char[findComp(alpha, bwt_pos)];
}

/*
  Returns SA[i]. Call index.supportsLocate() first.
*/
template<class Index>
size_type
locate(const Index& index, size_type i)
{
  if(i >= index.size()) { return index.size(); }

  size_type steps = 0;
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
  if(Range::empty(range)) { return std::string(); }

  size_type bwt_pos = index.inverse(range.second), text_pos = range.second;
  size_type c = findComp(index.alpha, bwt_pos);

  // Extract the sequence.
  std::string result(Range::length(range), 0);
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
size_type
Psi(const Index& index, size_type i, size_type k, size_type threshold)
{
  if(2 * k >= threshold) { return index.inverse(index.locate(i) + k); }
  else
  {
    for(size_type j = 0; j < k; j++) { i = index.Psi(i); }
    return i;
  }
}

//------------------------------------------------------------------------------

/*
  Generic in-memory construction from sdsl::int_vector_buffer<8> and size. Not very space-efficient, as it
  duplicates the data.
*/
template<class Type>
void
directConstruct(Type& structure, const sdsl::int_vector<8>& data)
{
  std::string ramfile = sdsl::ram_file_name(sdsl::util::to_string(&structure));
  sdsl::store_to_file(data, ramfile);
  {
    sdsl::int_vector_buffer<8> buffer(ramfile); // Must remove the buffer before removing the ramfile.
    Type temp(buffer, data.size());
    structure.swap(temp);
  }
  sdsl::ram_fs::remove(ramfile);
}

/*
  Extracts the given range from source, overwriting target.
*/
template<class VectorType>
void
extractBits(const VectorType& source, range_type range, sdsl::bit_vector& target)
{
  if(Range::empty(range) || range.second >= source.size()) { return; }

  sdsl::util::assign(target, sdsl::bit_vector(Range::length(range), 0));
  for(size_type i = 0; i < target.size(); i += 64)
  {
    size_type len = std::min((size_type)64, target.size() - i);
    target.set_int(i, source.get_int(range.first + i, len), len);
  }
}

//------------------------------------------------------------------------------

/*
  Determines the number of runs and the gap0, gap1, and run measures in bits.
*/
template<class VectorType>
void
countRuns(const VectorType& vec, size_type& runs, size_type& gap0, size_type& gap1, size_type& run, size_type& delta)
{
  bool in_run = false;
  size_type prev0 = 0, prev1 = 0, prevr = 0, onebits = 0;
  for(size_type i = 0; i < vec.size(); i++)
  {
    if(vec[i] == 1)
    {
      gap1 += bit_length(i + 1 - prev1);
      if(!in_run)
      {
        run += bit_length(i + 1 - prevr);
        runs++; prevr = i + 1;
      }
      in_run = true; prev1 = i + 1;
      onebits++;
    }
    else
    {
      gap0 += bit_length(i + 1 - prev0);
      if(in_run)
      {
        run += bit_length(i + 1 - prevr);
        prevr = i + 1;
      }
      in_run = false; prev0 = i + 1;
    }
  }

  sdsl::int_vector<64> buffer(onebits);
  for(size_type i = 0, j = 0; i < vec.size(); i++)
  {
    if(vec[i] == 1) { buffer[j] = i; j++; }
  }
  sdsl::enc_vector<> delta_vec(buffer);
  delta = sdsl::size_in_bytes(delta_vec);
}

//------------------------------------------------------------------------------

/*
  This encodes (val - prev) as a non-negative integer.
*/
struct DiffEncoder
{
  inline static size_type encode(size_type val, size_type prev)
  {
    if(val >= prev) { return 2 * (val - prev); }
    else            { return 2 * (prev - val) - 1; }
  }

  inline static size_type decode(size_type code, size_type prev)
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
  inline static size_type encode(size_type val, size_type prev)
  {
    if(val >= prev) { return 2 * (val - prev) + 1; }
    else            { return 2 * (prev - val); }
  }

  inline static size_type decode(size_type code, size_type prev)
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
differentialArray(const IntVector& source, sdsl::int_vector<0>& target, bool endmarker)
{
  sdsl::util::assign(target, sdsl::int_vector<0>(source.size() + endmarker, 0, target.width()));
  for(size_type i = 0, prev = 0; i < source.size(); i ++)
  {
    size_type curr = source[i];
    target[i] = Encoder::encode(curr, prev);
    prev = curr;
  }
}

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_UTILS_H
