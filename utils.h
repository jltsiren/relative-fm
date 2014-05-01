#ifndef _RELATIVE_FM_UTILS_H
#define _RELATIVE_FM_UTILS_H


#include <fstream>
#include <iostream>

#include <sdsl/csa_alphabet_strategy.hpp>
#include <sdsl/wavelet_trees.hpp>


using namespace sdsl;

//------------------------------------------------------------------------------

typedef std::pair<uint64_t, uint64_t> range_type;

template<class A, class B>
std::ostream& operator<<(std::ostream& stream, const std::pair<A,B>& data)
{
  return stream << "(" << data.first << ", " << data.second << ")";
}

inline bool
isEmpty(range_type range)
{
  return (range.first > range.second);
}

inline uint64_t
length(range_type range)
{
  return range.second + 1 - range.first;
}

//------------------------------------------------------------------------------

template<class A>
void
readInteger(std::ifstream& input, A& a)
{
  input.read((char*)&a, sizeof(A));
}

template<class A>
void
writeInteger(std::ofstream& output, A& a)
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

void printSize(const std::string& header, uint64_t bytes, uint64_t data_size, uint indent = 16);
void printTime(const std::string& header, uint64_t found, uint64_t matches, uint64_t bytes, double seconds, uint indent = 16);

//------------------------------------------------------------------------------

double readTimer();
uint64_t memoryUsage(); // Peak memory usage in bytes.

//------------------------------------------------------------------------------

// Returns the total length of the rows, excluding line ends.
uint64_t readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows);

//------------------------------------------------------------------------------

typedef wt_huff<bit_vector, rank_support_v5<> > bwt_type;
typedef byte_alphabet alphabet_type;

inline uint64_t
cumulative(const alphabet_type& alpha, uint8_t c)
{
  return alpha.C[alpha.char2comp[c]];
}

inline bool
hasChar(const alphabet_type& alpha, uint8_t c)
{
  return (c == 0 || alpha.char2comp[c] > 0);
}

inline range_type
charRange(const bwt_type& bwt, const alphabet_type& alpha, uint8_t c)
{
  uint16_t temp = alpha.char2comp[c];
  if(temp < alpha.sigma) { return range_type(alpha.C[temp], alpha.C[temp + 1] - 1); }
  else { return range_type(alpha.C[temp], bwt.size() - 1); }
}

inline range_type
LF(const bwt_type& bwt, const alphabet_type& alpha, range_type rng, uint8_t c)
{
  uint64_t begin = cumulative(alpha, c);
  return range_type(begin + bwt.rank(rng.first, c), begin + bwt.rank(rng.second + 1, c) - 1);
}

//------------------------------------------------------------------------------

#endif // _RELATIVE_FM_UTILS_H
