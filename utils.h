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

void printSize(const std::string& header, uint64_t bytes, uint64_t data_size, uint indent = 16);
void printTime(const std::string& header, uint64_t found, uint64_t matches, uint64_t bytes, double seconds, uint indent = 16);

//------------------------------------------------------------------------------

double readTimer();
uint64_t memoryUsage(); // Peak memory usage in bytes.

// Returns the total length of the rows, excluding line ends.
uint64_t readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows);

//------------------------------------------------------------------------------

const std::string ALPHA_EXTENSION = ".alpha";
const std::string BWT_EXTENSION = ".bwt";

typedef wt_huff<bit_vector, rank_support_v5<> > bwt_type;
typedef byte_alphabet alphabet_type;

/*
  All the following functions assume that c is a real character, not a comp value in a
  contiguous alphabet.
*/

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

template<class RankStructure>
inline range_type
charRange(const RankStructure& bwt, const alphabet_type& alpha, uint8_t c)
{
  uint64_t comp = alpha.char2comp[c];
  if(comp < alpha.sigma) { return range_type(alpha.C[comp], alpha.C[comp + 1] - 1); }
  else { return range_type(alpha.C[comp], bwt.size() - 1); }
}

template<class RankStructure>
inline range_type
LF(const RankStructure& bwt, const alphabet_type& alpha, range_type rng, uint8_t c)
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
  int_vector_buffer<8> buffer(ramfile);
  Type temp(buffer, data.size());
  structure.swap(temp);
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

#endif // _RELATIVE_FM_UTILS_H
