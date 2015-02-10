#ifndef _RELATIVE_FM_SIMPLE_FM_H
#define _RELATIVE_FM_SIMPLE_FM_H


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
    - Basic operations: swap(), size(), load()
    - With the structure as a parameter: serialize(bwt, .., ..)
*/
template<class RankStructure = bwt_type>
class SimpleFM
{
public:
  typedef RankStructure seq_type;

  // This constructor has been specialized for SimpleFM<RLSequence>.
  explicit SimpleFM(const std::string& base_name, LoadMode mode = mode_plain)
  {
    if(mode == mode_native)
    {
      std::string filename = base_name + NATIVE_BWT_EXTENSION;
      std::ifstream in(filename.c_str(), std::ios_base::binary);
      if(!in)
      {
        std::cerr << "SimpleFM::SimpleFM(): Cannot open BWT file " << filename << " (native format)" << std::endl;
        return;
      }
      this->bwt.load(in); in.close();
    }
    else if(mode == mode_ropebwt2)
    {
      std::cerr << "SimpleFM::SimpleFM(): Invalid sequence type for mode_ropebwt2!" << std::endl;
      return;
    }
    else
    {
      int_vector_buffer<8> buffer(base_name + BWT_EXTENSION);
      RankStructure temp(buffer, buffer.size());
      this->bwt.swap(temp);
    }
    this->loadAlphabet(base_name);
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

  // This function has been specialized for SimpleFM<RLSequence>.
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
  Alphabet alpha;

private:

  void loadAlphabet(const std::string& base_name)
  {
    std::string filename = base_name + ALPHA_EXTENSION;
    std::ifstream in(filename.c_str(), std::ios_base::binary);
    if(in) { this->alpha.load(in); in.close(); }
    else  // Try the default alphabet.
    {
      Alphabet temp(this->bwt, this->size());
      if(temp.sigma != SIMPLE_FM_DEFAULT_ALPHABET.length())
      {
        std::cerr << "SimpleFM()::SimpleFM(): Alphabet file " << filename
                  << " does not exist and the default alphabet cannot be used" << std::endl;
        std::cerr << "SimpleFM()::SimpleFM(): BWT alphabet size is " << temp.sigma
                  << ", while the default is " << SIMPLE_FM_DEFAULT_ALPHABET.length() << std::endl;
        return;
      }
      if(temp.assign(SIMPLE_FM_DEFAULT_ALPHABET)) { this->alpha = temp; }
    }
  }

  SimpleFM();
  SimpleFM(const SimpleFM&);
  SimpleFM(SimpleFM&&);
  SimpleFM& operator=(const SimpleFM&);
  SimpleFM& operator=(SimpleFM&&);
};

//------------------------------------------------------------------------------

#endif // _RELATIVE_FM_SIMPLE_FM_H
