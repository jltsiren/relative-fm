#ifndef _RELATIVE_FM_SIMPLE_FM_H
#define _RELATIVE_FM_SIMPLE_FM_H

#include <fstream>
#include <iostream>

#include <sdsl/csa_alphabet_strategy.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "utils.h"

namespace relative
{

//------------------------------------------------------------------------------

/*
  A simple FM-index for byte alphabet.

  The storage format is an int_vector<8> containing the BWT (with a contiguous alphabet)
  in base_name.bwt, and a byte_alphabet mapping the original alphabet to the contiguous
  alphabet in base_name.alpha.

  RankStructure must support the following operations:
    - Queries: [i], rank(i, c), inverse_select(i)
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
    this->sample_rate = 0;

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
      std::cerr << "SimpleFM::SimpleFM(): Invalid sequence type for mode_ropebwt2" << std::endl;
      return;
    }
    else
    {
      int_vector_buffer<8> buffer(base_name + BWT_EXTENSION);
      RankStructure temp(buffer, buffer.size());
      this->bwt.swap(temp);
    }

    this->loadAlphabet(base_name);
    this->loadSamples(base_name);
  }

  ~SimpleFM()
  {
  }

//------------------------------------------------------------------------------

  inline uint64_t size() const { return this->bwt.size(); }
  inline uint64_t sequences() const { return this->alpha.C[1]; }

  uint64_t reportSize(bool print = false) const
  {
    uint64_t bwt_bytes = size_in_bytes(this->bwt);
    uint64_t sample_bytes = sizeof(this->sample_rate) + size_in_bytes(this->samples);
    uint64_t bytes = bwt_bytes + size_in_bytes(this->alpha) + sample_bytes;

    if(print)
    {
#ifdef VERBOSE_OUTPUT
      printSize("BWT", bwt_bytes, this->size());
      if(this->sample_rate > 0) { printSize("Samples", sample_bytes, this->size()); }
#endif
      printSize("Simple FM", bytes, this->size());
      std::cout << std::endl;
    }

    return bytes;
  }

  void writeTo(const std::string& base_name, bool native_format = false) const
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

    if(native_format)
    {
      std::string filename = base_name + NATIVE_BWT_EXTENSION;
      std::ofstream output(filename.c_str(), std::ios_base::binary);
      if(!output)
      {
        std::cerr << "SimpleFM::writeTo(): Cannot open native BWT file " << filename << std::endl;
        return;
      }
      this->bwt.serialize(output); output.close();
    }
    else
    {
      std::string filename = base_name + BWT_EXTENSION;
      int_vector_buffer<8> output(filename, std::ios::out);
      int_vector<8> buffer(MEGABYTE, 0);
      for(uint64_t i = 0; i < this->size(); i += MEGABYTE)
      {
        range_type range(i, std::min(i + MEGABYTE, this->size()) - 1);
        this->extractBWT(range, buffer);
        for(uint64_t j = range.first; j <= range.second; j++) { output[j] = buffer[j - i]; }
      }
      output.close();
    }

    if(this->sample_rate > 0)
    {
      std::string filename = base_name + SAMPLE_EXTENSION;
      std::ofstream output(filename.c_str(), std::ios_base::binary);
      if(!output)
      {
        std::cerr << "SimpleFM::writeTo(): Cannot open SA sample file " << filename << std::endl;
        return;
      }
      write_member(this->sample_rate, output); this->samples.serialize(output); output.close();
    }
  }

//------------------------------------------------------------------------------

  // Use length(res) == 0 to check whether the range is empty.
  inline range_type LF(range_type range, uint8_t c) const
  {
    if(!hasChar(this->alpha, c)) { return range_type(1, 0); }
    return relative::LF(this->bwt, this->alpha, range, c);
  }

  // Returns (LF(i), BWT[i]) = (LF(i), char(LF(i))); the character is a comp value.
  inline range_type LF(uint64_t i) const
  {
    auto temp = this->bwt.inverse_select(i);
    return range_type(temp.first + this->alpha.C[temp.second], temp.second);
  }

  template<class Iter> range_type find(Iter begin, Iter end) const
  {
    range_type res(0, this->size() - 1);
    while(begin != end)
    {
      --end;
      res = this->LF(res, *end);
      if(length(res) == 0) { return range_type(1, 0); } // isEmpty() does not work with (0, -1).
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

  bool supportsLocate(bool print = false) const
  {
    if(this->sample_rate == 0)
    {
      if(print)
      {
        std::cerr << "SimpleFM::supportsLocate(): The index does not contain SA samples." << std::endl;
      }
      return false;
    }
    if(this->sequences() > 1)
    {
      if(print)
      {
        std::cerr << "SimpleFM::supportsLocate(): The index contains more than one sequence." << std::endl;
      }
      return false;
    }
    return true;
  }

  // Call supportsLocate() first.
  uint64_t locate(uint64_t i) const
  {
    if(i >= this->size()) { return 0; }

    uint64_t steps = 0;
    while(i % this->sample_rate != 0)
    {
      range_type res = this->LF(i);
      if(res.second == 0) { return steps; }
      i = res.first; steps++;
    }

    return this->samples[i / this->sample_rate] + steps;
  }

//------------------------------------------------------------------------------

  RankStructure bwt;
  Alphabet      alpha;
  uint64_t      sample_rate;
  int_vector<0> samples;

//------------------------------------------------------------------------------

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
        std::cerr << "SimpleFM()::loadAlphabet(): Alphabet file " << filename
                  << " does not exist and the default alphabet cannot be used" << std::endl;
        std::cerr << "SimpleFM()::loadAlphabet(): BWT alphabet size is " << temp.sigma
                  << ", while the default is " << SIMPLE_FM_DEFAULT_ALPHABET.length() << std::endl;
        return;
      }
      if(temp.assign(SIMPLE_FM_DEFAULT_ALPHABET)) { this->alpha = temp; }
    }
  }

  void loadSamples(const std::string& base_name)
  {
    std::string filename = base_name + SAMPLE_EXTENSION;
    std::ifstream in(filename.c_str(), std::ios_base::binary);
    if(in)
    {
      read_member(this->sample_rate, in); this->samples.load(in); in.close();
    }
  }

  SimpleFM();
  SimpleFM(const SimpleFM&);
  SimpleFM(SimpleFM&&);
  SimpleFM& operator=(const SimpleFM&);
  SimpleFM& operator=(SimpleFM&&);
};

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_SIMPLE_FM_H
