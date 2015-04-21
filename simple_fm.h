#ifndef _RELATIVE_FM_SIMPLE_FM_H
#define _RELATIVE_FM_SIMPLE_FM_H

#include <fstream>
#include <iostream>

#include "support.h"

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
    this->sa_sample_rate = 0; this->isa_sample_rate = 0;

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
    else if(mode == mode_ropebwt)
    {
      std::cerr << "SimpleFM::SimpleFM(): Invalid sequence type for mode_ropebwt" << std::endl;
      return;
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
    uint64_t sa_bytes = size_in_bytes(this->sa_samples);
    uint64_t isa_bytes = size_in_bytes(this->isa_samples);
    uint64_t bytes = bwt_bytes + size_in_bytes(this->alpha) +
      sizeof(this->sa_sample_rate) + sa_bytes + sizeof(this->isa_sample_rate) + isa_bytes;

    if(print)
    {
#ifdef VERBOSE_OUTPUT
      printSize("BWT", bwt_bytes, this->size());
      if(this->sa_sample_rate > 0) { printSize("SA samples", sa_bytes, this->size()); }
      if(this->isa_sample_rate > 0) { printSize("ISA samples", isa_bytes, this->size()); }
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

    if(this->sa_sample_rate > 0)
    {
      std::string filename = base_name + SAMPLE_EXTENSION;
      std::ofstream output(filename.c_str(), std::ios_base::binary);
      if(!output)
      {
        std::cerr << "SimpleFM::writeTo(): Cannot open SA sample file " << filename << std::endl;
        return;
      }
      write_member(this->sa_sample_rate, output); this->sa_samples.serialize(output);
      write_member(this->isa_sample_rate, output); this->isa_samples.serialize(output);
      output.close();
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

  uint64_t Psi(uint64_t i) const
  {
    if(i >= this->size()) { return this->size(); }
    uint64_t comp = relative::findComp(this->alpha, i);
    return this->bwt.select(i + 1 - this->alpha.C[comp], comp);
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
    if(this->sa_sample_rate == 0)
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
  inline uint64_t locate(uint64_t i) const
  {
    return relative::locate(*this, i);
  }

  bool supportsExtract(bool print = false) const
  {
    if(this->isa_sample_rate == 0)
    {
      if(print)
      {
        std::cerr << "SimpleFM::supportsExtract(): The index does not contain ISA samples." << std::endl;
      }
      return false;
    }
    if(this->sequences() > 1)
    {
      if(print)
      {
        std::cerr << "SimpleFM::supportsExtract(): The index contains more than one sequence." << std::endl;
      }
      return false;
    }
    return true;
  }

  // Call supportsExtract() first.
  inline std::string extract(range_type range) const
  {
    return relative::extract(*this, range);
  }

  // Call supportsExtract() first.
  inline std::string extract(uint64_t from, uint64_t to) const
  {
    return relative::extract(*this, range_type(from, to));
  }

  // Returns ISA[i]. Call supportsExtract() first.
  inline uint64_t inverse(uint64_t i) const
  {
    uint64_t bwt_pos = 0, text_pos =
      ((i + this->isa_sample_rate - 1) / this->isa_sample_rate) * this->isa_sample_rate;
    if(text_pos >= this->size()) { text_pos = this->size() - 1; }
    else { bwt_pos = this->isa_samples[text_pos / this->isa_sample_rate]; }

    // Move backwards from the sample.
    while(text_pos > i)
    {
      bwt_pos = this->LF(bwt_pos).first; text_pos--;
    }

    return bwt_pos;
  }

//------------------------------------------------------------------------------

  RankStructure bwt;
  Alphabet      alpha;
  uint64_t      sa_sample_rate, isa_sample_rate;
  int_vector<0> sa_samples, isa_samples;

//------------------------------------------------------------------------------

private:

  void loadAlphabet(const std::string& base_name)
  {
    std::string filename = base_name + ALPHA_EXTENSION;
    std::ifstream in(filename.c_str(), std::ios_base::binary);
    if(in) { this->alpha.load(in); in.close(); }
    else  // Try the default alphabet.
    {
      Alphabet temp(this->bwt);
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
      read_member(this->sa_sample_rate, in); this->sa_samples.load(in);
      read_member(this->isa_sample_rate, in); this->isa_samples.load(in);
      in.close();
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
