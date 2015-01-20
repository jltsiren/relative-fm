#ifndef _RELATIVE_FM_SEQUENCE_H
#define _RELATIVE_FM_SEQUENCE_H


#include "relative_fm.h"

//------------------------------------------------------------------------------

class Sequence
{
public:
  const static uint64_t SAMPLE_RATE = 128;
  typedef uint64_t size_type;

  Sequence();
  Sequence(int_vector_buffer<8>& buffer, uint64_t _size, uint64_t _sigma = 0);  // Set sigma in advance to save memory.
  Sequence(const Sequence& s);
  Sequence(Sequence&& s);
  ~Sequence();

  void swap(Sequence& v);
  Sequence& operator=(const Sequence& v);
  Sequence& operator=(Sequence&& v);

  uint64_t serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  inline uint64_t size() const { return this->data.size(); }

  inline uint64_t rank(uint64_t i, uint8_t c) const
  {
    if(c >= this->sigma) { return 0; }
    if(i > this->size()) { i = this->size(); }

    uint64_t block = i / SAMPLE_RATE;
    uint64_t res = this->samples[block * this->sigma + c];
    for(uint64_t pos = block * SAMPLE_RATE; pos < i; pos++) { if(this->data[pos] == c) { res++; } }
    return res;
  }

  inline uint64_t operator[](uint64_t i) const { return this->data[i]; }

private:
  int_vector<0> data, samples;
  uint64_t      sigma;

  void copy(const Sequence& v);
  void buildRank();
};  // class Sequence

template<>
SimpleFM<Sequence>::SimpleFM(const std::string& base_name);

//------------------------------------------------------------------------------


#endif // _RELATIVE_FM_SEQUENCE_H
