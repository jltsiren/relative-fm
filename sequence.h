#ifndef _RELATIVE_FM_SEQUENCE_H
#define _RELATIVE_FM_SEQUENCE_H


#include "simple_fm.h"

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

/*
  A basic run-length encoded sequence for DNA alphabet (\0, A, C, G, N, T)
  packed as values 0-5.
*/
class RLSequence
{
public:
  const static uint64_t SAMPLE_RATE = 64;
  const static uint64_t SIGMA = 6;
  const static uint64_t MAX_RUN = 256 / SIGMA;
  typedef uint64_t size_type;

  RLSequence();
  RLSequence(int_vector_buffer<8>& buffer, uint64_t _size);
  RLSequence(const RLSequence& s);
  RLSequence(RLSequence&& s);
  ~RLSequence();

  void swap(RLSequence& v);
  RLSequence& operator=(const RLSequence& v);
  RLSequence& operator=(RLSequence&& v);

  uint64_t serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  inline uint64_t size() const { return this->block_boundaries.size(); }
  inline uint64_t runs() const { return this->data.size(); }

  inline uint64_t rank(uint64_t i, uint8_t c) const
  {
    if(c >= SIGMA) { return 0; }
    if(i > this->size()) { i = this->size(); }

    uint64_t block = this->block_rank(i);
    uint64_t res = this->samples[block * SIGMA + c];
    uint64_t rle_pos = block * SAMPLE_RATE;
    uint64_t seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);
    while(seq_pos < i)
    {
      seq_pos += this->data[rle_pos] / SIGMA + 1; // The starting position of the next run.
      if(this->data[rle_pos] % SIGMA == c)
      {
        res += this->data[rle_pos] / SIGMA + 1; // Number of c's before the next run.
        if(seq_pos > i) { res -= seq_pos - i; }
      }
      rle_pos++;
    }

    return res;
  }

  inline uint64_t operator[](uint64_t i) const
  {
    if(i >= this->size()) { return 0; }

    uint64_t block = this->block_rank(i);
    uint64_t rle_pos = block * SAMPLE_RATE;
    uint64_t seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);
    while(true)
    {
      seq_pos += this->data[rle_pos] / SIGMA; // The last position in the run.
      if(seq_pos >= i) { return this->data[rle_pos] % SIGMA; }
      seq_pos++; rle_pos++; // Move to the first position in the next run.
    }
  }

private:
  int_vector<8> data;
  int_vector<0> samples;
  sd_vector<>   block_boundaries; // Marks the last sequence position in each block.
  sd_vector<>::rank_1_type block_rank;
  sd_vector<>::select_1_type block_select;

  void copy(const RLSequence& v);
  void buildRank();
};  // class RLSequence

//------------------------------------------------------------------------------


#endif // _RELATIVE_FM_SEQUENCE_H
