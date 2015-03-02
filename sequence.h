#ifndef _RELATIVE_FM_SEQUENCE_H
#define _RELATIVE_FM_SEQUENCE_H

#include "simple_fm.h"

namespace relative
{

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

  void swap(Sequence& s);
  Sequence& operator=(const Sequence& s);
  Sequence& operator=(Sequence&& s);

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

  inline range_type inverse_select(uint64_t i) const
  {
    range_type res(0, (*this)[i]);
    res.first = this->rank(i, res.second);
    return res;
  }

private:
  int_vector<0> data, samples;
  uint64_t      sigma;

  void copy(const Sequence& s);
  void buildRank();
};  // class Sequence

//------------------------------------------------------------------------------

/*
  A basic run-length encoded sequence for alphabet 0-5.
*/
class RLSequence
{
public:
  const static uint64_t SAMPLE_RATE = 64;
  const static uint64_t SIGMA = 6;
  const static uint64_t MAX_RUN = 256 / SIGMA;  // 42; encoded as 6 * 41
  typedef uint64_t size_type;

  inline static uint64_t charValue(uint8_t code) { return code % SIGMA; }
  inline static uint64_t runLength(uint8_t code) { return code / SIGMA + 1; }
  inline static uint8_t encode(uint64_t c, uint64_t run) { return c + SIGMA * (run - 1); }

  RLSequence();
  RLSequence(int_vector_buffer<8>& buffer, uint64_t _size);
  RLSequence(const RLSequence& s);
  RLSequence(RLSequence&& s);
  ~RLSequence();

  void swap(RLSequence& v);
  RLSequence& operator=(const RLSequence& v);
  RLSequence& operator=(RLSequence&& v);

  uint64_t serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in, bool rebuild_samples = false);

  inline uint64_t size() const { return this->block_boundaries.size(); }
  inline uint64_t runs() const { return this->data.size(); }
  inline uint64_t count(uint8_t c) const { return this->samples[c].sum(); }

  inline uint64_t rank(uint64_t i, uint8_t c) const
  {
    if(c >= SIGMA) { return 0; }
    if(i > this->size()) { i = this->size(); }

    uint64_t block = this->block_rank(i);
    uint64_t res = this->samples[c].sum(block);
    uint64_t rle_pos = block * SAMPLE_RATE;
    uint64_t seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);
    while(seq_pos < i)
    {
      seq_pos += runLength(this->data[rle_pos]); // The starting position of the next run.
      if(charValue(this->data[rle_pos]) == c)
      {
        res += runLength(this->data[rle_pos]); // Number of c's before the next run.
        if(seq_pos > i) { res -= seq_pos - i; }
      }
      rle_pos++;
    }

    return res;
  }

  inline uint64_t select(uint64_t i, uint8_t c) const
  {
    if(c >= SIGMA) { return 0; }
    if(i == 0) { return 0; }
    if(i > this->count(c)) { return this->size(); }

    uint64_t block = this->samples[c].inverse(i - 1);
    uint64_t count = this->samples[c].sum(block);
    uint64_t rle_pos = block * SAMPLE_RATE;
    uint64_t seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);
    while(true)
    {
      seq_pos += runLength(this->data[rle_pos]) - 1; // The last position in the run.
      if(charValue(this->data[rle_pos]) == c)
      {
        count += runLength(this->data[rle_pos]);  // Number of c's up to the end of the run.
        if(count >= i) { return seq_pos + i - count; }
      }
      seq_pos++; rle_pos++; // Move to the first position in the next run.
    }
  }

  inline uint64_t operator[](uint64_t i) const
  {
    if(i >= this->size()) { return 0; }

    uint64_t block = this->block_rank(i);
    uint64_t rle_pos = block * SAMPLE_RATE;
    uint64_t seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);
    while(true)
    {
      seq_pos += runLength(this->data[rle_pos]) - 1; // The last position in the run.
      if(seq_pos >= i) { return charValue(this->data[rle_pos]); }
      seq_pos++; rle_pos++; // Move to the first position in the next run.
    }
  }

  inline range_type inverse_select(uint64_t i) const
  {
    range_type res(0, 0);
    if(i >= this->size()) { return res; }

    uint64_t block = this->block_rank(i);
    uint64_t rle_pos = block * SAMPLE_RATE;
    uint64_t seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);

    // Determine the character.
    uint64_t block_start = seq_pos;
    while(true)
    {
      seq_pos += runLength(this->data[rle_pos]) - 1;  // The last position in the run.
      if(seq_pos >= i) { res.second = charValue(this->data[rle_pos]); }
      seq_pos++; rle_pos++; // Move to the first position in the next run.
    }

    // Determine the rank.
    res.first = this->samples[res.second].sum(block);
    rle_pos = block * SAMPLE_RATE; seq_pos = block_start;
    while(seq_pos < i)  // Determine the rank.
    {
      seq_pos += runLength(this->data[rle_pos]); // The starting position of the next run.
      if(charValue(this->data[rle_pos]) == res.second)
      {
        res.first += runLength(this->data[rle_pos]); // Number of c's before the next run.
        if(seq_pos > i) { res.first -= seq_pos - i; }
      }
      rle_pos++;
    }

    return res;
  }

  template<class ByteVector>
  void extract(range_type range, ByteVector& buffer) const
  {
    if(isEmpty(range) || range.second >= this->size()) { return; }
    buffer.resize(length(range));

    // Find the first character.
    uint64_t block = this->block_rank(range.first);
    uint64_t rle_pos = block * SAMPLE_RATE;
    uint64_t seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);
    while(true)
    {
      seq_pos += runLength(this->data[rle_pos]) - 1; // The last position in the run.
      if(seq_pos >= range.first) { break; }
      seq_pos++; rle_pos++; // Move to the first position in the next run.
    }

    // Fill the buffer.
    for(uint64_t i = range.first; i <= range.second; i++)
    {
      if(i > seq_pos) { rle_pos++; seq_pos += runLength(this->data[rle_pos]); }
      buffer[i - range.first] = charValue(this->data[rle_pos]);
    }
  }

  inline uint8_t rawData(uint64_t i) const { return this->data[i]; }

  uint64_t hash(const Alphabet& alpha) const;

private:
  std::vector<uint8_t> data;
  CumulativeArray      samples[SIGMA];

  sd_vector<>                block_boundaries; // Marks the last sequence position in each block.
  sd_vector<>::rank_1_type   block_rank;
  sd_vector<>::select_1_type block_select;

  void copy(const RLSequence& v);
  void buildRank();
};  // class RLSequence

template<>
SimpleFM<RLSequence>::SimpleFM(const std::string& base_name, LoadMode mode);

template<>
template<class ByteVector>
void
SimpleFM<RLSequence>::extractBWT(range_type range, ByteVector& buffer) const
{
  this->bwt.extract(range, buffer);
  for(uint64_t i = 0; i < buffer.size(); i++) { buffer[i] = this->alpha.comp2char[buffer[i]]; }
}

template<>
void characterCounts(const RLSequence& sequence, uint64_t size, int_vector<64>& counts);

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_SEQUENCE_H
