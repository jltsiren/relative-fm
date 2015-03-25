#ifndef _RELATIVE_FM_SEQUENCE_H
#define _RELATIVE_FM_SEQUENCE_H

#include "simple_fm.h"

namespace relative
{

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

  RLSequence();
  RLSequence(int_vector_buffer<8>& buffer, uint64_t _size);
  RLSequence(const RLSequence& s);
  RLSequence(RLSequence&& s);
  ~RLSequence();

  void swap(RLSequence& v);
  RLSequence& operator=(const RLSequence& v);
  RLSequence& operator=(RLSequence&& v);

  uint64_t serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in, LoadMode mode = mode_native);

  inline uint64_t size() const { return this->block_boundaries.size(); }
  inline uint64_t bytes() const { return this->data.size(); }
  inline uint64_t count(uint8_t c) const { return this->samples[c].sum(); }

  // FIXME implement
  uint64_t countRuns() const;

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
      range_type run = this->readRun(rle_pos);
      seq_pos += run.second;  // The starting position of the next run.
      if(run.first == c)
      {
        res += run.second;  // Number of c's before the next run.
        if(seq_pos > i) { res -= seq_pos - i; }
      }
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
      range_type run = this->readRun(rle_pos);
      seq_pos += run.second - 1;  // The last position in the run.
      if(run.first == c)
      {
        count += run.second;  // Number of c's up to the end of the run.
        if(count >= i) { return seq_pos + i - count; }
      }
      seq_pos++;  // Move to the first position in the next run.
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
      range_type run = this->readRun(rle_pos);
      seq_pos += run.second - 1;  // The last position in the run.
      if(seq_pos >= i) { return run.first; }
      seq_pos++;  // Move to the first position in the next run.
    }
  }

  // returns (rank(i, seq[i]), seq[i])
  inline range_type inverse_select(uint64_t i) const
  {
    range_type res(0, 0);
    if(i >= this->size()) { return res; }

    uint64_t block = this->block_rank(i);
    uint64_t rle_pos = block * SAMPLE_RATE;
    uint64_t seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);

    range_type run(0, 0);
    uint64_t ranks[SIGMA] = {};
    while(seq_pos <= i)
    {
      run = this->readRun(rle_pos);
      seq_pos += run.second;  // The starting position of the next run.
      ranks[run.first] += run.second; // Number of c's before the next run.
    }

    return range_type(this->samples[run.first].sum(block) + ranks[run.first] - (seq_pos - i), run.first);
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
    range_type run(0, 0);

    while(true)
    {
      run = this->readRun(rle_pos);
      seq_pos += run.second - 1;  // The last position in the run.
      if(seq_pos >= range.first) { break; }
      seq_pos++;  // Move to the first position in the next run.
    }

    // Fill the buffer.
    for(uint64_t i = range.first; i <= range.second; i++)
    {
      if(i > seq_pos)
      {
        run = this->readRun(rle_pos);
        seq_pos += run.second;
      }
      buffer[i - range.first] = run.first;
    }
  }

  inline uint8_t rawData(uint64_t i) const { return this->data[i]; }

  uint64_t hash(const Alphabet& alpha) const;

  // Returns (character, length).
  inline range_type readRun(uint64_t& rle_pos) const
  {
    range_type run(this->data[rle_pos] % SIGMA, this->data[rle_pos] / SIGMA + 1); rle_pos++;
    if(run.second >= MAX_RUN)
    {
      uint64_t temp = this->data[rle_pos] & 0x7F, offset = 7;
      while(this->data[rle_pos] & 0x80)
      {
        rle_pos++; temp |= ((uint64_t)(this->data[rle_pos] & 0x7F)) << offset; offset += 7;
      }
      run.second += temp; rle_pos++;
    }
    return run;
  }

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
void characterCounts(const RLSequence& sequence, int_vector<64>& counts);

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_SEQUENCE_H
