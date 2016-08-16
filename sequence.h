/*
  Copyright (c) 2015, 2016 Genome Research Ltd.

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
  typedef relative::size_type size_type;
  const static size_type SAMPLE_RATE = 64;
  const static size_type SIGMA = 6;
  const static size_type MAX_RUN = 256 / SIGMA;  // 42; encoded as 6 * 41

  RLSequence();
  RLSequence(sdsl::int_vector_buffer<8>& buffer, size_type _size);
  RLSequence(const RLSequence& s);
  RLSequence(RLSequence&& s);
  ~RLSequence();

  void swap(RLSequence& v);
  RLSequence& operator=(const RLSequence& v);
  RLSequence& operator=(RLSequence&& v);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in, LoadMode mode = mode_native);

  inline size_type size() const { return this->block_boundaries.size(); }
  inline size_type bytes() const { return this->data.size(); }
  inline size_type count(uint8_t c) const { return this->samples[c].sum(); }

  // FIXME implement
  size_type countRuns() const;

  inline size_type rank(size_type i, uint8_t c) const
  {
    if(c >= SIGMA) { return 0; }
    if(i > this->size()) { i = this->size(); }

    size_type block = this->block_rank(i);
    size_type res = this->samples[c].sum(block);
    size_type rle_pos = block * SAMPLE_RATE;
    size_type seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);

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

  inline size_type select(size_type i, uint8_t c) const
  {
    if(c >= SIGMA) { return 0; }
    if(i == 0) { return 0; }
    if(i > this->count(c)) { return this->size(); }

    size_type block = this->samples[c].inverse(i - 1);
    size_type count = this->samples[c].sum(block);
    size_type rle_pos = block * SAMPLE_RATE;
    size_type seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);
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

  inline size_type operator[](size_type i) const
  {
    if(i >= this->size()) { return 0; }

    size_type block = this->block_rank(i);
    size_type rle_pos = block * SAMPLE_RATE;
    size_type seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);
    while(true)
    {
      range_type run = this->readRun(rle_pos);
      seq_pos += run.second - 1;  // The last position in the run.
      if(seq_pos >= i) { return run.first; }
      seq_pos++;  // Move to the first position in the next run.
    }
  }

  // returns (rank(i, seq[i]), seq[i])
  inline range_type inverse_select(size_type i) const
  {
    range_type res(0, 0);
    if(i >= this->size()) { return res; }

    size_type block = this->block_rank(i);
    size_type rle_pos = block * SAMPLE_RATE;
    size_type seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);

    range_type run(0, 0);
    size_type ranks[SIGMA] = {};
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
    if(Range::empty(range) || range.second >= this->size()) { return; }
    buffer.resize(Range::length(range));

    // Find the first character.
    size_type block = this->block_rank(range.first);
    size_type rle_pos = block * SAMPLE_RATE;
    size_type seq_pos = (block > 0 ? this->block_select(block) + 1 : 0);
    range_type run(0, 0);

    while(true)
    {
      run = this->readRun(rle_pos);
      seq_pos += run.second - 1;  // The last position in the run.
      if(seq_pos >= range.first) { break; }
      seq_pos++;  // Move to the first position in the next run.
    }

    // Fill the buffer.
    for(size_type i = range.first; i <= range.second; i++)
    {
      if(i > seq_pos)
      {
        run = this->readRun(rle_pos);
        seq_pos += run.second;
      }
      buffer[i - range.first] = run.first;
    }
  }

  inline uint8_t rawData(size_type i) const { return this->data[i]; }

  size_type hash(const Alphabet& alpha) const;

  // Returns (character, length).
  inline range_type readRun(size_type& rle_pos) const
  {
    range_type run(this->data[rle_pos] % SIGMA, this->data[rle_pos] / SIGMA + 1); rle_pos++;
    if(run.second >= MAX_RUN)
    {
      size_type temp = this->data[rle_pos] & 0x7F, offset = 7;
      while(this->data[rle_pos] & 0x80)
      {
        rle_pos++; temp |= ((size_type)(this->data[rle_pos] & 0x7F)) << offset; offset += 7;
      }
      run.second += temp; rle_pos++;
    }
    return run;
  }

private:
  std::vector<uint8_t> data;
  CumulativeArray      samples[SIGMA];

  sdsl::sd_vector<>                block_boundaries; // Marks the last sequence position in each block.
  sdsl::sd_vector<>::rank_1_type   block_rank;
  sdsl::sd_vector<>::select_1_type block_select;

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
  for(size_type i = 0; i < buffer.size(); i++) { buffer[i] = this->alpha.comp2char[buffer[i]]; }
}

template<>
void characterCounts(const RLSequence& sequence, sdsl::int_vector<64>& counts);

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_SEQUENCE_H
