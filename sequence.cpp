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

#include <sstream>

#include "sequence.h"

namespace relative
{

//------------------------------------------------------------------------------

RLSequence::RLSequence()
{
}

inline void
addBasicRun(size_type c, size_type run, std::vector<uint8_t>& runs)
{
  runs.push_back(c + RLSequence::SIGMA * (run - 1));
}

void
addRun(size_type c, size_type run, std::vector<uint8_t>& runs)
{
  while(run > 0)
  {
    // Runs of length <= 41 are encoded in a single byte.
    if(run < RLSequence::MAX_RUN)
    {
      addBasicRun(c, run, runs);
      return;
    }

    size_type bytes_remaining = RLSequence::SAMPLE_RATE - (runs.size() % RLSequence::SAMPLE_RATE);
    size_type l = (bytes_remaining < 2 ? RLSequence::MAX_RUN - 1 : RLSequence::MAX_RUN);
    addBasicRun(c, l, runs); run -= l;
    if(bytes_remaining < 2) { continue; } // No room for the additional run length in the current block.
    bytes_remaining--;

    // Write the rest of the run with 7 bits/byte, least significant byte first.
    if(bit_length(run) > 7 * bytes_remaining)  // Cannot encode the entire run in current block.
    {
      size_type temp = 0;
      for(size_type i = 1; i < bytes_remaining; i++) { runs.push_back(0xFF); temp = (temp << 7) | 0x7F; }
      runs.push_back(0x7F); run -= (temp << 7) | 0x7F;
    }
    else
    {
      do
      {
        uint8_t val = run & 0x7F; run >>= 7;
        if(run > 0) { val |= 0x80; }
        runs.push_back(val);
      }
      while(run > 0);
    }
  }
}

RLSequence::RLSequence(sdsl::int_vector_buffer<8>& buffer, size_type _size)
{
  // Process the input.
  size_type c = 0, run = 0;
  for(size_type i = 0; i < _size; i++)
  {
    if(buffer[i] == c) { run++; }
    else
    {
      addRun(c, run, this->data);
      c = buffer[i]; run = 1;
    }
  }
  addRun(c, run, this->data);

  // Blocks and samples.
  this->buildRank();
}

RLSequence::RLSequence(const RLSequence& s)
{
  this->copy(s);
}

RLSequence::RLSequence(RLSequence&& s)
{
  *this = std::move(s);
}

RLSequence::~RLSequence()
{
}

void
RLSequence::copy(const RLSequence& s)
{
  this->data = s.data;
  for(size_type c = 0; c < SIGMA; c++) { this->samples[c] = s.samples[c]; }
  this->block_boundaries = s.block_boundaries;
  sdsl::util::init_support(this->block_rank, &(this->block_boundaries));
  sdsl::util::init_support(this->block_select, &(this->block_boundaries));
}

void
RLSequence::swap(RLSequence& s)
{
  if(this != &s)
  {
    this->data.swap(s.data);
    for(size_type c = 0; c < SIGMA; c++) { this->samples[c].swap(s.samples[c]); }
    this->block_boundaries.swap(s.block_boundaries);
    sdsl::util::swap_support(this->block_rank, s.block_rank, &(this->block_boundaries), &(s.block_boundaries));
    sdsl::util::swap_support(this->block_select, s.block_select, &(this->block_boundaries), &(s.block_boundaries));
  }
}

RLSequence&
RLSequence::operator=(const RLSequence& s)
{
  if(this != &s) { this->copy(s); }
  return *this;
}

RLSequence&
RLSequence::operator=(RLSequence&& s)
{
  if(this != &s)
  {
    this->data = std::move(s.data);
    for(size_type c = 0; c < SIGMA; c++) { this->samples[c] = std::move(s.samples[c]); }
    this->block_boundaries = std::move(s.block_boundaries);
    sdsl::util::init_support(this->block_rank, &(this->block_boundaries));
    sdsl::util::init_support(this->block_select, &(this->block_boundaries));
  }
  return *this;
}

size_type
RLSequence::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += write_vector(this->data, out, child, "data");
  for(size_type c = 0; c < SIGMA; c++)
  {
    std::stringstream ss; ss << "samples_" << c;
    written_bytes += this->samples[c].serialize(out, child, ss.str());
  }
  written_bytes += this->block_boundaries.serialize(out, child, "block_boundaries");
  written_bytes += this->block_rank.serialize(out, child, "block_rank");
  written_bytes += this->block_select.serialize(out, child, "block_select");
  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
RLSequence::load(std::istream& in, LoadMode mode)
{
  if(mode == mode_plain)
  {
    // FIXME implement?
    std::cerr << "RLSequence::load(): Use the constructor for mode_plain." << std::endl;
    return;
  }
  else if(mode == mode_ropebwt)
  {
    size_t _sequences = 0, _symbols = 0, _runs = 0;
    {
      uint16_t magic_constant = 0;
      sdsl::read_member(magic_constant, in);
      if(magic_constant != 0xCACA) { std::cerr << "RLSequence::load(): Invalid magic constant!" << std::endl; }
      sdsl::read_member(_sequences, in);
      sdsl::read_member(_symbols, in);
      sdsl::read_member(_runs, in);
      std::cout << "Sequences: " << _sequences
                << ", symbols: " << _symbols
                << ", runs: " << _runs << std::endl;
      int32_t bwflag = 0;
      sdsl::read_member(bwflag, in);
      if(bwflag != 0) { std::cerr << "RLSequence::load(): Invalid BWFlag: " << bwflag << std::endl; }
    }
    /*
      FIXME Read body, output runs, bytes
      runs: 1 byte each
        3 high bits for the character
        5 low bits for the actual run length
    */
    std::cerr << "RLSequence::load(): Loading in mode_ropebwt has not been implemented." << std::endl;
    return;
  }
  else
  {
    read_vector(this->data, in);
  }

  if(mode == mode_native)
  {
    for(size_type c = 0; c < SIGMA; c++) { this->samples[c].load(in); }
    this->block_boundaries.load(in);
    this->block_rank.load(in, &(this->block_boundaries));
    this->block_select.load(in, &(this->block_boundaries));
  }
  else
  {
    this->buildRank();
  }
}

void
RLSequence::buildRank()
{
  std::vector<size_type> block_ends;
  size_type seq_pos = 0, rle_pos = 0;
  while(rle_pos < this->bytes())
  {
    range_type run = this->readRun(rle_pos); seq_pos += run.second;
    if(rle_pos >= this->bytes() || rle_pos % SAMPLE_RATE == 0) { block_ends.push_back(seq_pos - 1); }
  }

  // Block boundaries.
  size_type blocks = block_ends.size();
  {
    sdsl::sd_vector<> temp(block_ends.begin(), block_ends.end());
    sdsl::util::clear(block_ends);
    this->block_boundaries.swap(temp);
    sdsl::util::init_support(this->block_rank, &(this->block_boundaries));
    sdsl::util::init_support(this->block_select, &(this->block_boundaries));
  }

  sdsl::int_vector<0> counts[SIGMA];
  for(size_type c = 0; c < SIGMA; c++) { sdsl::util::assign(counts[c], sdsl::int_vector<0>(blocks, 0, bit_length(this->size()))); }
  for(size_type block = 0; block < blocks; block++)
  {
    size_type rle_pos = block * SAMPLE_RATE, limit = std::min(this->bytes(), (block + 1) * SAMPLE_RATE);
    while(rle_pos < limit)
    {
      range_type run = this->readRun(rle_pos);
      counts[run.first][block] += run.second;
    }
  }
  for(size_type c = 0; c < SIGMA; c++) { sdsl::util::assign(this->samples[c], CumulativeArray(counts[c])); }
}

size_type
RLSequence::hash(const Alphabet& alpha) const
{
  size_type rle_pos = 0, val = FNV_OFFSET_BASIS;
  while(rle_pos < this->bytes())
  {
    range_type run = this->readRun(rle_pos);
    uint8_t c = alpha.comp2char[run.first];
    for(size_type i = 0; i < run.second; i++) { val = fnv1a_hash(c, val); }
  }
  return val;
}

//------------------------------------------------------------------------------

template<>
SimpleFM<RLSequence>::SimpleFM(const std::string& base_name, LoadMode mode)
{
  this->sa_sample_rate = 0; this->isa_sample_rate = 0;

  if(mode == mode_plain)
  {
    sdsl::int_vector_buffer<8> buffer(base_name + BWT_EXTENSION);
    RLSequence temp(buffer, buffer.size());
    this->bwt.swap(temp);
  }
  else
  {
    std::string filename = base_name + NATIVE_BWT_EXTENSION;
    std::ifstream in(filename.c_str(), std::ios_base::binary);
    if(!in)
    {
      std::cerr << "SimpleFM::SimpleFM(): Cannot open BWT file " << filename << " (native format)" << std::endl;
      return;
    }
    this->bwt.load(in, mode); in.close();
  }

  this->loadAlphabet(base_name);
  this->loadSamples(base_name);
}

template<>
void
characterCounts(const RLSequence& sequence, sdsl::int_vector<64>& counts)
{
  for(size_type c = 0; c < counts.size(); c++) { counts[c] = 0; }

  size_type rle_pos = 0, seq_pos = 0;
  while(rle_pos < sequence.bytes())
  {
    range_type run = sequence.readRun(rle_pos);
    counts[run.first] += run.second; seq_pos += run.second;
    if(seq_pos >= sequence.size()) { counts[run.first] -= seq_pos - sequence.size(); break; }
  }
}

//------------------------------------------------------------------------------

} // namespace relative
