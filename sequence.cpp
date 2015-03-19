#include <sstream>

#include "sequence.h"

namespace relative
{

//------------------------------------------------------------------------------

RLSequence::RLSequence()
{
}

inline void
addBasicRun(uint64_t c, uint64_t run, std::vector<uint8_t>& runs)
{
  runs.push_back(c + RLSequence::SIGMA * (run - 1));
}

void
addRun(uint64_t c, uint64_t run, std::vector<uint8_t>& runs)
{
  while(run > 0)
  {
    // Runs of length <= 41 are encoded in a single byte.
    if(run < RLSequence::MAX_RUN)
    {
      addBasicRun(c, run, runs);
      return;
    }

    uint64_t bytes_remaining = RLSequence::SAMPLE_RATE - (runs.size() % RLSequence::SAMPLE_RATE);
    uint64_t l = (bytes_remaining < 2 ? RLSequence::MAX_RUN - 1 : RLSequence::MAX_RUN);
    addBasicRun(c, l, runs); run -= l;
    if(bytes_remaining < 2) { continue; } // No room for the additional run length in the current block.
    bytes_remaining--;

    // Write the rest of the run with 7 bits/byte, least significant byte first.
    if(bitlength(run) > 7 * bytes_remaining)  // Cannot encode the entire run in current block.
    {
      uint64_t temp = 0;
      for(uint64_t i = 1; i < bytes_remaining; i++) { runs.push_back(0xFF); temp = (temp << 7) | 0x7F; }
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

RLSequence::RLSequence(int_vector_buffer<8>& buffer, uint64_t _size)
{
  // Process the input.
  uint64_t c = 0, run = 0;
  for(uint64_t i = 0; i < _size; i++)
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
  for(uint64_t c = 0; c < SIGMA; c++) { this->samples[c] = s.samples[c]; }
  this->block_boundaries = s.block_boundaries;
  util::init_support(this->block_rank, &(this->block_boundaries));
  util::init_support(this->block_select, &(this->block_boundaries));
}

void
RLSequence::swap(RLSequence& s)
{
  if(this != &s)
  {
    this->data.swap(s.data);
    for(uint64_t c = 0; c < SIGMA; c++) { this->samples[c].swap(s.samples[c]); }
    this->block_boundaries.swap(s.block_boundaries);
    util::swap_support(this->block_rank, s.block_rank, &(this->block_boundaries), &(s.block_boundaries));
    util::swap_support(this->block_select, s.block_select, &(this->block_boundaries), &(s.block_boundaries));
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
    for(uint64_t c = 0; c < SIGMA; c++) { this->samples[c] = std::move(s.samples[c]); }
    this->block_boundaries = std::move(s.block_boundaries);
    util::init_support(this->block_rank, &(this->block_boundaries));
    util::init_support(this->block_select, &(this->block_boundaries));
  }
  return *this;
}

uint64_t
RLSequence::serialize(std::ostream& out, structure_tree_node* s, std::string name) const
{
  structure_tree_node* child = structure_tree::add_child(s, name, util::class_name(*this));
  uint64_t written_bytes = 0;
  written_bytes += write_vector(this->data, out, child, "data");
  for(uint64_t c = 0; c < SIGMA; c++)
  {
    std::stringstream ss; ss << "samples_" << c;
    written_bytes += this->samples[c].serialize(out, child, ss.str());
  }
  written_bytes += this->block_boundaries.serialize(out, child, "block_boundaries");
  written_bytes += this->block_rank.serialize(out, child, "block_rank");
  written_bytes += this->block_select.serialize(out, child, "block_select");
  structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
RLSequence::load(std::istream& in, bool rebuild_samples)
{
  read_vector(this->data, in);
  if(rebuild_samples)
  {
    this->buildRank();
  }
  else
  {
    for(uint64_t c = 0; c < SIGMA; c++) { this->samples[c].load(in); }
    this->block_boundaries.load(in);
    this->block_rank.load(in, &(this->block_boundaries));
    this->block_select.load(in, &(this->block_boundaries));
  }
}

void
RLSequence::buildRank()
{
  std::vector<uint64_t> block_ends;
  uint64_t seq_pos = 0, rle_pos = 0;
  while(rle_pos < this->runs())
  {
    range_type run = this->readRun(rle_pos); seq_pos += run.second;
    if(rle_pos >= this->runs() || rle_pos % SAMPLE_RATE == 0) { block_ends.push_back(seq_pos - 1); }
  }

  // Block boundaries.
  uint64_t blocks = block_ends.size();
  {
    sd_vector<> temp(block_ends.begin(), block_ends.end());
    util::clear(block_ends);
    this->block_boundaries.swap(temp);
    util::init_support(this->block_rank, &(this->block_boundaries));
    util::init_support(this->block_select, &(this->block_boundaries));
  }

  int_vector<0> counts[SIGMA];
  for(uint64_t c = 0; c < SIGMA; c++) { util::assign(counts[c], int_vector<0>(blocks, 0, bitlength(this->size()))); }
  for(uint64_t block = 0; block < blocks; block++)
  {
    uint64_t rle_pos = block * SAMPLE_RATE, limit = std::min(this->runs(), (block + 1) * SAMPLE_RATE);
    while(rle_pos < limit)
    {
      range_type run = this->readRun(rle_pos);
      counts[run.first][block] += run.second;
    }
  }
  for(uint64_t c = 0; c < SIGMA; c++) { util::assign(this->samples[c], CumulativeArray(counts[c])); }
}

uint64_t
RLSequence::hash(const Alphabet& alpha) const
{
  uint64_t rle_pos = 0, val = FNV_OFFSET_BASIS;
  while(rle_pos < this->runs())
  {
    range_type run = this->readRun(rle_pos);
    uint8_t c = alpha.comp2char[run.first];
    for(uint64_t i = 0; i < run.second; i++) { val = fnv1a_hash(c, val); }
  }
  return val;
}

//------------------------------------------------------------------------------

template<>
SimpleFM<RLSequence>::SimpleFM(const std::string& base_name, LoadMode mode)
{
  this->sa_sample_rate = 0; this->isa_sample_rate = 0;

  if(mode == mode_native || mode == mode_ropebwt2)
  {
    std::string filename = base_name + NATIVE_BWT_EXTENSION;
    std::ifstream in(filename.c_str(), std::ios_base::binary);
    if(!in)
    {
      std::cerr << "SimpleFM::SimpleFM(): Cannot open BWT file " << filename << " (native format)" << std::endl;
      return;
    }
    this->bwt.load(in, mode == mode_ropebwt2); in.close();
  }
  else
  {
    int_vector_buffer<8> buffer(base_name + BWT_EXTENSION);
    RLSequence temp(buffer, buffer.size());
    this->bwt.swap(temp);
  }

  this->loadAlphabet(base_name);
  this->loadSamples(base_name);
}

template<>
void
characterCounts(const RLSequence& sequence, uint64_t size, int_vector<64>& counts)
{
  for(uint64_t c = 0; c < counts.size(); c++) { counts[c] = 0; }

  uint64_t rle_pos = 0, seq_pos = 0;
  while(rle_pos < sequence.runs())
  {
    range_type run = sequence.readRun(rle_pos);
    counts[run.first] += run.second; seq_pos += run.second;
    if(seq_pos >= size) { counts[run.first] -= seq_pos - size; break; }
  }
}

//------------------------------------------------------------------------------

} // namespace relative
