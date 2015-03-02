#include <sstream>

#include "sequence.h"

namespace relative
{

//------------------------------------------------------------------------------

Sequence::Sequence()
{
  this->sigma = 0;
}

Sequence::Sequence(int_vector_buffer<8>& buffer, uint64_t _size, uint64_t _sigma)
{
  this->sigma = _sigma;
  if(this->sigma > 0) { this->data.width(bitlength(this->sigma - 1)); }
  else { this->data.width(8); }
  this->data.resize(_size);
  for(uint64_t i = 0; i < this->size(); i++) { this->data[i] = buffer[i]; }

  if(this->sigma == 0)
  {
    util::bit_compress(this->data);
    this->sigma = *std::max_element(this->data.begin(), this->data.end()) + 1;
  }

  this->buildRank();
}

Sequence::Sequence(const Sequence& s)
{
  this->copy(s);
}

Sequence::Sequence(Sequence&& s)
{
  *this = std::move(s);
}

Sequence::~Sequence()
{
}

void
Sequence::copy(const Sequence& s)
{
  this->data = s.data;
  this->samples = s.samples;
  this->sigma = s.sigma;
}

void
Sequence::swap(Sequence& s)
{
  if(this != &s)
  {
    this->data.swap(s.data);
    this->samples.swap(s.samples);
    std::swap(this->sigma, s.sigma);
  }
}

Sequence&
Sequence::operator=(const Sequence& s)
{
  if(this != &s) { this->copy(s); }
  return *this;
}

Sequence&
Sequence::operator=(Sequence&& s)
{
  if(this != &s)
  {
    this->data = std::move(s.data);
    this->samples = std::move(s.samples);
    this->sigma = s.sigma;
  }
  return *this;
}

uint64_t
Sequence::serialize(std::ostream& out, structure_tree_node* s, std::string name) const
{
  structure_tree_node* child = structure_tree::add_child(s, name, util::class_name(*this));
  uint64_t written_bytes = 0;
  written_bytes += this->data.serialize(out, child, "data");
  written_bytes += this->samples.serialize(out, child, "samples");
  written_bytes += write_member(this->sigma, out, child, "sigma");
  structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
Sequence::load(std::istream& in)
{
  this->data.load(in);
  this->samples.load(in);
  read_member(this->sigma, in);
}

void
Sequence::buildRank()
{
  uint64_t blocks = (this->size() + SAMPLE_RATE - 1) / SAMPLE_RATE;
  int_vector<0> temp((blocks + 1) * this->sigma, 0, bitlength(this->size()));
  for(uint64_t block = 0; block < blocks; block++)
  {
    for(uint64_t c = 0; c < this->sigma; c++) { temp[(block + 1) * this->sigma + c] = temp[block * this->sigma + c]; }
    uint64_t limit = std::min(this->size(), (block + 1) * SAMPLE_RATE);
    for(uint64_t i = block * SAMPLE_RATE; i < limit; i++) { temp[(block + 1) * this->sigma + this->data[i]]++; }
  }
  util::bit_compress(temp); this->samples.swap(temp);
}

//------------------------------------------------------------------------------

RLSequence::RLSequence()
{
}

void
addRun(uint64_t c, uint64_t run, std::vector<uint8_t>& runs)
{
  while(run > 0)
  {
    uint64_t temp = (run > RLSequence::MAX_RUN ? RLSequence::MAX_RUN : run); run -= temp;
    runs.push_back(RLSequence::encode(c, temp));
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
  for(uint64_t i = 1, pos = 0; i <= this->runs(); i++)
  {
    pos += runLength(this->data[i - 1]);
    if(i == this->runs() || i % SAMPLE_RATE == 0) { block_ends.push_back(pos - 1); }
  }

  // Block boundaries.
  uint64_t blocks = block_ends.size();
  {
    sd_vector<> temp(block_ends.begin(), block_ends.end());
    { std::vector<uint64_t> temp_v; temp_v.swap(block_ends); }
    this->block_boundaries.swap(temp);
    util::init_support(this->block_rank, &(this->block_boundaries));
    util::init_support(this->block_select, &(this->block_boundaries));
  }

  // Samples. SIGMA passes vs. higher space usage? Multi-threaded?
  for(uint64_t c = 0; c < SIGMA; c++)
  {
    int_vector<0> temp(blocks, 0, bitlength(this->size()));
    for(uint64_t block = 0; block < blocks; block++)
    {
      uint64_t limit = std::min(this->runs(), (block + 1) * SAMPLE_RATE);
      for(uint64_t i = block * SAMPLE_RATE; i < limit; i++)
      {
        if(charValue(this->data[i]) == c) { temp[block] += runLength(this->data[i]); }
      }
    }
    util::assign(this->samples[c], CumulativeArray(temp, temp.size()));
  }
}

uint64_t
RLSequence::hash(const Alphabet& alpha) const
{
  uint64_t val = FNV_OFFSET_BASIS;
  for(uint64_t i = 0; i < this->runs(); i++)
  {
    uint64_t c = alpha.comp2char[charValue(this->data[i])], l = runLength(this->data[i]);
    for(uint64_t j = 0; j < l; j++) { val = fnv1a_hash(c, val); }
  }
  return val;
}

//------------------------------------------------------------------------------

template<>
SimpleFM<RLSequence>::SimpleFM(const std::string& base_name, LoadMode mode)
{
  this->sample_rate = 0;

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

  for(uint64_t i = 0; i < sequence.runs(); i++)
  {
    uint8_t temp = sequence.rawData(i);
    counts[RLSequence::charValue(temp)] += RLSequence::runLength(temp);
  }
}

//------------------------------------------------------------------------------

} // namespace relative
