#include "sequence.h"


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

template<>
SimpleFM<Sequence>::SimpleFM(const std::string& base_name)
{
  {
    std::string filename = base_name + ALPHA_EXTENSION;
    std::ifstream in(filename.c_str(), std::ios_base::binary);
    if(!in)
    {
      std::cerr << "SimpleFM()::SimpleFM(): Cannot open alphabet file " << filename << std::endl;
      return;
    }
    this->alpha.load(in);
    in.close();
  }

  {
    int_vector_buffer<8> buffer(base_name + BWT_EXTENSION);
    Sequence temp(buffer, buffer.size(), this->alpha.sigma);
    this->bwt.swap(temp);
  }
}

//------------------------------------------------------------------------------

RLSequence::RLSequence()
{
}

void
addRun(uint64_t c, uint64_t run, uint64_t i,
       std::vector<uint8_t>& runs, std::vector<uint64_t>& block_ends, bool force_block)
{
  while(run > 0)
  {
    uint64_t temp = (run > RLSequence::MAX_RUN ? RLSequence::MAX_RUN : run); run -= temp;
    runs.push_back(c + RLSequence::SIGMA * (temp - 1));
    if((force_block && run == 0) || runs.size() % RLSequence::SAMPLE_RATE == 0)
    {
      block_ends.push_back(i - run - 1);
    }
  }
}

RLSequence::RLSequence(int_vector_buffer<8>& buffer, uint64_t _size)
{
  // Process the input.
  uint64_t c = 0, run = 0;
  std::vector<uint64_t> block_ends;
  for(uint64_t i = 0; i < _size; i++)
  {
    if(buffer[i] == c) { run++; }
    else
    {
      addRun(c, run, i, this->data, block_ends, false);
      c = buffer[i]; run = 1;
    }
  }
  addRun(c, run, _size, this->data, block_ends, true);

  // Block boundaries.
  sd_vector<> temp(block_ends.begin(), block_ends.end());
  this->block_boundaries.swap(temp);
  util::init_support(this->block_rank, &(this->block_boundaries));
  util::init_support(this->block_select, &(this->block_boundaries));

  // Rank.
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
  this->samples = s.samples;
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
    this->samples.swap(s.samples);
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
    this->samples = std::move(s.samples);
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
  written_bytes += this->samples.serialize(out, child, "samples");
  written_bytes += this->block_boundaries.serialize(out, child, "block_boundaries");
  written_bytes += this->block_rank.serialize(out, child, "block_rank");
  written_bytes += this->block_select.serialize(out, child, "block_select");
  structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
RLSequence::load(std::istream& in)
{
  read_vector(this->data, in);
  this->samples.load(in);
  this->block_boundaries.load(in);
  this->block_rank.load(in, &(this->block_boundaries));
  this->block_select.load(in, &(this->block_boundaries));
}

void
RLSequence::buildRank()
{
  uint64_t blocks = this->block_rank(this->size());
  int_vector<0> temp((blocks + 1) * SIGMA, 0, bitlength(this->size()));
  for(uint64_t block = 0; block < blocks; block++)
  {
    for(uint64_t c = 0; c < SIGMA; c++) { temp[(block + 1) * SIGMA + c] = temp[block * SIGMA + c]; }
    uint64_t limit = std::min(this->runs(), (block + 1) * SAMPLE_RATE);
    for(uint64_t i = block * SAMPLE_RATE; i < limit; i++)
    {
      temp[(block + 1) * SIGMA + this->data[i] % SIGMA] += this->data[i] / SIGMA + 1;
    }
  }
  util::bit_compress(temp); this->samples.swap(temp);
}

//------------------------------------------------------------------------------
