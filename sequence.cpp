#include "sequence.h"


namespace sdsl
{

//------------------------------------------------------------------------------

Sequence::Sequence()
{
  this->sigma = 0;
}

Sequence::Sequence(int_vector_buffer<8>& buffer, uint64_t _size)
{
  this->data.width(8); this->data.resize(_size);
  for(uint64_t i = 0; i < this->size(); i++) { this->data[i] = buffer[i]; }
  util::bit_compress(this->data);
  this->sigma = *std::max_element(this->data.begin(), this->data.end()) + 1;
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

} // namespace sdsl
