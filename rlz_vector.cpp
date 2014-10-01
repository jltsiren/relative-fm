#include "rlz_vector.h"

//------------------------------------------------------------------------------

RLZVector::RLZVector(const bit_vector& text, const bit_vector& _reference,
  const bit_vector::rank_1_type& _ref_rank,
  const bit_vector::select_1_type& _ref_select) :
  reference(_reference), ref_rank(_ref_rank), ref_select(_ref_select)
{
  std::vector<range_type> ranges;
  relativeLZ(text, this->reference, ranges);

  bit_vector len_vec(text.size() + 1); len_vec[0] = 1;
  bit_vector one_vec(util::cnt_one_bits(text) + 1); one_vec[0] = 1;
  this->phrases.width(bits::hi(this->reference.size() - 1) + 1);
  this->phrases.resize(ranges.size());

  uint64_t phrase = 0, bits = 0, onebits = 0;
  for(auto range : ranges)
  {
    this->phrases[phrase] = range.first; phrase++;        // Starting points.
    bits += range.second; len_vec[bits] = 1;              // Bits after this phrase.
    onebits += this->oneBits(range.first, range.second);
    one_vec[onebits] = 1;                                 // 1-bits after this phrase.
  }
  this->lengths = len_vec;
  this->ones = one_vec;

  this->buildRankSelect();
}

RLZVector::RLZVector(std::ifstream& input, const bit_vector& _reference,
  const bit_vector::rank_1_type& _ref_rank,
  const bit_vector::select_1_type& _ref_select) :
  reference(_reference), ref_rank(_ref_rank), ref_select(_ref_select)
{
  this->phrases.load(input);
  this->lengths.load(input);
  this->ones.load(input);
  this->buildRankSelect();
}

RLZVector::~RLZVector()
{
}

//------------------------------------------------------------------------------

uint64_t
RLZVector::reportSize() const
{
  uint64_t bytes = sizeof(*this);
  bytes += size_in_bytes(this->phrases);
  bytes += size_in_bytes(this->lengths) + size_in_bytes(this->length_rank) + size_in_bytes(this->length_select);
  bytes += size_in_bytes(this->ones) + size_in_bytes(this->one_rank) + size_in_bytes(this->one_select);
  return bytes;
}


void
RLZVector::writeTo(std::ofstream& output) const
{
  this->phrases.serialize(output);
  this->lengths.serialize(output);
  this->ones.serialize(output);
}

void
RLZVector::buildRankSelect()
{
  util::init_support(this->length_rank, &(this->lengths));
  util::init_support(this->length_select, &(this->lengths));
  util::init_support(this->one_rank, &(this->ones));
  util::init_support(this->one_select, &(this->ones));
}

//------------------------------------------------------------------------------

uint64_t
RLZVector::rank(uint64_t i) const
{
  if(i >= this->size()) { return this->items(); }

  uint64_t phrase = this->length_rank(i + 1);  // Phrase is 1-based.
  uint64_t offset = i - this->length_select(phrase);

  return this->one_select(phrase) + this->oneBits(this->phrases[phrase - 1], offset);
}

uint64_t
RLZVector::select(uint64_t i) const
{
  // FIXME what happens for i == 0?
  if(i > this->items()) { return this->size(); }

  uint64_t phrase = this->one_rank(i);  // Both i and phrase are 1-based.
  uint64_t offset = i - this->one_select(phrase);

  return this->length_select(phrase) + this->findBit(this->phrases[phrase - 1], offset);
}

bool
RLZVector::operator[](uint64_t i) const
{
  if(i >= this->size()) { return false; }

  uint64_t phrase = this->length_rank(i + 1);  // Phrase is 1-based.
  uint64_t offset = i - this->length_select(phrase);

  return this->bitAt(this->phrases[phrase - 1], offset);
}

//------------------------------------------------------------------------------
