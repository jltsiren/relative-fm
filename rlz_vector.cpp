#include "rlz_vector.h"

//------------------------------------------------------------------------------

RLZVector::RLZVector(const bit_vector& text, const bit_vector& _reference,
  const bit_vector::rank_1_type& _ref_rank,
  const bit_vector::select_1_type& _ref_select) :
  reference(_reference), ref_rank(_ref_rank), ref_select(_ref_select)
{
  std::vector<uint64_t> phrase_starts, phrase_lengths;
  relativeLZ(text, this->reference, phrase_starts, phrase_lengths, this->mismatches);

  std::vector<uint64_t> phrase_buffer;
  this->phrase_rle.resize(phrase_starts.size()); util::set_to_value(this->phrase_rle, 0);
  bit_vector len_vec(text.size());
  bit_vector one_vec(util::cnt_one_bits(text));

  uint64_t bits = 0, onebits = 0, prev = ~(uint64_t)0, max_val = 0;
  for(uint64_t phrase = 0; phrase < phrase_starts.size(); phrase++)
  {
    uint64_t temp = relativeEncoding(phrase_starts[phrase], bits);
    if(temp != prev && phrase_lengths[phrase] > 1)
    {
      phrase_buffer.push_back(temp); prev = temp;
      this->phrase_rle[phrase] = 1;
      max_val = std::max(max_val, temp);
    }
    bits += phrase_lengths[phrase]; len_vec[bits - 1] = 1;  // Last bit in this phrase.
    onebits += this->oneBits(phrase_starts[phrase], phrase_lengths[phrase] - 1);
    if(this->mismatches[phrase]) { onebits++; }
    one_vec[onebits - 1] = 1; // Last 1-bit in this phrase.
  }
  this->phrases.width(bits::hi(max_val) + 1); this->phrases.resize(phrase_buffer.size());
  for(uint64_t i = 0; i < phrase_buffer.size(); i++) { this->phrases[i] = phrase_buffer[i]; }
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
  this->phrase_rle.load(input);
  this->lengths.load(input);
  this->mismatches.load(input);
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
  bytes += size_in_bytes(this->phrases) + size_in_bytes(this->phrase_rle) + size_in_bytes(this->phrase_rank);
  bytes += size_in_bytes(this->lengths) + size_in_bytes(this->length_rank) + size_in_bytes(this->length_select);
  bytes += size_in_bytes(this->mismatches);
  bytes += size_in_bytes(this->ones) + size_in_bytes(this->one_rank) + size_in_bytes(this->one_select);
  return bytes;
}


void
RLZVector::writeTo(std::ofstream& output) const
{
  this->phrases.serialize(output);
  this->phrase_rle.serialize(output);
  this->lengths.serialize(output);
  this->mismatches.serialize(output);
  this->ones.serialize(output);
}

void
RLZVector::buildRankSelect()
{
  util::init_support(this->phrase_rank, &(this->phrase_rle));
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

  uint64_t phrase = this->length_rank(i);               // Phrase is 0-based.
  if(phrase == 0) { return this->oneBits(this->refPos(0, 0), i); }
  uint64_t text_pos = this->length_select(phrase) + 1;  // Starting position of the phrase in the text.

  return this->one_select(phrase) + 1 + this->oneBits(this->refPos(phrase, text_pos), i - text_pos);
}

uint64_t
RLZVector::select(uint64_t i) const
{
  // FIXME what happens for i == 0?
  if(i > this->items()) { return this->size(); }

  uint64_t phrase = this->one_rank(i - 1);              // i is 1-based, phrase is 0-based.
  // Check if the requested 1-bit is the mismatching bit at the end of the phrase.
  if(this->ones[i - 1] && this->mismatches[phrase]) { return this->length_select(phrase + 1); }
  if(phrase == 0) { return this->findBit(this->refPos(0, 0), i); }
  i -= this->one_select(phrase) + 1;                    // i is now relative to the phrase.
  uint64_t text_pos = this->length_select(phrase) + 1;  // Starting position of the phrase in the text.

  return text_pos + this->findBit(this->refPos(phrase, text_pos), i);
}

bool
RLZVector::operator[](uint64_t i) const
{
  if(i >= this->size()) { return false; }

  uint64_t phrase = this->length_rank(i);                   // Phrase is 0-based.
  if(this->lengths[i]) { return this->mismatches[phrase]; } // We wanted the mismatching bit.

  return this->reference[this->refPos(phrase, i)];
}

//------------------------------------------------------------------------------
