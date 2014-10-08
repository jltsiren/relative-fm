#include "rlz.h"
#include "rlz_vector.h"


namespace sdsl
{

//------------------------------------------------------------------------------

rlz_vector::rlz_vector()
{
}

rlz_vector::rlz_vector(const rlz_vector& v)
{
  this->copy(v);
}

rlz_vector::rlz_vector(rlz_vector&& v)
{
  *this = std::move(v);
}

rlz_vector::rlz_vector(const bit_vector& v)
{
  this->m_size = v.size();
  this->m_plain = v;
  this->init_support();
}

rlz_vector::rlz_vector(bit_vector&& v)
{
  this->m_size = v.size();
  this->m_plain = std::move(v);
  this->init_support();
}

void
rlz_vector::copy(const rlz_vector& v)
{
  this->m_size = v.m_size;
  this->m_plain = v.m_plain;
  this->m_plain_rank = v.m_plain_rank; this->m_plain_rank.set_vector(&(this->m_plain));
  this->m_plain_select_1 = v.m_plain_select_1; this->m_plain_select_1.set_vector(&(this->m_plain));
  this->m_plain_select_0 = v.m_plain_select_0; this->m_plain_select_0.set_vector(&(this->m_plain));
}

void
rlz_vector::init_support()
{
  util::init_support(this->m_plain_rank, &(this->m_plain));
  util::init_support(this->m_plain_select_1, &(this->m_plain));
  util::init_support(this->m_plain_select_0, &(this->m_plain));
}

void
rlz_vector::swap(rlz_vector& v)
{
  if(this != &v)
  {
    std::swap(this->m_size, v.m_size);
    this->m_plain.swap(v.m_plain);
    util::swap_support(this->m_plain_rank, v.m_plain_rank, &(this->m_plain), &(v.m_plain));
    util::swap_support(this->m_plain_select_1, v.m_plain_select_1, &(this->m_plain), &(v.m_plain));
    util::swap_support(this->m_plain_select_0, v.m_plain_select_0, &(this->m_plain), &(v.m_plain));
  }
}

rlz_vector&
rlz_vector::operator=(const rlz_vector& v)
{
  if(this != &v) { this->copy(v); }
  return *this;
}

rlz_vector&
rlz_vector::operator=(rlz_vector&& v)
{
  if(this != &v)
  {
    this->m_size = v.m_size;
    this->m_plain = std::move(v.m_plain);
    this->m_plain_rank = std::move(v.m_plain_rank); this->m_plain_rank.set_vector(&(this->m_plain));
    this->m_plain_select_1 = std::move(v.m_plain_select_1); this->m_plain_select_1.set_vector(&(this->m_plain));
    this->m_plain_select_0 = std::move(v.m_plain_select_0); this->m_plain_select_0.set_vector(&(this->m_plain));
  }
  return *this;
}

rlz_vector::size_type
rlz_vector::serialize(std::ostream& out, structure_tree_node* v, std::string name) const
{
  structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += write_member(this->m_size, out, child, "size");
  written_bytes += this->m_plain.serialize(out, child, "plain");
  written_bytes += this->m_plain_rank.serialize(out, child, "plain_rank");
  written_bytes += this->m_plain_select_1.serialize(out, child, "plain_select_1");
  written_bytes += this->m_plain_select_0.serialize(out, child, "plain_select_0");
  structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
rlz_vector::load(std::istream& in)
{
  read_member(this->m_size, in);
  this->m_plain.load(in);
  this->m_plain_rank.load(in, &(this->m_plain));
  this->m_plain_select_1.load(in, &(this->m_plain));
  this->m_plain_select_0.load(in, &(this->m_plain));
}

//------------------------------------------------------------------------------

RLZVector::RLZVector(const bit_vector& text, const bit_vector& _reference,
  const bit_vector::rank_1_type& _ref_rank,
  const bit_vector::select_1_type& _ref_select) :
  reference(_reference), ref_rank(_ref_rank), ref_select(_ref_select)
{
  std::vector<uint64_t> phrase_starts, phrase_lengths;
  relativeLZSuccinct(text, this->reference, phrase_starts, phrase_lengths, this->mismatches);

  std::vector<uint64_t> phrase_buffer;
  this->phrase_rle.resize(phrase_starts.size()); util::set_to_value(this->phrase_rle, 0);

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
    bits += phrase_lengths[phrase];
    onebits += this->oneBits(phrase_starts[phrase], phrase_lengths[phrase] - 1);
    if(this->mismatches[phrase]) { onebits++; }
    phrase_starts[phrase] = bits - 1;     // The last bit in this phrase.
    phrase_lengths[phrase] = onebits - 1; // The last 1-bit in this phrase.
  }
  this->phrases.width(bits::hi(max_val) + 1); this->phrases.resize(phrase_buffer.size());
  for(uint64_t i = 0; i < phrase_buffer.size(); i++) { this->phrases[i] = phrase_buffer[i]; }
  this->lengths = sd_vector<>(phrase_starts.begin(), phrase_starts.end());
  this->ones = sd_vector<>(phrase_lengths.begin(), phrase_lengths.end());

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

} // namespace sdsl
