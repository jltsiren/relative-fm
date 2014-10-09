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
  const bit_vector::select_1_type& _ref_select_1,
  const bit_vector::select_0_type& _ref_select_0) :
  reference(_reference), ref_rank(_ref_rank), ref_select_1(_ref_select_1), ref_select_0(_ref_select_0)
{
  std::vector<uint64_t> phrase_starts, phrase_lengths;
  relativeLZSuccinct(text, this->reference, phrase_starts, phrase_lengths, this->mismatches);

  std::vector<uint64_t> phrase_buffer;
  this->phrase_rle.resize(phrase_starts.size()); util::set_to_value(this->phrase_rle, 0);

  uint64_t bits = 0, prev = ~(uint64_t)0, max_val = 0;
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
    phrase_starts[phrase] = this->oneBits(phrase_starts[phrase], phrase_lengths[phrase] - 1);
    if(this->mismatches[phrase]) { phrase_starts[phrase]++; }
  }
  this->phrases.width(bits::hi(max_val) + 1); this->phrases.resize(phrase_buffer.size());
  for(uint64_t i = 0; i < phrase_buffer.size(); i++) { this->phrases[i] = phrase_buffer[i]; }
  this->blocks.init(phrase_lengths);
  this->ones.init(phrase_starts);
  for(uint64_t i = 0; i < phrase_starts.size(); i++) { phrase_lengths[i] -= phrase_starts[i]; }
  this->zeros.init(phrase_lengths);

  util::init_support(this->phrase_rank, &(this->phrase_rle));
}

RLZVector::RLZVector(std::ifstream& input, const bit_vector& _reference,
  const bit_vector::rank_1_type& _ref_rank,
  const bit_vector::select_1_type& _ref_select_1,
  const bit_vector::select_0_type& _ref_select_0) :
  reference(_reference), ref_rank(_ref_rank), ref_select_1(_ref_select_1), ref_select_0(_ref_select_0)
{
  this->phrases.load(input);
  this->phrase_rle.load(input);
  this->blocks.load(input);
  this->ones.load(input);
  this->zeros.load(input);
  this->mismatches.load(input);

  util::init_support(this->phrase_rank, &(this->phrase_rle));
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
  bytes += this->blocks.reportSize();
  bytes += this->ones.reportSize();
  bytes += this->zeros.reportSize();
  bytes += size_in_bytes(this->mismatches);
  return bytes;
}


void
RLZVector::writeTo(std::ofstream& output) const
{
  this->phrases.serialize(output);
  this->phrase_rle.serialize(output);
  this->blocks.writeTo(output);
  this->ones.writeTo(output);
  this->zeros.writeTo(output);
  this->mismatches.serialize(output);
}

//------------------------------------------------------------------------------

uint64_t
RLZVector::rank(uint64_t i) const
{
  if(i >= this->size()) { return this->items(); }

  uint64_t phrase = this->blocks.blockFor(i); // Phrase is 0-based.
  if(phrase == 0) { return this->oneBits(this->refPos(0, 0), i); }
  uint64_t text_pos = this->blocks.itemsAfter(phrase - 1); // Starting position of the phrase in the text.

  return this->ones.itemsAfter(phrase - 1) +
    this->oneBits(this->refPos(phrase, text_pos), i - text_pos);
}

uint64_t
RLZVector::select_1(uint64_t i) const
{
  // FIXME what happens for i == 0?
  if(i > this->items()) { return this->size(); }

  uint64_t phrase = this->ones.blockFor(i - 1); // Phrase is 0-based.
  // Check if the requested 1-bit is the mismatching bit at the end of the phrase.
  if(this->ones.isLast(i - 1) && this->mismatches[phrase]) { return this->blocks.itemsAfter(phrase) - 1; }
  if(phrase == 0) { return this->findBit(this->refPos(0, 0), i); }
  i -= this->ones.itemsAfter(phrase - 1); // i is now relative to the phrase.
  uint64_t text_pos = this->blocks.itemsAfter(phrase - 1);  // Starting position of the phrase in the text.

  return text_pos + this->findBit(this->refPos(phrase, text_pos), i);
}

uint64_t
RLZVector::select_0(uint64_t i) const
{
  // FIXME what happens for i == 0?
  if(i > this->size() - this->items()) { return this->size(); }

  // FIXME find the correct phrase
  uint64_t phrase = this->zeros.blockFor(i - 1);  // Phrase is 0-based.
  // Check if the requested 0-bit is the mismatching bit at the end of the phrase.
  if(this->zeros.isLast(i - 1) && !(this->mismatches[phrase])) { return this->blocks.itemsAfter(phrase) - 1; }
  if(phrase == 0) { return this->findZero(this->refPos(0, 0), i); }
  i -= this->zeros.itemsAfter(phrase - 1);  // i is now relative to the phrase.
  uint64_t text_pos = this->blocks.itemsAfter(phrase - 1);  // Starting position of the phrase in the text.

  return text_pos + this->findZero(this->refPos(phrase, text_pos), i);
}

bool
RLZVector::operator[](uint64_t i) const
{
  if(i >= this->size()) { return false; }

  uint64_t phrase = this->blocks.blockFor(i); // Phrase is 0-based.
  if(this->blocks.isLast(i)) { return this->mismatches[phrase]; } // The mismatching bit.

  return this->reference[this->refPos(phrase, i)];
}

//------------------------------------------------------------------------------

void
rlz_helper::init(const std::vector<uint64_t>& values)
{
  if(values.size() == 0) { return; }

  this->nonzero.resize(values.size()); util::set_to_value(this->nonzero, 1);
  uint64_t zeros = 0;
  for(uint64_t i = 0; i < values.size(); i++)
  {
    if(values[i] == 0) { this->nonzero[i] = 0; zeros++; }
  }
  if(zeros > 0)
  {
    util::init_support(this->nz_rank, &(this->nonzero));
    util::init_support(this->nz_select, &(this->nonzero));
  }
  else
  {
    util::clear(this->nonzero);
    util::clear(this->nz_rank);
    util::clear(this->nz_select);
  }

  std::vector<uint64_t> real_values(values.size() - zeros);
  for(uint64_t i = 0, j = 0, sum = 0; i < values.size(); i++)
  {
    sum += values[i];
    if(values[i] != 0) { real_values[j] = sum - 1; j++; }
  }
  this->v = sd_vector<>(real_values.begin(), real_values.end());
  util::init_support(this->v_rank, &(this->v));
  util::init_support(this->v_select, &(this->v));
}

uint64_t
rlz_helper::reportSize() const
{
  uint64_t bytes = 0;
  bytes += size_in_bytes(this->v) + size_in_bytes(this->v_rank) + size_in_bytes(this->v_select);
  bytes += size_in_bytes(this->nonzero) + size_in_bytes(this->nz_rank) + size_in_bytes(this->nz_select);
  return bytes;
}

void
rlz_helper::load(std::istream& input)
{
  this->v.load(input);
  util::init_support(this->v_rank, &(this->v));
  util::init_support(this->v_select, &(this->v));

  this->nonzero.load(input);
  if(this->nonzero.size() > 0)
  {
    util::init_support(this->nz_rank, &(this->nonzero));
    util::init_support(this->nz_select, &(this->nonzero));
  }
  else
  {
    util::clear(this->nz_rank);
    util::clear(this->nz_select);
  }
}

void
rlz_helper::writeTo(std::ostream& output) const
{
  this->v.serialize(output);
  this->nonzero.serialize(output);
}

//------------------------------------------------------------------------------

} // namespace sdsl
