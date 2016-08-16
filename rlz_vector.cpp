/*
  Copyright (c) 2015, 2016 Genome Research Ltd.
  Copyright (c) 2014 Jouni Siren

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

#include "rlz.h"
#include "rlz_vector.h"

namespace relative
{

//------------------------------------------------------------------------------

rlz_vector::rlz_vector()
{
  this->m_compressed = 0;
  this->compressed = 0;
}

rlz_vector::rlz_vector(const rlz_vector& v)
{
  this->m_compressed = 0;
  this->compressed = 0;
  this->copy(v);
}

rlz_vector::rlz_vector(rlz_vector&& v)
{
  this->m_compressed = 0;
  this->compressed = 0;
  *this = std::move(v);
}

rlz_vector::rlz_vector(const sdsl::bit_vector& v)
{
  this->m_compressed = 0;
  this->compressed = 0;
  this->m_size = v.size();
  this->m_plain = v;
  this->init_support();
}

rlz_vector::rlz_vector(sdsl::bit_vector&& v)
{
  this->m_compressed = 0;
  this->compressed = 0;
  this->m_size = v.size();
  this->m_plain = std::move(v);
  this->init_support();
}

rlz_vector::~rlz_vector()
{
  this->clear_compressed();
}

void
rlz_vector::compress(const sdsl::bit_vector& reference,
  const sdsl::bit_vector::rank_1_type& ref_rank,
  const sdsl::bit_vector::select_1_type& ref_select_1,
  const sdsl::bit_vector::select_0_type& ref_select_0,
  const bv_fmi* fmi)
{
  if(this->isCompressed())  // FIXME maybe decompress and recompress
  {
    std::cerr << "rlz_vector::compress(): The vector is already compressed" << std::endl;
    return;
  }

  this->m_compressed = new RLZVector(this->plain, reference, ref_rank, ref_select_1, ref_select_0, fmi);
  this->compressed = this->m_compressed;
  this->clear_plain();
}

void
rlz_vector::copy(const rlz_vector& v)
{
  this->m_size = v.m_size;
  this->clear_compressed();

  if(v.isCompressed())
  {
    this->clear_plain();
    this->m_compressed = new RLZVector(*(v.compressed));
    this->compressed = this->m_compressed;
  }
  else
  {
    this->m_plain = v.m_plain;
    this->m_plain_rank = v.m_plain_rank; this->m_plain_rank.set_vector(&(this->m_plain));
    this->m_plain_select_1 = v.m_plain_select_1; this->m_plain_select_1.set_vector(&(this->m_plain));
    this->m_plain_select_0 = v.m_plain_select_0; this->m_plain_select_0.set_vector(&(this->m_plain));
  }
}

void
rlz_vector::init_support()
{
  sdsl::util::init_support(this->m_plain_rank, &(this->m_plain));
  sdsl::util::init_support(this->m_plain_select_1, &(this->m_plain));
  sdsl::util::init_support(this->m_plain_select_0, &(this->m_plain));
}

void
rlz_vector::clear_plain()
{
  sdsl::util::clear(this->m_plain);
  sdsl::util::clear(this->m_plain_rank);
  sdsl::util::clear(this->m_plain_select_1);
  sdsl::util::clear(this->m_plain_select_0);
}

void
rlz_vector::clear_compressed()
{
  delete this->m_compressed; this->m_compressed = 0;
  this->compressed = 0;
}


void
rlz_vector::swap(rlz_vector& v)
{
  if(this != &v)
  {
    std::swap(this->m_size, v.m_size);
    this->m_plain.swap(v.m_plain);
    sdsl::util::swap_support(this->m_plain_rank, v.m_plain_rank, &(this->m_plain), &(v.m_plain));
    sdsl::util::swap_support(this->m_plain_select_1, v.m_plain_select_1, &(this->m_plain), &(v.m_plain));
    sdsl::util::swap_support(this->m_plain_select_0, v.m_plain_select_0, &(this->m_plain), &(v.m_plain));
    std::swap(this->m_compressed, v.m_compressed);
    std::swap(this->compressed, v.compressed);
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

    this->clear_compressed();
    this->m_compressed = v.m_compressed; v.m_compressed = 0;
    this->compressed = v.compressed; v.compressed = 0;
  }
  return *this;
}

rlz_vector::size_type
rlz_vector::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += sdsl::write_member(this->m_size, out, child, "size");
  bool compressed = this->isCompressed();
  written_bytes += sdsl::write_member(compressed, out, child, "compressed");

  if(compressed)
  {
    written_bytes += this->compressed->writeTo(out);
  }
  else
  {
    written_bytes += this->m_plain.serialize(out, child, "plain");
    written_bytes += this->m_plain_rank.serialize(out, child, "plain_rank");
    written_bytes += this->m_plain_select_1.serialize(out, child, "plain_select_1");
    written_bytes += this->m_plain_select_0.serialize(out, child, "plain_select_0");
  }

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
rlz_vector::load(std::istream& in)
{
  this->clear_plain();
  this->clear_compressed();

  size_type new_size = 0;
  sdsl::read_member(new_size, in);
  bool compressed = false;
  sdsl::read_member(compressed, in);

  if(compressed)
  {
    std::cerr << "rlz_vector::load(): Cannot load a compressed bitvector without a reference" << std::endl;
    return;
  }
  else
  {
    this->m_size = new_size;
    this->m_plain.load(in);
    this->m_plain_rank.load(in, &(this->m_plain));
    this->m_plain_select_1.load(in, &(this->m_plain));
    this->m_plain_select_0.load(in, &(this->m_plain));
  }
}

void
rlz_vector::load(std::istream& in, const sdsl::bit_vector& reference,
    const sdsl::bit_vector::rank_1_type& ref_rank,
    const sdsl::bit_vector::select_1_type& ref_select_1,
    const sdsl::bit_vector::select_0_type& ref_select_0)
{
  this->clear_plain();
  this->clear_compressed();

  sdsl::read_member(this->m_size, in);
  bool compressed = false;
  sdsl::read_member(compressed, in);

  if(compressed)
  {
    this->m_compressed = new RLZVector(in, reference, ref_rank, ref_select_1, ref_select_0);
    this->compressed = this->m_compressed;
  }
  else  // We expected a compressed bitvector, but it's not dangerous.
  {
    this->m_plain.load(in);
    this->m_plain_rank.load(in, &(this->m_plain));
    this->m_plain_select_1.load(in, &(this->m_plain));
    this->m_plain_select_0.load(in, &(this->m_plain));
  }
}

//------------------------------------------------------------------------------

RLZVector::RLZVector(const sdsl::bit_vector& text, const sdsl::bit_vector& _reference,
  const sdsl::bit_vector::rank_1_type& _ref_rank,
  const sdsl::bit_vector::select_1_type& _ref_select_1,
  const sdsl::bit_vector::select_0_type& _ref_select_0,
  const bv_fmi* fmi) :
  reference(_reference), ref_rank(_ref_rank), ref_select_1(_ref_select_1), ref_select_0(_ref_select_0)
{
  std::vector<size_type> phrase_starts, phrase_lengths;
  if(fmi != 0)
  {
    relativeLZSuccinct(text, *fmi, phrase_starts, phrase_lengths, this->mismatches);
  }
  else
  {
    relativeLZSuccinct(text, this->reference, phrase_starts, phrase_lengths, this->mismatches);
  }

  // Initialize phrases and blocks.
  this->phrases.init(phrase_starts, phrase_lengths);
  sdsl::util::assign(this->blocks, CumulativeArray(phrase_lengths));

  // Initialize ones and zeros, using phrase_starts for ones and phrase_lengths for zeros.
  for(size_type i = 0; i < phrase_starts.size(); i++)
  {
    phrase_starts[i] = this->oneBits(phrase_starts[i], phrase_lengths[i] - 1);
    if(this->mismatches[i]) { phrase_starts[i]++; }
    phrase_lengths[i] -= phrase_starts[i];
  }
  sdsl::util::assign(this->ones, CumulativeArray(phrase_starts));
  sdsl::util::assign(this->zeros, CumulativeArray(phrase_lengths));
}

RLZVector::RLZVector(std::istream& input, const sdsl::bit_vector& _reference,
  const sdsl::bit_vector::rank_1_type& _ref_rank,
  const sdsl::bit_vector::select_1_type& _ref_select_1,
  const sdsl::bit_vector::select_0_type& _ref_select_0) :
  reference(_reference), ref_rank(_ref_rank), ref_select_1(_ref_select_1), ref_select_0(_ref_select_0)
{
  this->phrases.load(input);
  this->blocks.load(input);
  this->ones.load(input);
  this->zeros.load(input);
  this->mismatches.load(input);
}

RLZVector::RLZVector(const RLZVector& v) :
  reference(v.reference), ref_rank(v.ref_rank), ref_select_1(v.ref_select_1), ref_select_0(v.ref_select_0)
{
  this->phrases = v.phrases;

  this->blocks = v.blocks;
  this->ones = v.ones;
  this->zeros = v.zeros;

  this->mismatches = v.mismatches;
}

RLZVector::RLZVector(RLZVector&& v) :
  reference(v.reference), ref_rank(v.ref_rank), ref_select_1(v.ref_select_1), ref_select_0(v.ref_select_0)
{
  this->phrases = std::move(v.phrases);

  this->blocks = std::move(v.blocks);
  this->ones = std::move(v.ones);
  this->zeros = std::move(v.zeros);

  this->mismatches = std::move(v.mismatches);
}

RLZVector::~RLZVector()
{
}

//------------------------------------------------------------------------------

size_type
RLZVector::reportSize() const
{
  size_type bytes = sizeof(*this);
  bytes += this->phrases.reportSize();
  bytes += sdsl::size_in_bytes(this->blocks);
  bytes += sdsl::size_in_bytes(this->ones);
  bytes += sdsl::size_in_bytes(this->zeros);
  bytes += sdsl::size_in_bytes(this->mismatches);
  return bytes;
}


size_type
RLZVector::writeTo(std::ostream& output) const
{
  size_type bytes = 0;
  bytes += this->phrases.serialize(output);
  bytes += this->blocks.serialize(output);
  bytes += this->ones.serialize(output);
  bytes += this->zeros.serialize(output);
  bytes += this->mismatches.serialize(output);
  return bytes;
}

//------------------------------------------------------------------------------

size_type
RLZVector::rank(size_type i) const
{
  if(i >= this->size()) { return this->items(); }

  size_type phrase = this->blocks.inverse(i); // Phrase is 0-based.
  if(phrase == 0) { return this->oneBits(this->phrases.decode(0, 0), i); }
  size_type text_pos = this->blocks.sum(phrase); // Starting position of the phrase in the text.

  return this->ones.sum(phrase) + this->oneBits(this->phrases.decode(phrase, text_pos), i - text_pos);
}

size_type
RLZVector::select_1(size_type i) const
{
  // FIXME what happens for i == 0?
  if(i > this->items()) { return this->size(); }

  size_type phrase = this->ones.inverse(i - 1); // Phrase is 0-based.
  // Check if the requested 1-bit is the mismatching bit at the end of the phrase.
  if(this->ones.isLast(i - 1) && this->mismatches[phrase]) { return this->blocks.sum(phrase + 1) - 1; }
  if(phrase == 0) { return this->findBit(this->phrases.decode(0, 0), i); }
  i -= this->ones.sum(phrase); // i is now relative to the phrase.
  size_type text_pos = this->blocks.sum(phrase);  // Starting position of the phrase in the text.

  return text_pos + this->findBit(this->phrases.decode(phrase, text_pos), i);
}

size_type
RLZVector::select_0(size_type i) const
{
  // FIXME what happens for i == 0?
  if(i > this->size() - this->items()) { return this->size(); }

  size_type phrase = this->zeros.inverse(i - 1);  // Phrase is 0-based.
  // Check if the requested 0-bit is the mismatching bit at the end of the phrase.
  if(this->zeros.isLast(i - 1) && !(this->mismatches[phrase])) { return this->blocks.sum(phrase + 1) - 1; }
  if(phrase == 0) { return this->findZero(this->phrases.decode(0, 0), i); }
  i -= this->zeros.sum(phrase);  // i is now relative to the phrase.
  size_type text_pos = this->blocks.sum(phrase);  // Starting position of the phrase in the text.

  return text_pos + this->findZero(this->phrases.decode(phrase, text_pos), i);
}

bool
RLZVector::operator[](size_type i) const
{
  if(i >= this->size()) { return false; }

  size_type phrase = this->blocks.inverse(i); // Phrase is 0-based.
  if(this->blocks.isLast(i)) { return this->mismatches[phrase]; } // The mismatching bit.

  return this->reference[this->phrases.decode(phrase, i)];
}

//------------------------------------------------------------------------------

} // namespace relative
