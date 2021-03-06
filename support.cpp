/*
  Copyright (c) 2015, 2016, 2017 Genome Research Ltd.
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

#include "support.h"

namespace relative
{

//------------------------------------------------------------------------------

Alphabet::Alphabet() :
  char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma)
{
  this->m_sigma = 0;
}

Alphabet::Alphabet(const Alphabet& a) :
  char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma)
{
  this->copy(a);
}

Alphabet::Alphabet(Alphabet&& a) :
  char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma)
{
  *this = std::move(a);
}

Alphabet::~Alphabet()
{
}

void
Alphabet::copy(const Alphabet& a)
{
  this->m_char2comp = a.m_char2comp;
  this->m_comp2char = a.m_comp2char;
  this->m_C = a.m_C;
  this->m_sigma = a.m_sigma;
}

void
Alphabet::swap(Alphabet& a)
{
  if(this != &a)
  {
    this->m_char2comp.swap(a.m_char2comp);
    this->m_comp2char.swap(a.m_comp2char);
    this->m_C.swap(a.m_C);
    std::swap(this->m_sigma, a.m_sigma);
  }
}

Alphabet&
Alphabet::operator=(const Alphabet& a)
{
  if(this != &a) { this->copy(a); }
  return *this;
}

Alphabet&
Alphabet::operator=(Alphabet&& a)
{
  if(this != &a)
  {
    this->m_char2comp = std::move(a.m_char2comp);
    this->m_comp2char = std::move(a.m_comp2char);
    this->m_C = std::move(a.m_C);
    this->m_sigma = a.m_sigma;
  }
  return *this;
}

Alphabet::size_type
Alphabet::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += this->m_char2comp.serialize(out, child, "char2comp");
  written_bytes += this->m_comp2char.serialize(out, child, "comp2char");
  written_bytes += this->m_C.serialize(out, child, "C");
  written_bytes += sdsl::write_member(this->m_sigma, out, child, "sigma");
  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
Alphabet::load(std::istream& in)
{
  this->m_char2comp.load(in);
  this->m_comp2char.load(in);
  this->m_C.load(in);
  sdsl::read_member(this->m_sigma, in);
}

bool
Alphabet::assign(const std::string& alphabet_string)
{
  if(alphabet_string.length() != this->m_sigma)
  {
    std::cerr << "Alphabet::assign(): Alphabet string has length " << alphabet_string.length()
              << " (should be " << this->m_sigma << ")" << std::endl;
    return false;
  }

  sdsl::bit_vector exists(MAX_SIGMA);
  sdsl::util::assign(this->m_char2comp, sdsl::int_vector<8>(MAX_SIGMA, 0));
  for(size_type i = 0; i < alphabet_string.length(); i++)
  {
    char_type c = alphabet_string[i];
    if(exists[c])
    {
      std::cerr << "Alphabet::assign(): The alphabet string contains multiple occurrences of " << c << std::endl;
      return false;
    }
    exists[c] = 1;
    this->m_char2comp[c] = i;
    this->m_comp2char[i] = c;
  }
  return true;
}

//------------------------------------------------------------------------------

LCS::LCS()
{
  this->lcs_size = 0;
}

LCS::LCS(const sdsl::bit_vector& a, const sdsl::bit_vector& b, LCS::size_type _lcs_size)
{
  this->ref = a; this->seq = b;
  this->lcs_size = _lcs_size;

  sdsl::util::init_support(this->ref_rank, &(this->ref));
#ifndef USE_HYBRID_BITVECTORS
  sdsl::util::init_support(this->ref_select, &(this->ref));
#endif

  sdsl::util::init_support(this->seq_rank, &(this->seq));
#ifndef USE_HYBRID_BITVECTORS
  sdsl::util::init_support(this->seq_select, &(this->seq));
#endif
}

LCS::LCS(const LCS& l)
{
  this->copy(l);
}

LCS::LCS(LCS&& l)
{
  *this = std::move(l);
}

LCS::~LCS()
{
}

void
LCS::copy(const LCS& l)
{
  this->ref = l.ref;
  this->ref_rank = l.ref_rank;
#ifndef USE_HYBRID_BITVECTORS
  this->ref_select = l.ref_select;
#endif

  this->seq = l.seq;
  this->seq_rank = l.seq_rank;
#ifndef USE_HYBRID_BITVECTORS
  this->seq_select = l.seq_select;
#endif

  this->lcs_size = l.lcs_size;

  this->set_vectors();
}

void
LCS::swap(LCS& l)
{
  if(this != &l)
  {
    this->ref.swap(l.ref);
    sdsl::util::swap_support(this->ref_rank, l.ref_rank, &(this->ref), &(l.ref));
#ifndef USE_HYBRID_BITVECTORS
    sdsl::util::swap_support(this->ref_select, l.ref_select, &(this->ref), &(l.ref));
#endif

    this->seq.swap(l.seq);
    sdsl::util::swap_support(this->seq_rank, l.seq_rank, &(this->seq), &(l.seq));
#ifndef USE_HYBRID_BITVECTORS
    sdsl::util::swap_support(this->seq_select, l.seq_select, &(this->seq), &(l.seq));
#endif

    std::swap(this->lcs_size, l.lcs_size);
  }
}

LCS&
LCS::operator=(const LCS& l)
{
  if(this != &l) { this->copy(l); }
  return *this;
}

LCS&
LCS::operator=(LCS&& l)
{
  if(this != &l)
  {
    this->ref = std::move(l.ref);
    this->ref_rank = std::move(l.ref_rank);
#ifndef USE_HYBRID_BITVECTORS
    this->ref_select = std::move(l.ref_select);
#endif

    this->seq = std::move(l.seq);
    this->seq_rank = std::move(l.seq_rank);
#ifndef USE_HYBRID_BITVECTORS
    this->seq_select = std::move(l.seq_select);
#endif

    this->lcs_size = l.lcs_size;

    this->set_vectors();
  }
  return *this;
}

size_type
LCS::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->ref.serialize(out, child, "ref");
  written_bytes += this->ref_rank.serialize(out, child, "ref_rank");
#ifndef USE_HYBRID_BITVECTORS
  written_bytes += this->ref_select.serialize(out, child, "ref_select");
#endif

  written_bytes += this->seq.serialize(out, child, "seq");
  written_bytes += this->seq_rank.serialize(out, child, "seq_rank");
#ifndef USE_HYBRID_BITVECTORS
  written_bytes += this->seq_select.serialize(out, child, "seq_select");
#endif

  written_bytes += sdsl::write_member(this->lcs_size, out, child, "lcs_size");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
LCS::load(std::istream& in)
{
  this->ref.load(in);
  this->ref_rank.load(in, &(this->ref));
#ifndef USE_HYBRID_BITVECTORS
  this->ref_select.load(in, &(this->ref));
#endif

  this->seq.load(in);
  this->seq_rank.load(in, &(this->seq));
#ifndef USE_HYBRID_BITVECTORS
  this->seq_select.load(in, &(this->seq));
#endif

  sdsl::read_member(this->lcs_size, in);
}

void
LCS::set_vectors()
{
  this->ref_rank.set_vector(&(this->ref));
#ifndef USE_HYBRID_BITVECTORS
  this->ref_select.set_vector(&(this->ref));
#endif

  this->seq_rank.set_vector(&(this->seq));
#ifndef USE_HYBRID_BITVECTORS
  this->seq_select.set_vector(&(this->seq));
#endif
}

#ifdef USE_HYBRID_BITVECTORS
size_type
LCS::select(const LCS::vector_type& vec, const LCS::vector_type::rank_1_type& rank, size_type i) const
{
  if(i == 0) { return 0; }

  // Find the last position, where rank < i.
  size_type low = 0, high = vec.size();
  while(low < high)
  {
    size_type mid = low + (high - low + 1) / 2;
    if(rank(mid) >= i) { high = mid - 1; }
    else { low = mid; }
  }

  return low;
}
#endif

//------------------------------------------------------------------------------

SLArray::SLArray() :
  small(sdsl::int_vector<8>(1, 0))
{
}

SLArray::SLArray(const SLArray& s)
{
  this->copy(s);
}

SLArray::SLArray(SLArray&& s)
{
  *this = std::move(s);
}

SLArray::~SLArray()
{
}

SLArray::SLArray(sdsl::int_vector_buffer<0>& source)
{
  this->buildFrom(source);
}

void
SLArray::copy(const SLArray& s)
{
  this->small = s.small;
  this->large = s.large;
  this->samples = s.samples;
}

void
SLArray::swap(SLArray& s)
{
  if(this != &s)
  {
    this->small.swap(s.small);
    this->large.swap(s.large);
    this->samples.swap(s.samples);
  }
}

SLArray&
SLArray::operator=(const SLArray& s)
{
  if(this != &s) { this->copy(s); }
  return *this;
}

SLArray&
SLArray::operator=(SLArray&& s)
{
  if(this != &s)
  {
    this->small = std::move(s.small);
    this->large = std::move(s.large);
    this->samples = std::move(s.samples);
  }
  return *this;
}

size_type
SLArray::serialize(std::ostream& out, sdsl::structure_tree_node* s, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->small.serialize(out, child, "small");
  written_bytes += this->large.serialize(out, child, "large");
  written_bytes += this->samples.serialize(out, child, "samples");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
SLArray::load(std::istream& in)
{
  this->small.load(in);
  this->large.load(in);
  this->samples.load(in);
}

sdsl::int_vector<64>
SLArray::extract(size_type from, size_type to) const
{
  from = std::min(from, this->size()); to = std::min(to, this->size());

  sdsl::int_vector<64> result(to - from, 0);
  iterator iter(this, from);
  for(size_type i = 0; i < result.size(); i++, ++iter) { result[i] = *iter; }

  return result;
}

//------------------------------------------------------------------------------

} // namespace relative
