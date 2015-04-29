/*
  Copyright (c) 2015 Genome Research Ltd.
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

#ifndef _RELATIVE_FM_RLZ_VECTOR_H
#define _RELATIVE_FM_RLZ_VECTOR_H

#include "rlz.h"
#include "support.h"

namespace relative
{

//------------------------------------------------------------------------------

class RLZVector
{
public:
  // The optional bv_fmi allows reusing an already built index for the reference,
  // making parsing much faster.
  RLZVector(const bit_vector& text, const bit_vector& _reference,
    const bit_vector::rank_1_type& _ref_rank,
    const bit_vector::select_1_type& _ref_select_1,
    const bit_vector::select_0_type& _ref_select_0,
    const bv_fmi* fmi = 0);
  RLZVector(std::istream& input, const bit_vector& _reference,
    const bit_vector::rank_1_type& _ref_rank,
    const bit_vector::select_1_type& _ref_select_1,
    const bit_vector::select_0_type& _ref_select_0);
  RLZVector(const RLZVector& v);
  RLZVector(RLZVector&& v);
  ~RLZVector();

  uint64_t reportSize() const;
  uint64_t writeTo(std::ostream& output) const;

  uint64_t size() const { return this->blocks.sum(); }
  uint64_t items() const { return this->ones.sum(); }

  /*
    These follow SDSL conventions.
  */
  uint64_t rank(uint64_t i) const;
  uint64_t select_1(uint64_t i) const;
  uint64_t select_0(uint64_t i) const;
  bool operator[](uint64_t i) const;

private:
  const bit_vector&                 reference;
  const bit_vector::rank_1_type&    ref_rank;
  const bit_vector::select_1_type&  ref_select_1;
  const bit_vector::select_0_type&  ref_select_0;

  relative_encoder                  phrases;

  CumulativeArray                   blocks;
  CumulativeArray                   ones;
  CumulativeArray                   zeros;

  bit_vector                        mismatches;

  // Counts the number of 1-bits in reference[ref_pos, ref_pos + phrase_length - 1].
  inline uint64_t oneBits(uint64_t ref_pos, uint64_t phrase_length) const
  {
    if(phrase_length == 0) { return 0; }
    return this->ref_rank(ref_pos + phrase_length) - this->ref_rank(ref_pos);
  }

  // Returns the offset of the i'th 1-bit in reference[ref_pos, ...] (i is 1-based).
  inline uint64_t findBit(uint64_t ref_pos, uint64_t i) const
  {
    return this->ref_select_1(this->ref_rank(ref_pos) + i) - ref_pos;
  }

  // Returns the offset of the i'th 0-bit in reference[ref_pos, ...] (i is 1-based).
  inline uint64_t findZero(uint64_t ref_pos, uint64_t i) const
  {
    return this->ref_select_0(ref_pos - this->ref_rank(ref_pos) + i) - ref_pos;
  }

  RLZVector();
  RLZVector& operator=(const RLZVector& v);
  RLZVector& operator=(RLZVector&& v);
};

//------------------------------------------------------------------------------

// Forward declarations for rank/select support.
template<uint8_t t_b = 1> class rank_support_rlz;
template<uint8_t t_b = 1> class select_support_rlz;

/*
  A wrapper class that contains either a bit_vector or an RLZVector. Use
  compress()/decompress() to switch between the encodings. Loading a compressed
  bitvector from disk requires the special version of load(). It is recommended
  to decompress a bitvector used as a part of another structure before serializing
  the structure.

  FIXME: Template?
*/
class rlz_vector
{
public:
  typedef bit_vector::size_type                     size_type;
  typedef bool                                      value_type;
  typedef bit_vector::difference_type               difference_type;
  typedef random_access_const_iterator<rlz_vector>  iterator;
  typedef bv_tag                                    index_category;

  typedef rank_support_rlz<1>                       rank_1_type;
  typedef rank_support_rlz<0>                       rank_0_type;
  typedef select_support_rlz<1>                     select_1_type;
  typedef select_support_rlz<0>                     select_0_type;

private:
  size_type                 m_size = 0;
  bit_vector                m_plain;
  bit_vector::rank_1_type   m_plain_rank;
  bit_vector::select_1_type m_plain_select_1;
  bit_vector::select_0_type m_plain_select_0;

  RLZVector*                m_compressed;

  void copy(const rlz_vector& v);
  void init_support();
  void clear_plain();
  void clear_compressed();

public:
  const bit_vector&                 plain           = m_plain;
  const bit_vector::rank_1_type&    plain_rank      = m_plain_rank;
  const bit_vector::select_1_type&  plain_select_1  = m_plain_select_1;
  const bit_vector::select_0_type&  plain_select_0  = m_plain_select_0;

  const RLZVector*                  compressed;

  rlz_vector();
  rlz_vector(const rlz_vector& v);
  rlz_vector(rlz_vector&& v);
  rlz_vector(const bit_vector& v);
  rlz_vector(bit_vector&& v);
  // FIXME construction from iterators.
  ~rlz_vector();

  // Reusing an existing index for the reference makes compression much faster.
  void compress(const bit_vector& reference,
    const bit_vector::rank_1_type& ref_rank,
    const bit_vector::select_1_type& ref_select_1,
    const bit_vector::select_0_type& ref_select_0,
    const bv_fmi* fmi = 0);
  void decompress();  // FIXME implement

  void swap(rlz_vector& v);
  rlz_vector& operator=(const rlz_vector& v);
  rlz_vector& operator=(rlz_vector&& v);

  size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);
  void load(std::istream& in, const bit_vector& reference,
    const bit_vector::rank_1_type& ref_rank,
    const bit_vector::select_1_type& ref_select_1,
    const bit_vector::select_0_type& ref_select_0);

  inline size_type size() const { return this->m_size; }
  inline iterator begin() const { return iterator(this, 0); }
  inline iterator end() const { return iterator(this, this->size()); }
  inline bool isCompressed() const { return (this->m_compressed != 0); }

  inline value_type operator[](size_type i) const
  {
    return (this->isCompressed() ? (*(this->compressed))[i] : this->plain[i]);
  }

  // This only works for non-compressed bitvectors.
  inline uint64_t get_int(size_type idx, const uint8_t len = 64) const
  {
    return this->plain.get_int(idx, len);
  }
};  // class rlz_vector

//------------------------------------------------------------------------------

template<uint8_t t_b>
struct rank_support_rlz_trait
{
  typedef bit_vector::size_type size_type;
  static size_type adjust_rank(size_type r, size_type) { return r; }
};

template<>
struct rank_support_rlz_trait<0>
{
  typedef bit_vector::size_type size_type;
  static size_type adjust_rank(size_type r, size_type n) { return n - r; }
};

template<uint8_t t_b>
class rank_support_rlz
{
  static_assert(t_b == 1u or t_b == 0u , "rank_support_rlz: bit pattern must be `0` or `1`");

public:
  typedef bit_vector::size_type size_type;
  typedef rlz_vector            bit_vector_type;
  enum { bit_pat = t_b };

private:
  const bit_vector_type* m_v;

public:
  explicit rank_support_rlz(const bit_vector_type* v = nullptr)
  {
    this->set_vector(v);
  }

  size_type rank(size_type i) const
  {
    assert(this->m_v != nullptr);
    assert(i <= this->m_v->size());
    uint64_t r = (this->m_v->isCompressed() ? this->m_v->compressed->rank(i) : this->m_v->plain_rank(i));
    return rank_support_rlz_trait<t_b>::adjust_rank(r, this->m_v->size());
  }

  size_type operator()(size_type i) const { return this->rank(i); }
  size_type size() const { return (this->m_v == 0 ? 0 : this->m_v->size()); }

  void set_vector(const bit_vector_type* v = nullptr) { this->m_v = v; }

  rank_support_rlz& operator=(const rank_support_rlz& rs)
  {
    if(this != &rs) { this->set_vector(rs.m_v); }
    return *this;
  }

  void swap(rank_support_rlz& rs)
  {
    if(this != &rs) { std::swap(this->m_v, rs.m_v); }
  }

  size_type serialize(std::ostream&, structure_tree_node* v = nullptr, std::string name = "") const
  {
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    structure_tree::add_size(child, 0);
    return 0;
  }

  void load(const std::istream&, const bit_vector_type* v = nullptr)
  {
    this->set_vector(v);
  }
};

//------------------------------------------------------------------------------

template<uint8_t t_b, class t_rlz_vec>
struct select_support_rlz_trait
{
  typedef bit_vector::size_type size_type;
  static size_type select(size_type i, const t_rlz_vec* v)
  {
    return (v->isCompressed() ? v->compressed->select_1(i) : v->plain_select_1(i));
  }
};

template<class t_rlz_vec>
struct select_support_rlz_trait<0, t_rlz_vec>
{
  typedef bit_vector::size_type size_type;
  static size_type select(size_type i, const t_rlz_vec* v)
  {
    return (v->isCompressed() ? v->compressed->select_0(i) : v->plain_select_0(i));
  }
};

template<uint8_t t_b>
class select_support_rlz
{
  static_assert(t_b == 1u or t_b == 0u , "select_support_rlz: bit pattern must be `0` or `1`");

public:
  typedef bit_vector::size_type size_type;
  typedef rlz_vector            bit_vector_type;
  enum { bit_pat = t_b };

private:
  const bit_vector_type* m_v;

public:
  explicit select_support_rlz(const bit_vector_type* v = nullptr)
  {
    this->set_vector(v);
  }

  size_type select(size_type i) const
  {
    return select_support_rlz_trait<t_b, bit_vector_type>::select(i, m_v);
  }

  size_type operator()(size_type i) const { return this->select(i); }
  size_type size() const { return (this->m_v == 0 ? 0 : this->m_v->size()); }

  void set_vector(const bit_vector_type* v = nullptr) { this->m_v = v; }

  select_support_rlz& operator=(const select_support_rlz& rs)
  {
    if(this != &rs) { this->set_vector(rs.m_v); }
    return *this;
  }

  void swap(select_support_rlz& rs)
  {
    if(this != &rs) { std::swap(this->m_v, rs.m_v); }
  }

  size_type serialize(std::ostream&, structure_tree_node* v = nullptr, std::string name = "") const
  {
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    structure_tree::add_size(child, 0);
    return 0;
  }

  void load(const std::istream&, const bit_vector_type* v = nullptr)
  {
    this->set_vector(v);
  }
};

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_RLZ_VECTOR_H
