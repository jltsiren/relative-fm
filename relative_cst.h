/*
  Copyright (c) 2015 Genome Research Ltd.

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

#ifndef _RELATIVE_FM_RELATIVE_CST_H
#define _RELATIVE_FM_RELATIVE_CST_H

#include <sdsl/suffix_trees.hpp>

#include "relative_fm.h"
#include "relative_lcp.h"

namespace relative
{

//------------------------------------------------------------------------------

struct rcst_node
{
  size_type sp, ep;
  size_type left_lcp, right_lcp; // lcp[sp], lcp[ep + 1]

  rcst_node() : sp(0), ep(0), left_lcp(0), right_lcp(0) {}

  rcst_node(size_type _sp, size_type _ep, size_type left, size_type right) :
    sp(_sp), ep(_ep), left_lcp(left), right_lcp(right)
  {
  }

  inline bool operator== (const rcst_node& node) const
  {
    return (this->sp == node.sp && this->ep == node.ep);
  }

  inline bool operator!= (const rcst_node& node) const
  {
    return (this->sp != node.sp || this->ep != node.ep);
  }
};

std::ostream&
operator<<(std::ostream& stream, const rcst_node& node)
{
  return stream << range_type(node.sp, node.ep);
}

//------------------------------------------------------------------------------

template<class IndexType = RelativeFM<>>
class RelativeCST
{
public:
  typedef rcst_node           node_type;
  typedef relative::size_type size_type;
  typedef relative::char_type char_type;

  typedef sdsl::cst_dfs_const_forward_iterator<RelativeCST> const_iterator;

//------------------------------------------------------------------------------
  RelativeCST(const IndexType& _index, const RelativeLCP& _lcp) :
    index(_index), lcp(_lcp)
  {
  }

  ~RelativeCST() {}

  inline size_type size() const { return this->index.size(); }

  size_type reportSize(bool print = false) const
  {
    size_type index_bytes = this->index.reportSize();
    size_type lcp_bytes = this->lcp.reportSize();
    size_type bytes = index_bytes + lcp_bytes;

    if(print)
    {
      printSize("FM-index", index_bytes, this->size());
      printSize("Relative LCP", lcp_bytes, this->size());
      printSize("Relative CST", bytes, this->size());
      std::cout << std::endl;
    }

    return bytes;
  }

  const IndexType&   index;
  const RelativeLCP& lcp;

//------------------------------------------------------------------------------

  /*
    This uses forward searching, which is slower than backward searching using index.find().
  */
  template<class Iter>
  range_type find(Iter from, Iter to) const
  {
    node_type next = this->root();
    size_type depth = 0, next_depth = 0;
    size_type bwt_pos = this->size();

    while(from != to)
    {
      if(!(this->forward_search(next, next_depth, depth, *from, bwt_pos))) { return range_type(1, 0); }
      ++from; depth++;
    }

    return range_type(next.sp, next.ep);
  }

  /*
    One step of forward searching. Returns false if not found. To start from the root, set
    next = root(), next_depth = 0, and bwt_pos = size().
  */
  bool
  forward_search(node_type& next, size_type& next_depth, size_type depth, char_type c, size_type& bwt_pos) const
  {
    char_type comp = this->index.alpha.char2comp[c];
    if(depth >= next_depth) // Next node reached, follow a new edge.
    {
      node_type temp = this->child(next, c, bwt_pos);
      if(temp == this->root()) { return false; }
      next = temp; next_depth = this->depth(next);
    }
    else  // Continue in the edge.
    {
      bwt_pos = this->index.Psi(bwt_pos);
      if(bwt_pos < this->index.alpha.C[comp] || bwt_pos >= this->index.alpha.C[comp + 1])
      {
        return false;
      }
    }
    return true;
  }

//------------------------------------------------------------------------------

  /*
    SDSL-compatible CST operations. Not the full set yet.

    Tree operations from Fischer2009a pages 7-8 with fixes/optimizations from Abeliuk2013 pages 14-15.
  */

  inline node_type root() const { return node_type(0, this->size() - 1, 0, 0); }
  inline bool is_leaf(const node_type& v) const { return (v.sp == v.ep); }

  node_type parent(const node_type& v) const
  {
    if(v == this->root()) { return this->root(); }

    size_type k = (v.left_lcp > v.right_lcp ? v.sp : v.ep + 1);
    range_type left = this->lcp.psv(k), right = this->lcp.nsv(k);
    if(left.first >= this->size()) { left.first = 0; }  // No psv found.

    return node_type(left.first, right.first - 1, left.second, right.second);
  }

  // i is 1-based.
  node_type select_child(const node_type& v, size_type i) const
  {
    if(this->is_leaf(v) || i == 0) { return this->root(); }

    node_type res = this->first_child(v);
    for(size_type j = 1; j < i && res != this->root(); j++) { res = this->sibling(res); }

    return res;
  }

  node_type first_child(const node_type& v) const
  {
    if(this->is_leaf(v)) { return this->root(); }

    range_type right = this->lcp.rmq(v.sp + 1, v.ep);

    return node_type(v.sp, right.first - 1, v.left_lcp, right.second);
  }

  node_type sibling(const node_type& v) const
  {
    if(v.ep + 1 >= this->size()) { return this->root(); }
    if(v.left_lcp > v.right_lcp) { return this->root(); } // v is the last child of its parent.

    range_type right = this->lcp.nsev(v.ep + 1);

    return node_type(v.ep + 1, right.first - 1, v.right_lcp, right.second);
  }

  node_type child(const node_type& v, char_type c, size_type& bwt_pos) const
  {
    if(!hasChar(this->index.alpha, c)) { return this->root(); }
    c = this->index.alpha.char2comp[c];

    size_type comp = 0;  // The next comp value to check.
    for(node_type curr = this->first_child(v); curr != this->root(); curr = this->sibling(curr))
    {
      size_type d = std::max(curr.left_lcp, curr.right_lcp);
      bwt_pos = this->index.Psi(curr.sp, d);
      if(bwt_pos >= this->size()) { continue; }

      while(this->index.alpha.C[comp + 1] <= bwt_pos) { comp++; }
      if(comp == c) { return curr; }
      else if(comp > c) { break; }
    }

    return this->root();
  }

  inline node_type child(const node_type& v, char_type c) const
  {
    size_type bwt_pos;
    return this->child(v, c, bwt_pos);
  }

  // i is 1-based; no sanity checking.
  node_type select_leaf(size_type i) const
  {
    i--;
    return node_type(i, i, this->lcp[i], this->lcp[i + 1]);
  }

  // Suffix link.
  node_type sl(const node_type& v) const
  {
    if(v == this->root()) { return this->root(); }
    if(this->is_leaf(v))
    {
      if(v.sp < this->index.sequences()) { return this->root(); }  // v was an empty suffix.
      else { return this->select_leaf(this->index.Psi(v.sp) + 1); }
    }

    size_type sp = this->index.Psi(v.sp), ep = this->index.Psi(v.ep);
    size_type k = this->lcp.rmq(sp + 1, ep).first;
    range_type left = this->lcp.psv(k), right = this->lcp.nsv(k);

    return node_type(left.first, right.first - 1, left.second, right.second);
  }

  size_type depth(const node_type& v) const
  {
    if(this->is_leaf(v)) { return this->size() - this->index.locate(v.sp) - 1; }
    else if(v == this->root()) { return 0; }
    else { return this->lcp.rmq(v.sp + 1, v.ep).second; }
  }

  inline size_type lb(const node_type& v) const { return v.sp; }
  inline size_type rb(const node_type& v) const { return v.ep; }

//------------------------------------------------------------------------------

  const_iterator begin() const
  {
    if(this->size() == 0) { return this->end(); }
    return const_iterator(this, this->root(), false, true);
  }

  const_iterator begin(const node_type& v) const
  {
    if(this->size() == 0 && v == this->root()) { return this->end(); }
    return const_iterator(this, v, false, true);
  }

  const_iterator end() const
  {
    return const_iterator(this, this->root(), true, false);
  }

  const_iterator end(const node_type& v) const
  {
    if(v == this->root()) { return this->end(); }
    return ++const_iterator(this, v, true, true);
  }

//------------------------------------------------------------------------------

private:
  RelativeCST();
  RelativeCST(const RelativeCST&);
  RelativeCST(RelativeCST&&);
  RelativeCST& operator=(const RelativeCST&);
  RelativeCST& operator==(RelativeCST&);
};

//------------------------------------------------------------------------------

} // namespace relative

#endif  // _RELATIVE_FM_RELATIVE_CST_H
