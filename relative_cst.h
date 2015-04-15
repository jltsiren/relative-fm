#ifndef _RELATIVE_FM_RELATIVE_CST_H
#define _RELATIVE_FM_RELATIVE_CST_H

#include <sdsl/cst_iterators.hpp>

#include "relative_fm.h"
#include "relative_lcp.h"

namespace relative
{

//------------------------------------------------------------------------------

struct rcst_node
{
  uint64_t sp, ep;
  uint64_t left_lcp, right_lcp; // lcp[sp], lcp[ep + 1]

  rcst_node(uint64_t _sp, uint64_t _ep, uint64_t left, uint64_t right) :
    sp(_sp), ep(_ep), left_lcp(left), right_lcp(right)
  {
  }
};

//------------------------------------------------------------------------------

template<class IndexType = RelativeFM<>>
class RelativeCST
{
public:
  typedef rcst_node node_type;
  typedef uint64_t  size_type;

  typedef cst_dfs_const_forward_iterator<RelativeCST> const_iterator;

//------------------------------------------------------------------------------
  RelativeCST(const IndexType& _index, const RelativeLCP& _lcp) :
    index(_index), lcp(_lcp)
  {
  }

  ~RelativeCST() {}

  inline uint64_t size() const { return this->index.size(); }

  uint64_t reportSize(bool print = false) const
  {
    uint64_t index_bytes = this->index.reportSize();
    uint64_t lcp_bytes = this->lcp.reportSize();
    uint64_t bytes = index_bytes + lcp_bytes;

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
    SDSL-compatible CST operations. Not the full set yet, just enough for const_iterator.

    Tree operations from Fischer2009a pages 7-8 with fixes/optimizations from Canovas2010 pages 8-9.
  */

  inline node_type root() const { return node_type(0, this->size() - 1, 0, 0); }
  inline bool is_leaf(const node_type& v) const { return (v.sp == v.ep); }

  node_type parent(const node_type& v) const
  {
    if(v == this->root()) { return this->root(); }

    uint64_t k = (v.left_lcp > v.right_lcp ? v.sp : v.ep + 1);
    range_type left = this->lcp.psv(k), right = this->lcp.nsv(k);
    if(left.first >= this->size()) { left.first = 0; }  // No psv found.

    return node_type(left.first, right.first - 1, left.second, right.second);
  }

  // i is 1-based.
  node_type select_child(const node_type& v, size_type i) const
  {
    if(this->is_leaf(v) || i == 0) { return this->root(); }

    node_type res = this->first_child(v);
    for(uint64_t j = 1; j < i && res != this->root(); j++) { res = this->sibling(res); }

    return res;
  }

  node_type first_child(const node_type& v) const
  {
    if(this->is_leaf(v) { return this->root(); }

    range_type right = this->rmq(v.sp + 1, v.ep);

    return node_type(v.sp, right.first - 1, v.left_lcp, right.second);
  }

  node_type sibling(const node_type& v) const
  {
    if(v == this->root()) { return this->root(); }
    if(v.left_lcp > v.right_lcp) { return this->root(); } // v is the last child of its parent.

    range_type right = this->nsev(v.ep + 1);

    return node_type(v.sp + 1, right.first - 1, v.right_lcp, right.second);
  }

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
