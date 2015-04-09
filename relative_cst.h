#ifndef _RELATIVE_FM_RELATIVE_CST_H
#define _RELATIVE_FM_RELATIVE_CST_H

#include "relative_fm.h"
#include "relative_lcp.h"

namespace relative
{

//------------------------------------------------------------------------------

template<class IndexType = RelativeFM<>>
class RelativeCST
{
public:
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
    Tree operations from Fischer2009a pages 7-8 with fixes/optimizations from Canovas2010 pages 8-9.
    Implement a traversal iterator similar to the SDSL CSTs.
  */

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
