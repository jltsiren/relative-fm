#ifndef _RELATIVE_FM_RELATIVE_LCP_H
#define _RELATIVE_FM_RELATIVE_LCP_H

#include "rlz.h"
#include "support.h"

namespace relative
{

//------------------------------------------------------------------------------

class RelativeLCP
{
public:
  typedef SLArray       lcp_type;
  typedef int_vector<0> dlcp_type;
  typedef int_vector<0> index_type;

  const static uint64_t BRANCHING_FACTOR = 4; // For the minimal LCP value tree.
  const static std::string EXTENSION; // .rlcp

//------------------------------------------------------------------------------

  // Reference is an LCP array.
  RelativeLCP(const lcp_type& ref, const lcp_type& seq,
              const index_type& ref_sa, bool print = false);
  RelativeLCP(const lcp_type& ref, const std::string& base_name);
  RelativeLCP(const lcp_type& ref, std::istream& input);
  ~RelativeLCP();

  inline uint64_t size() const { return this->blocks.sum(); }

  uint64_t reportSize(bool print = false) const;
  void writeTo(const std::string& base_name) const;
  void writeTo(std::ostream& output) const;

//------------------------------------------------------------------------------

  // FIXME implement extracting multiple cells?
  inline uint64_t operator[] (uint64_t i) const
  {
    if(i >= this->size()) { return 0; }

    uint64_t phrase = this->blocks.inverse(i); // Phrase is 0-based.
    if(this->blocks.isLast(i)) { return this->samples[phrase + 1]; }

    uint64_t seq_pos = this->seqPos(phrase), ref_pos = this->refPos(phrase);
    uint64_t res = this->samples[phrase];
    if(ref_pos > 0) { res -= this->reference[ref_pos - 1]; }
    res += this->reference[ref_pos + i - seq_pos];
    return res;
  }

  /*
    Range minimum queries always return the leftmost position with minimum value.
  */
  inline uint64_t rmq(range_type range) const
  {
    return this->rmq(range.first, range.second);
  }

  uint64_t rmq(uint64_t from, uint64_t to) const;

  /*
    Return value >= size means "not found".
  */
  uint64_t psv(uint64_t pos) const;
  uint64_t nsv(uint64_t pos) const;

//------------------------------------------------------------------------------

  const lcp_type&   reference;
  int_vector<0>     phrases;
  CumulativeNZArray blocks;
  SLArray           samples;
  SLArray           tree;
  int_vector<64>    offsets;  // Offsets for each layer of the tree.

//------------------------------------------------------------------------------

private:
  void loadFrom(std::istream& input);
  void absoluteSamples(const lcp_type& lcp, const std::vector<uint64_t>& lengths);

  /*
    RMQ in [from, to] within the given phrase. Returns (lcp[i], i).
    from == 0 is interpreted as the beginning of the phrase.
    to >= sample_pos is interpreted as the sampled position (to + 1 should be a valid number).
  */
  range_type rmq(uint64_t phrase, uint64_t from, uint64_t to) const;

  /*
    PSV in the given phrase. Returns (psv_pos, LCP[pos]) or (psv_pos, val).
    If pos <= sample_pos, finds the psv_pos < pos with LCP[psv_pos] < LCP[pos].
    If pos > sample_pos, the entire phrase is considered, and the comparison is based on val.
    If the PSV does not exist, psv_pos will be >= size.
  */
  range_type psv(uint64_t phrase, uint64_t pos, uint64_t val) const;

//------------------------------------------------------------------------------

  inline uint64_t refPos(uint64_t phrase) const { return this->phrases[phrase]; }
  inline uint64_t seqPos(uint64_t phrase) const { return this->blocks.sum(phrase); }
  inline uint64_t samplePos(uint64_t phrase) const { return this->blocks.sum(phrase + 1) - 1; }

//------------------------------------------------------------------------------

  /*
    Tree operations. Node is identified by its offset in this->tree, while level refers
    to the level of the current node (leaves are at level 0).
  */

  inline uint64_t offset(uint64_t level) const { return this->offsets[level]; }
  inline uint64_t nodes(uint64_t level) const { return this->offsets[level + 1] - this->offsets[level]; }

  inline uint64_t findLevel(uint64_t node) const
  {
    uint64_t level = 0;
    while(this->offset(level + 1) <= node) { level++; }
    return level;
  }

  inline uint64_t parent(uint64_t node, uint64_t level) const
  {
    return this->offset(level + 1) + (node - this->offset(level)) / BRANCHING_FACTOR;
  }

  inline uint64_t child(uint64_t node, uint64_t level) const
  {
    return this->offset(level - 1) + (node - this->offset(level)) * BRANCHING_FACTOR;
  }

  inline uint64_t lastSibling(uint64_t first_child, uint64_t level) const
  {
    return std::min(this->offset(level + 1) - 1, first_child + BRANCHING_FACTOR - 1);
  }

  inline uint64_t firstSibling(uint64_t last_child, uint64_t level) const
  {
    return std::max(this->offset(level), last_child - (BRANCHING_FACTOR - 1));
  }

  inline uint64_t lastChild(uint64_t node, uint64_t level) const
  {
    return this->lastSibling(this->child(node, level), level - 1);
  }

  inline uint64_t root() const { return this->tree.size() - 1; }

//------------------------------------------------------------------------------

  RelativeLCP();
  RelativeLCP(const RelativeLCP&);
  RelativeLCP(RelativeLCP&&);
  RelativeLCP& operator=(const RelativeLCP&);
  RelativeLCP& operator==(RelativeLCP&);
};

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_RELATIVE_LCP_H
