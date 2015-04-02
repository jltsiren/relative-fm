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

  inline uint64_t operator[] (uint64_t i) const
  {
    if(i >= this->size()) { return 0; }

    uint64_t phrase = this->blocks.inverse(i); // Phrase is 0-based.
    if(this->blocks.isLast(i)) { return this->samples[phrase + 1]; }

    // Starting position of the phrase in seq and ref.
    uint64_t seq_pos = this->blocks.sum(phrase);
    uint64_t ref_pos = this->phrases[phrase];

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

  // FIXME implement extracting multiple cells?

  // FIXME what these should return?
  // FIXME implement
  uint64_t previousSmaller(uint64_t i) const;
  uint64_t nextSmaller(uint64_t i) const;

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
  inline range_type rmq(uint64_t phrase, uint64_t from, uint64_t to) const
  {
    // Starting position of the phrase in seq and ref.
    uint64_t seq_pos = this->blocks.sum(phrase);
    uint64_t ref_pos = this->phrases[phrase];

    // Determine the range and handle the special case where the phrase body is not needed.
    uint64_t sample_pos = this->blocks.sum(phrase + 1) - 1;
    if(from == 0) { from = seq_pos; }
    if(from == sample_pos) { return range_type(this->samples[phrase + 1], from); }
    uint64_t limit = std::min(to + 1, sample_pos);

    // Determine the minimum within the phrase body.
    uint64_t prev = this->reference[ref_pos + from - seq_pos];
    uint64_t curr = this->samples[phrase] + prev;
    if(ref_pos > 0) { curr -= this->reference[ref_pos - 1]; }
    range_type res(curr, from);
    for(uint64_t i = from + 1; i < limit; i++)
    {
      uint64_t next = this->reference[ref_pos + i - seq_pos];
      curr = curr + next - prev; prev = next;
      if(curr < res.first) { res.first = curr; res.second = i; }
    }

    // Check the sample if it falls within the range.
    if(to >= sample_pos)
    {
      uint64_t temp = this->samples[phrase + 1];
      if(temp < res.first) { res.first = temp; res.second = sample_pos; }
    }

    return res;
  }

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
