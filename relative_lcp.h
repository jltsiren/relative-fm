/*
  Copyright (c) 2015, 2016 Genome Research Ltd.

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
  typedef sdsl::int_vector<0> dlcp_type;
  typedef sdsl::int_vector<0> index_type;

  const static size_type BRANCHING_FACTOR = 4; // For the minimal LCP value tree.
  const static std::string EXTENSION; // .rlcp

//------------------------------------------------------------------------------

  // Reference is an LCP array.
  RelativeLCP(const lcp_type& ref, const lcp_type& seq,
              const index_type& ref_sa, bool print = false);
  RelativeLCP(const lcp_type& ref, const std::string& base_name);
  RelativeLCP(const lcp_type& ref, std::istream& input);
  ~RelativeLCP();

  inline size_type size() const { return this->blocks.sum(); }

  size_type reportSize(bool print = false) const;
  void writeTo(const std::string& base_name) const;
  void writeTo(std::ostream& output) const;

//------------------------------------------------------------------------------

  // FIXME implement extracting multiple cells?
  inline size_type operator[] (size_type i) const
  {
    if(i >= this->size()) { return 0; }

    size_type phrase = this->blocks.inverse(i); // Phrase is 0-based.
    if(this->blocks.isLast(i)) { return this->samples[phrase + 1]; }

    size_type seq_pos = this->seqPos(phrase), ref_pos = this->refPos(phrase);
    size_type res = this->samples[phrase];
    if(ref_pos > 0) { res -= this->reference[ref_pos - 1]; }
    res += this->reference[ref_pos + i - seq_pos];
    return res;
  }

  /*
    Range minimum queries always return the leftmost position with minimum value.
    The return value is (pos, LCP[pos]).
  */
  inline range_type rmq(range_type range) const
  {
    return this->rmq(range.first, range.second);
  }

  range_type rmq(size_type from, size_type to) const;

  /*
    The return value is (res, LCP[res]). If res >= size, no value was found.
  */
  range_type psv(size_type pos) const;
  range_type psev(size_type pos) const;
  range_type nsv(size_type pos) const;
  range_type nsev(size_type pos) const;

//------------------------------------------------------------------------------

  const lcp_type&   reference;
  sdsl::int_vector<0>     phrases;
  CumulativeNZArray blocks;
  SLArray           samples;
  SLArray           tree;
  sdsl::int_vector<64>    offsets;  // Offsets for each layer of the tree.

//------------------------------------------------------------------------------

  inline size_type refPos(size_type phrase) const { return this->phrases[phrase]; }
  inline size_type seqPos(size_type phrase) const { return this->blocks.sum(phrase); }
  inline size_type samplePos(size_type phrase) const { return this->blocks.sum(phrase + 1) - 1; }

//------------------------------------------------------------------------------

  /*
    Tree operations. Node is identified by its offset in this->tree, while level refers
    to the level of the current node (leaves are at level 0).
  */

  inline size_type offset(size_type level) const { return this->offsets[level]; }
  inline size_type nodes(size_type level) const { return this->offsets[level + 1] - this->offsets[level]; }

  inline size_type findLevel(size_type node) const
  {
    size_type level = 0;
    while(this->offset(level + 1) <= node) { level++; }
    return level;
  }

  inline size_type parent(size_type node, size_type level) const
  {
    return this->offset(level + 1) + (node - this->offset(level)) / BRANCHING_FACTOR;
  }

  inline size_type child(size_type node, size_type level) const
  {
    return this->offset(level - 1) + (node - this->offset(level)) * BRANCHING_FACTOR;
  }

  inline size_type lastSibling(size_type first_child, size_type level) const
  {
    return std::min(this->offset(level + 1) - 1, first_child + BRANCHING_FACTOR - 1);
  }

  inline size_type firstSibling(size_type last_child, size_type level) const
  {
    return last_child - (last_child - this->offset(level)) % BRANCHING_FACTOR;
  }

  inline size_type lastChild(size_type node, size_type level) const
  {
    return this->lastSibling(this->child(node, level), level - 1);
  }

  inline size_type root() const { return this->tree.size() - 1; }

//------------------------------------------------------------------------------

private:
  void loadFrom(std::istream& input);
  void absoluteSamples(const lcp_type& lcp, const std::vector<size_type>& lengths);

  /*
    RMQ in [from, to] within the given phrase. Returns (i, lcp[i]).
    from == 0 is interpreted as the beginning of the phrase.
    to >= sample_pos is interpreted as the sampled position (to + 1 should be a valid number).
  */
  range_type rmq(size_type phrase, size_type from, size_type to) const;

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
