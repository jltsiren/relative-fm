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

#include <functional>
#include <stack>

#include "relative_lcp.h"

namespace relative
{

//------------------------------------------------------------------------------

const std::string RelativeLCP::EXTENSION = ".rlcp";

RelativeLCP::RelativeLCP(const lcp_type& ref, const lcp_type& seq,
  const index_type& ref_sa, bool print) :
  reference(ref)
{
  std::vector<size_type> starts, lengths;
  relativeLZ<lcp_type, true>(seq, ref, ref_sa, starts, lengths, 0);
  if(print)
  {
    std::cout << "The RLZ parsing of the LCP array consists of " << starts.size() << " phrases" << std::endl;
  }

  this->absoluteSamples(seq, lengths);
  sdsl::util::assign(this->phrases, sdsl::int_vector<0>(starts.size(), 0, bit_length(ref.size() - 1)));
  for(size_type i = 0; i < starts.size(); i++) { this->phrases[i] = starts[i]; }
  sdsl::util::clear(starts);
  sdsl::util::assign(this->blocks, CumulativeNZArray(lengths)); sdsl::util::clear(lengths);
}

RelativeLCP::RelativeLCP(const lcp_type& ref, const std::string& base_name) :
  reference(ref)
{
  std::string filename = base_name + EXTENSION;
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "RelativeLCP::RelativeLCP(): Cannot open input file " << filename << std::endl;
    return;
  }
  this->loadFrom(input); input.close();
}

RelativeLCP::RelativeLCP(const lcp_type& ref, std::istream& input) :
  reference(ref)
{
  this->loadFrom(input);
}

RelativeLCP::~RelativeLCP()
{
}

//------------------------------------------------------------------------------

size_type
RelativeLCP::reportSize(bool print) const
{
  size_type phrase_bytes = sdsl::size_in_bytes(this->phrases);
  size_type block_bytes = sdsl::size_in_bytes(this->blocks);
  size_type sample_bytes = sdsl::size_in_bytes(this->samples);
  size_type tree_bytes = sdsl::size_in_bytes(this->tree) + sdsl::size_in_bytes(this->offsets);
  size_type bytes = phrase_bytes + block_bytes + sample_bytes + tree_bytes;

  if(print)
  {
#ifdef VERBOSE_OUTPUT
    printSize("Phrases", phrase_bytes, this->size());
    printSize("Blocks", block_bytes, this->size());
    printSize("Samples", sample_bytes, this->size());
    printSize("Tree", tree_bytes, this->size());
#endif
    printSize("Relative LCP", bytes, this->size());
    std::cout << std::endl;
  }

  return bytes;
}

void
RelativeLCP::writeTo(const std::string& base_name) const
{
  std::string filename = base_name + EXTENSION;
  std::ofstream output(filename.c_str(), std::ios_base::binary);
  if(!output)
  {
    std::cerr << "RelativeLCP::writeTo(): Cannot open output file " << filename << std::endl;
    return;
  }
  this->writeTo(output); output.close();
}

void
RelativeLCP::writeTo(std::ostream& output) const
{
  this->phrases.serialize(output);
  this->blocks.serialize(output);
  this->samples.serialize(output);
  this->tree.serialize(output);
  this->offsets.serialize(output);
}

void
RelativeLCP::loadFrom(std::istream& input)
{
  this->phrases.load(input);
  this->blocks.load(input);
  this->samples.load(input);
  this->tree.load(input);
  this->offsets.load(input);
}

//------------------------------------------------------------------------------

inline void
updateRes(const RelativeLCP& lcp, range_type& res, size_type i)
{
  size_type temp = lcp.tree[i];
  if(temp < res.second) { res.first = i; res.second = temp; }
}

range_type
RelativeLCP::rmq(size_type from, size_type to) const
{
  to = std::min(to, this->size() - 1);
  if(from > to) { return range_type(this->size(), 0); }

  size_type phrase_from = this->blocks.inverse(from), phrase_to = this->blocks.inverse(to);
  range_type res(this->size(), this->size());

  // Process the full blocks first.
  // FIXME The first and the last block could also be full...
  if(phrase_to == phrase_from + 2)  // Just a single full phrase.
  {
    res = this->rmq(phrase_from + 1, 0, this->size());
  }
  else if(phrase_to > phrase_from + 2)
  {
    /*
      Invariants:
        - left < right
        - nodes before subtree(left) are processed
        - nodes after subtree(right) are in tail
    */
    size_type level = 0, left = phrase_from + 1, right = phrase_to - 1;
    std::stack<range_type> tail;  // Process the right tail last.
    while(true)
    {
      size_type left_par = this->parent(left, level), right_par = this->parent(right, level);
      if(left_par == right_par)
      {
        for(size_type i = left; i <= right; i++) { updateRes(*this, res, i); }
        break;
      }
      else
      {
        size_type left_child = this->child(left_par, level + 1);
        if(left != left_child)
        {
          size_type last_child = this->lastSibling(left_child, level);
          for(size_type i = left; i <= last_child; i++) { updateRes(*this, res, i); }
          left_par++;
        }

        size_type right_child = this->lastChild(right_par, level + 1);
        if(right != right_child)
        {
          size_type first_child = this->firstSibling(right_child, level);
          for(size_type i = right; i >= first_child; i--) { tail.push(range_type(i, this->tree[i])); }
          right_par--;
        }
        if(left_par >= right_par)
        {
          if(left_par == right_par) { updateRes(*this, res, left_par); }
          break;
        }
      }
      left = left_par; right = right_par; level++;
    }

    // Process the tail and proceed to the leaf to find the correct block.
    while(!(tail.empty()))
    {
      range_type temp = tail.top(); tail.pop();
      if(temp.second < res.second) { res = temp; }
    }
    level = this->findLevel(res.first);
    while(level > 0)
    {
      res.first = this->child(res.first, level);
      while(this->tree[res.first] != res.second) { res.first++; }
      level--;
    }
    res = this->rmq(res.first, 0, this->size());
  }

  // Process the partial blocks.
  if(this->tree[phrase_from] <= res.second)
  {
    range_type temp = this->rmq(phrase_from, from, to);
    if(temp.second <= res.second) { res = temp; }
  }
  if(phrase_from < phrase_to && this->tree[phrase_to] < res.second)
  {
    range_type temp = this->rmq(phrase_to, 0, to);
    if(temp.second < res.second) { res = temp; }
  }

  return res;
}

range_type
RelativeLCP::rmq(size_type phrase, size_type from, size_type to) const
{
  size_type seq_pos = this->seqPos(phrase);
  size_type ref_pos = this->refPos(phrase);
  size_type sample_pos = this->samplePos(phrase);

  // Determine the range and handle the special case where the phrase body is not needed.
  if(from == 0) { from = seq_pos; }
  if(from == sample_pos) { return range_type(from, this->samples[phrase + 1]); }
  size_type limit = std::min(to + 1, sample_pos);

  // Determine the minimum within the phrase body.
  // If from is far from seq_pos, it is probably cheaper to access the reference directly.
  size_type r_pos = ref_pos + from - seq_pos;
  size_type rank = this->reference.initForward(r_pos);
  size_type prev = this->reference.accessForward(r_pos, rank);
  size_type curr = this->samples[phrase] + prev;
  if(ref_pos > 0) { curr -= this->reference[ref_pos - 1]; }
  range_type res(from, curr);
  for(size_type i = from + 1; i < limit; i++)
  {
    r_pos++;
    size_type next = this->reference.accessForward(r_pos, rank);
    curr = curr + next - prev; prev = next;
    if(curr < res.second) { res.first = i; res.second = curr; }
  }

  // Check the sample if it falls within the range.
  if(to >= sample_pos)
  {
    size_type temp = this->samples[phrase + 1];
    if(temp < res.second) { res.first = sample_pos; res.second = temp; }
  }

  return res;
}

//------------------------------------------------------------------------------

/*
  PSV in the given phrase. Returns (psv_pos, LCP[psv_pos]).
  Updates val, if it has not been specified.
  If pos <= sample_pos, finds the psv_pos < pos with comp(LCP[psv_pos], LCP[pos]).
  If pos > sample_pos, the entire phrase is considered, and the comparison is based on val.
  If the PSV does not exist, psv_pos will be >= size.
*/
template<class Comp>
range_type
psv(const RelativeLCP& lcp, size_type phrase, size_type pos, size_type& val, const Comp& comp)
{
  size_type seq_pos = lcp.seqPos(phrase);
  size_type ref_pos = lcp.refPos(phrase);
  size_type sample_pos = lcp.samplePos(phrase);

  // Set val to the actual LCP we are comparing against and pos to the upper bound
  // of considered positions. Also handle the special cases.
  if(pos > sample_pos)
  {
    size_type temp = lcp.samples[phrase + 1];
    if(comp(temp, val)) { return range_type(sample_pos, temp); }
    pos = sample_pos;
  }
  else
  {
    if(pos == sample_pos) { val = lcp.samples[phrase + 1]; }
    else
    {
      val = lcp.samples[phrase] + lcp.reference[ref_pos + pos - seq_pos];
      if(ref_pos > 0) { val -= lcp.reference[ref_pos - 1]; }
    }
    if(!comp(lcp.tree[phrase], val)) { return range_type(lcp.size(), 0); }
  }
  if(pos <= seq_pos) { return range_type(lcp.size(), 0); }

  // Handle the phrase body.
  size_type r_pos = ref_pos + pos - 1 - seq_pos;
  size_type rank = lcp.reference.initBackward(r_pos);
  size_type prev = lcp.reference.accessBackward(r_pos, rank);
  size_type curr = lcp.samples[phrase] + prev;
  if(ref_pos > 0) { curr -= lcp.reference[ref_pos - 1]; }
  while(pos > seq_pos)
  {
    pos--;
    if(comp(curr, val)) { return range_type(pos, curr); }
    size_type temp = (r_pos > 0 ? lcp.reference.accessBackward(r_pos - 1, rank) : 0); r_pos--;
    curr = curr + temp - prev; prev = temp;
  }

  return range_type(lcp.size(), 0);
}

template<class Comp>
range_type
psv(const RelativeLCP& lcp, size_type pos, Comp comp)
{
  if(pos == 0 || pos >= lcp.size()) { return range_type(lcp.size(), 0); }

  size_type phrase = lcp.blocks.inverse(pos);
  size_type val = 0; // lcp[pos]
  range_type res = psv(lcp, phrase, pos, val, comp);
  if(res.first < lcp.size() || phrase == 0) { return res; } // psv in the same phrase or does not exist.
  if(comp(lcp.tree[phrase - 1], val)) // psv is in the previous phrase.
  {
    return psv(lcp, phrase - 1, lcp.size(), val, comp);
  }

  // Go upward until psv is in the current subtree.
  size_type tree_node = phrase - 1, level = 0;
  bool found = false;
  while(!found)
  {
    if(tree_node == lcp.root()) { return range_type(lcp.size(), 0); }
    size_type parent_node = lcp.parent(tree_node, level); level++;
    if(comp(lcp.tree[parent_node], val))
    {
      // The smaller value may come from a subtree after leaf 'phrase'.
      size_type first_child = lcp.child(parent_node, level);
      for(size_type i = tree_node; i > first_child && !found; i--)
      {
        if(comp(lcp.tree[i - 1], val)) { tree_node = i - 1; level--; found = true; }
      }
    }
    if(!found) { tree_node = parent_node; }
  }

  // Find the leaf where the psv is and then locate it.
  while(level > 0)
  {
    tree_node = lcp.lastChild(tree_node, level); level--;
    while(!comp(lcp.tree[tree_node], val)) { tree_node--; }
  }
  return psv(lcp, tree_node, lcp.size(), val, comp);
}

range_type
RelativeLCP::psv(size_type pos) const
{
  return relative::psv(*this, pos, std::less<size_type>());
}

range_type
RelativeLCP::psev(size_type pos) const
{
  return relative::psv(*this, pos, std::less_equal<size_type>());
}

//------------------------------------------------------------------------------

/*
  NSV in the given phrase. Returns (nsv_pos, LCP[nsv_pos]).
  Updates val, if it has not been specified.
  If pos >= seq_pos, finds the nsv_pos > pos with comp(LCP[nsv_pos], LCP[pos]).
  If pos < seq_pos, the entire phrase is considered, and the comparison is based on val.
  If the NSV does not exist, nsv_pos will be >= size.
*/
template<class Comp>
range_type
nsv(const RelativeLCP& lcp, size_type phrase, size_type pos, size_type& val, const Comp& comp)
{
  size_type seq_pos = lcp.seqPos(phrase);
  size_type ref_pos = lcp.refPos(phrase);
  size_type sample_pos = lcp.samplePos(phrase);

  // Set val to the actual LCP we are comparing against and pos to the lower bound
  // of considered positions. Also handle the special cases.
  if(pos < seq_pos) { pos = seq_pos; }
  else
  {
    if(pos == sample_pos) { val = lcp.samples[phrase + 1]; return range_type(lcp.size(), 0); }
    val = lcp.samples[phrase] + lcp.reference[ref_pos + pos - seq_pos];
    if(ref_pos > 0) { val -= lcp.reference[ref_pos - 1]; }
    pos++;
    if(!comp(lcp.tree[phrase], val)) { return range_type(lcp.size(), 0); }
  }

  // Handle the phrase body.
  if(pos < sample_pos)
  {
    size_type r_pos = ref_pos + pos - seq_pos;
    size_type rank = lcp.reference.initForward(r_pos);
    size_type prev = lcp.reference.accessForward(r_pos, rank);
    size_type curr = lcp.samples[phrase] + prev;
    if(ref_pos > 0) { curr -= lcp.reference[ref_pos - 1]; }
    while(pos < sample_pos)
    {
      if(comp(curr, val)) { return range_type(pos, curr); }
      pos++; r_pos++;
      size_type temp = lcp.reference.accessForward(r_pos, rank);
      curr = curr + temp - prev; prev = temp;
    }
  }

  // Look at the sample.
  size_type temp = lcp.samples[phrase + 1];
  if(comp(temp, val)) { return range_type(sample_pos, temp); }

  return range_type(lcp.size(), 0);
}

template<class Comp>
range_type
nsv(const RelativeLCP& lcp, size_type pos, Comp comp)
{
  if(pos + 1 >= lcp.size()) { return range_type(lcp.size(), 0); }

  size_type phrase = lcp.blocks.inverse(pos);
  size_type val = 0; // lcp[pos]
  range_type res = nsv(lcp, phrase, pos, val, comp);
  // nsv in the same phrase or does not exist.
  if(res.first < lcp.size() || phrase + 1 >= lcp.phrases.size()) { return res; }
  if(comp(lcp.tree[phrase + 1], val)) // nsv is in the next phrase.
  {
    return nsv(lcp, phrase + 1, 0, val, comp);
  }

  // Go upward until nsv is in the current subtree.
  size_type tree_node = phrase + 1, level = 0;
  bool found = false;
  while(!found)
  {
    if(tree_node == lcp.root()) { return range_type(lcp.size(), 0); }
    size_type parent_node = lcp.parent(tree_node, level); level++;
    if(comp(lcp.tree[parent_node], val))
    {
      // The smaller value may come from a subtree before leaf 'phrase'.
      size_type last_child = lcp.lastChild(parent_node, level);
      for(size_type i = tree_node + 1; i <= last_child && !found; i++)
      {
        if(comp(lcp.tree[i], val)) { tree_node = i; level--; found = true; }
      }
    }
    if(!found) { tree_node = parent_node; }
  }

  // Find the leaf where the nsv is and then locate it.
  while(level > 0)
  {
    tree_node = lcp.child(tree_node, level); level--;
    while(!comp(lcp.tree[tree_node], val)) { tree_node++; }
  }
  return nsv(lcp, tree_node, 0, val, comp);
}

range_type
RelativeLCP::nsv(size_type pos) const
{
  return relative::nsv(*this, pos, std::less<size_type>());
}

range_type
RelativeLCP::nsev(size_type pos) const
{
  return relative::nsv(*this, pos, std::less_equal<size_type>());
}

//------------------------------------------------------------------------------

void
RelativeLCP::absoluteSamples(const lcp_type& lcp, const std::vector<size_type>& lengths)
{
  size_type tree_nodes = 0;
  std::vector<size_type> offset_buffer;
  for(size_type next_nodes = lengths.size(); next_nodes > 1;
    next_nodes = (next_nodes + BRANCHING_FACTOR - 1) / BRANCHING_FACTOR)
  {
    offset_buffer.push_back(tree_nodes);
    tree_nodes += next_nodes;
  }
  offset_buffer.push_back(tree_nodes); tree_nodes++;  // Add the root.
  offset_buffer.push_back(tree_nodes);                // Add a guard.
  sdsl::int_vector<0> tree_buffer(tree_nodes, 0, bit_length(lcp.size() - 1));
  sdsl::util::assign(this->offsets, sdsl::int_vector<64>(offset_buffer.size()));
  for(size_type i = 0; i < offset_buffer.size(); i++) { this->offsets[i] = offset_buffer[i]; }

  // We add a sample before the first phrase.
  sdsl::int_vector<0> sample_buffer(lengths.size() + 1, 0, bit_length(lcp.size() - 1));
  for(size_type phrase = 0, pos = 0; phrase < lengths.size(); phrase++)
  {
    size_type min_val = ~(size_type)0, limit = pos + lengths[phrase];
    sdsl::int_vector<64> buffer = lcp.extract(pos, limit);
    for(size_type i = 0; i < buffer.size(); i++) { min_val = std::min(min_val, buffer[i]); }
    tree_buffer[phrase] = min_val;
    sample_buffer[phrase + 1] = buffer[buffer.size() - 1];
    pos = limit;
  }

  size_type tree_offset = 0; tree_nodes = lengths.size();
  while(tree_nodes > 1)
  {
    size_type limit = tree_offset + tree_nodes; tree_nodes = 0;
    while(tree_offset < limit)
    {
      size_type block = std::min(tree_offset + BRANCHING_FACTOR, limit);
      size_type min_val = tree_buffer[tree_offset]; tree_offset++;
      while(tree_offset < block)
      {
        if(tree_buffer[tree_offset] < min_val) { min_val = tree_buffer[tree_offset]; }
        tree_offset++;
      }
      tree_buffer[limit + tree_nodes] = min_val; tree_nodes++;
    }
  }

  sdsl::util::assign(this->samples, SLArray(sample_buffer)); sdsl::util::clear(sample_buffer);
  sdsl::util::assign(this->tree, SLArray(tree_buffer)); sdsl::util::clear(tree_buffer);
}

//------------------------------------------------------------------------------

} // namespace relative
