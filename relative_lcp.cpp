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
  std::vector<uint64_t> starts, lengths;
  relativeLZ<lcp_type, true>(seq, ref, ref_sa, starts, lengths, 0);
  if(print)
  {
    std::cout << "The RLZ parsing of the LCP array consists of " << starts.size() << " phrases" << std::endl;
  }

  this->absoluteSamples(seq, lengths);
  util::assign(this->phrases, int_vector<0>(starts.size(), 0, bitlength(ref.size() - 1)));
  for(uint64_t i = 0; i < starts.size(); i++) { this->phrases[i] = starts[i]; }
  util::clear(starts);
  util::assign(this->blocks, CumulativeNZArray(lengths)); util::clear(lengths);
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

uint64_t
RelativeLCP::reportSize(bool print) const
{
  uint64_t phrase_bytes = size_in_bytes(this->phrases);
  uint64_t block_bytes = size_in_bytes(this->blocks);
  uint64_t sample_bytes = size_in_bytes(this->samples);
  uint64_t tree_bytes = size_in_bytes(this->tree) + size_in_bytes(this->offsets);
  uint64_t bytes = phrase_bytes + block_bytes + sample_bytes + tree_bytes;

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
updateRes(const RelativeLCP& lcp, range_type& res, uint64_t i)
{
  uint64_t temp = lcp.tree[i];
  if(temp < res.first) { res.first = temp; res.second = i; }
}

uint64_t
RelativeLCP::rmq(uint64_t from, uint64_t to) const
{
  to = std::min(to, this->size() - 1);
  if(from > to) { return this->size(); }

  uint64_t phrase_from = this->blocks.inverse(from), phrase_to = this->blocks.inverse(to);
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
    uint64_t level = 0, left = phrase_from + 1, right = phrase_to - 1;
    std::stack<range_type> tail;  // Process the right tail last.
    while(true)
    {
      uint64_t left_par = this->parent(left, level), right_par = this->parent(right, level);
      if(left_par == right_par)
      {
        for(uint64_t i = left; i <= right; i++) { updateRes(*this, res, i); }
        break;
      }
      else
      {
        uint64_t left_child = this->child(left_par, level + 1);
        if(left != left_child)
        {
          uint64_t last_child = this->lastSibling(left_child, level);
          for(uint64_t i = left; i <= last_child; i++) { updateRes(*this, res, i); }
          left_par++;
        }

        uint64_t right_child = this->lastChild(right_par, level + 1);
        if(right != right_child)
        {
          uint64_t first_child = this->firstSibling(right_child, level);
          for(uint64_t i = right; i >= first_child; i--) { tail.push(range_type(this->tree[i], i)); }
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
      if(temp.first < res.first) { res = temp; }
    }
    level = this->findLevel(res.second);
    while(level > 0)
    {
      res.second = this->child(res.second, level);
      while(this->tree[res.second] != res.first) { res.second++; }
      level--;
    }
    res = this->rmq(res.second, 0, this->size());
  }

  // Process the partial blocks.
  if(this->tree[phrase_from] <= res.first)
  {
    range_type temp = this->rmq(phrase_from, from, to);
    if(temp.first <= res.first) { res = temp; }
  }
  if(phrase_from < phrase_to && this->tree[phrase_to] < res.first)
  {
    range_type temp = this->rmq(phrase_to, 0, to);
    if(temp.first < res.first) { res = temp; }
  }

  return res.second;
}

//------------------------------------------------------------------------------

void
RelativeLCP::absoluteSamples(const lcp_type& lcp, const std::vector<uint64_t>& lengths)
{
  uint64_t tree_nodes = 0;
  std::vector<uint64_t> offset_buffer;
  for(uint64_t next_nodes = lengths.size(); next_nodes > 1;
    next_nodes = (next_nodes + BRANCHING_FACTOR - 1) / BRANCHING_FACTOR)
  {
    offset_buffer.push_back(tree_nodes);
    tree_nodes += next_nodes;
  }
  offset_buffer.push_back(tree_nodes); tree_nodes++;  // Add the root.
  offset_buffer.push_back(tree_nodes);                // Add a guard.
  int_vector<0> tree_buffer(tree_nodes, 0, bitlength(lcp.size() - 1));
  util::assign(this->offsets, int_vector<64>(offset_buffer.size()));
  for(uint64_t i = 0; i < offset_buffer.size(); i++) { this->offsets[i] = offset_buffer[i]; }

  // We add a sample before the first phrase.
  int_vector<0> sample_buffer(lengths.size() + 1, 0, bitlength(lcp.size() - 1));
  for(uint64_t phrase = 0, pos = 0; phrase < lengths.size(); phrase++)
  {
    uint64_t min_val = ~(uint64_t)0, limit = pos + lengths[phrase];
    int_vector<64> buffer = lcp.extract(pos, limit);
    for(uint64_t i = 0; i < buffer.size(); i++) { min_val = std::min(min_val, buffer[i]); }
    tree_buffer[phrase] = min_val;
    sample_buffer[phrase + 1] = buffer[buffer.size() - 1];
    pos = limit;
  }

  uint64_t tree_offset = 0; tree_nodes = lengths.size();
  while(tree_nodes > 1)
  {
    uint64_t limit = tree_offset + tree_nodes; tree_nodes = 0;
    while(tree_offset < limit)
    {
      uint64_t block = std::min(tree_offset + BRANCHING_FACTOR, limit);
      uint64_t min_val = tree_buffer[tree_offset]; tree_offset++;
      while(tree_offset < block)
      {
        if(tree_buffer[tree_offset] < min_val) { min_val = tree_buffer[tree_offset]; }
        tree_offset++;
      }
      tree_buffer[limit + tree_nodes] = min_val; tree_nodes++;
    }
  }

  util::assign(this->samples, SLArray(sample_buffer)); util::clear(sample_buffer);
  util::assign(this->tree, SLArray(tree_buffer)); util::clear(tree_buffer);
}

//------------------------------------------------------------------------------

} // namespace relative
