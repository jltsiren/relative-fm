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
  uint64_t tree_bytes = size_in_bytes(this->tree);
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
}

void
RelativeLCP::loadFrom(std::istream& input)
{
  this->phrases.load(input);
  this->blocks.load(input);
  this->samples.load(input);
  this->tree.load(input);
}

//------------------------------------------------------------------------------

void
RelativeLCP::absoluteSamples(const lcp_type& lcp, const std::vector<uint64_t>& lengths)
{
  uint64_t nodes = 1;
  for(uint64_t next_nodes = lengths.size(); next_nodes > 1;
    next_nodes = (next_nodes + BRANCHING_FACTOR - 1) / BRANCHING_FACTOR)
  {
    nodes += next_nodes;
  }
  int_vector<0> tree_buffer(nodes, 0, bitlength(lcp.size() - 1));

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

  uint64_t tree_offset = 0, tree_nodes = lengths.size();
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
