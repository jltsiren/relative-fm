/*
  Copyright (c) 2015, 2016, 2017 Genome Research Ltd.

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

#include "new_relative_lcp.h"

namespace relative
{

//------------------------------------------------------------------------------

const std::string NewRelativeLCP::EXTENSION = ".rlcp";

//------------------------------------------------------------------------------

/*
  Tree operations on the range minima tree. If level is required, it refers to the
  level of the parameter node.
*/

inline size_type
rmtRoot(const NewRelativeLCP& lcp)
{
  return lcp.values() - 1;
}

inline size_type
rmtParent(const NewRelativeLCP& lcp, size_type node, size_type level)
{
  return lcp.offsets[level + 1] + (node - lcp.offsets[level]) / lcp.branching();
}

inline bool
rmtIsFirst(const NewRelativeLCP& lcp, size_type node, size_type level)
{
  return ((node - lcp.offsets[level]) % lcp.branching() == 0);
}

inline size_type
rmtFirstSibling(const NewRelativeLCP& lcp, size_type node, size_type level)
{
  return node - (node - lcp.offsets[level]) % lcp.branching();
}

inline size_type
rmtLastSibling(const NewRelativeLCP& lcp, size_type first_child, size_type level)
{
  return std::min(lcp.offsets[level + 1], first_child + lcp.branching()) - 1;
}

inline size_type
rmtFirstChild(const NewRelativeLCP& lcp, size_type node, size_type level)
{
  return lcp.offsets[level - 1] + (node - lcp.offsets[level]) * lcp.branching();
}

inline size_type
rmtLastChild(const NewRelativeLCP& lcp, size_type node, size_type level)
{
  return rmtLastSibling(lcp, rmtFirstChild(lcp, node, level), level - 1);
}

inline size_type
rmtLevel(const NewRelativeLCP& lcp, size_type node)
{
  size_type level = 0;
  while(lcp.offsets[level + 1] <= node) { level++; }
  return level;
}

//------------------------------------------------------------------------------

NewRelativeLCP::NewRelativeLCP(const lcp_type& ref, const lcp_type& seq) :
  reference(ref), branching_factor(BRANCHING_FACTOR)
{
  this->array = rlz::lcp::index_build<value_type>(seq.begin(), seq.end(), ref.begin(), ref.end());

 // Determine the number of levels.
  size_type level_count = 1, level_size = this->phrases();
  while(level_size > 1)
  {
    level_count++; level_size = (level_size + this->branching() - 1) / this->branching();
  }

  // Initialize offsets.
  this->offsets = sdsl::int_vector<64>(level_count + 1, 0);
  level_size = this->phrases();
  size_type total_size = 0;
  for(size_type level = 0; level < this->levels(); level++)
  {
    total_size += level_size;
    this->offsets[level + 1] = total_size;
    level_size = (level_size + this->branching() - 1) / this->branching();
  }

  // Initialize array.
  std::vector<value_type> tree_buffer(total_size, ~(value_type)0);
  for(size_type i = 0; i < this->phrases(); i++)
  {
    rlcp_type::iter begin, end;
    std::tie(begin, end) = this->array.iterator_phrase_id(i);
    while(begin != end) { tree_buffer[i] = std::min(tree_buffer[i], *begin); ++begin; }
  }
  for(size_type level = 0; level < this->levels(); level++)
  {
    for(size_type i = this->offsets[level]; i < this->offsets[level + 1]; i++)
    {
      size_type par = rmtParent(*this, i, level);
      tree_buffer[par] = std::min(tree_buffer[par], tree_buffer[i]);
    }
  }
  this->tree = SLArray(tree_buffer);
}

NewRelativeLCP::NewRelativeLCP(const lcp_type& ref, const std::string& base_name) :
  reference(ref)
{
  std::string filename = base_name + EXTENSION;
  std::ifstream input(filename, std::ios_base::binary);
  if(!input)
  {
    std::cerr << "NewRelativeLCP::NewRelativeLCP(): Cannot open input file " << filename << std::endl;
    return;
  }
  this->loadFrom(input); input.close();
}

NewRelativeLCP::NewRelativeLCP(const lcp_type& ref, std::istream& input) :
  reference(ref)
{
  this->loadFrom(input);
}

NewRelativeLCP::~NewRelativeLCP()
{
}

//------------------------------------------------------------------------------

size_type
NewRelativeLCP::reportSize(bool print) const
{
  size_type lcp_bytes = sdsl::size_in_bytes(this->array);
  size_type tree_bytes = sizeof(this->branching_factor) + sdsl::size_in_bytes(this->tree) + sdsl::size_in_bytes(this->offsets);
  size_type bytes = lcp_bytes + tree_bytes;

  if(print)
  {
#ifdef VERBOSE_OUTPUT
    printSize("LCP array", lcp_bytes, this->size());
    printSize("Tree", tree_bytes, this->size());
#endif
    printSize("Relative LCP", bytes, this->size());
  }

  return bytes;
}

void
NewRelativeLCP::writeTo(const std::string& base_name) const
{
  std::string filename = base_name + EXTENSION;
  std::ofstream output(filename, std::ios_base::binary);
  if(!output)
  {
    std::cerr << "NewRelativeLCP::writeTo(): Cannot open output file " << filename << std::endl;
    return;
  }
  this->writeTo(output); output.close();
}

void
NewRelativeLCP::writeTo(std::ostream& output) const
{
  this->array.serialize(output);
  sdsl::write_member(this->branching_factor, output);
  this->tree.serialize(output);
  this->offsets.serialize(output);
}

void
NewRelativeLCP::loadFrom(std::istream& input)
{
  this->array.load(input);
  this->array.set_source(rlz::iterator_container<rlz::alphabet::lcp<value_type>, SLArray::iterator>(this->reference.begin(), this->reference.end()));

  sdsl::read_member(this->branching_factor, input);
  this->tree.load(input);
  this->offsets.load(input);
}

//------------------------------------------------------------------------------

} // namespace relative
