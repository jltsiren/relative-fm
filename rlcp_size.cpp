/*
  Copyright (c) 2017 Genome Research Ltd.

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

#include <cstdlib>
#include <unistd.h>

#include "new_relative_lcp.h"

using namespace relative;

//------------------------------------------------------------------------------

void
printRecursive(sdsl::structure_tree_node* v, size_type bytes)
{
  if(!(v->name.empty())) { printSize(v->name, v->size, bytes); }
  if(v->children.size() > 0 && v->name != "literals" && v->name != "parse")
  {
    for(const auto& child : v->children) { printRecursive(child.second.get(), bytes); }
  }
}

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 3)
  {
    std::cerr << "Usage: rlcp_size ref seq1 [seq2 ...]" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  std::cout << "RLCP size breakdown" << std::endl;
  std::cout << std::endl;

  std::string ref_name = argv[1];
  std::cout << "Reference: " << ref_name << std::endl;
  std::cout << std::endl;

  NewRelativeLCP::lcp_type reference;
  sdsl::load_from_file(reference, ref_name + LCP_EXTENSION);
  printSize("LCP", sdsl::size_in_bytes(reference), reference.size()); std::cout << std::endl;
  std::cout << std::endl;

  for(int arg = 2; arg < argc; arg++)
  {
    std::string seq_name = argv[arg];
    std::cout << "Target: " << seq_name << std::endl;
    std::cout << std::endl;

    NewRelativeLCP target(reference, seq_name);
    target.reportSize(true);
    std::cout << std::endl;

    std::unique_ptr<sdsl::structure_tree_node> st_node(new sdsl::structure_tree_node("", "type"));
    sdsl::nullstream ns;
    target.array.serialize(ns, st_node.get(), "RLZAP");
    printRecursive(st_node.get(), target.size());
    std::cout << std::endl;

    std::cout << std::endl;
  }

  std::cout << "Memory used: " << inMegabytes(memoryUsage()) << " MB" << std::endl;
  std::cout << std::endl;
  return 0;
}

//------------------------------------------------------------------------------
