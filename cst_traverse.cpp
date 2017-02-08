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

#include <sdsl/suffix_trees.hpp>

#include "relative_cst.h"

using namespace relative;

//------------------------------------------------------------------------------

//#define USE_HASH

template<class CST>
void buildCST(CST& cst, const std::string& base_name, const std::string& type);

template<class CST>
void traverseHash(CST& cst, const std::string& name, size_type indent = 18);

template<class CST>
void traverse(CST& cst, const std::string& name, size_type indent = 18);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 3)
  {
    std::cerr << "Usage: cst_traverse ref seq1 [seq2 ...]" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  std::cout << "DFS traversal in compressed suffix trees" << std::endl;
  std::cout << std::endl;

  std::string ref_name = argv[1];
  std::cout << "Reference: " << ref_name << std::endl;
  std::cout << std::endl;

  SimpleFM<> ref_fm(ref_name);
  size_type fm_bytes = ref_fm.reportSize();
  printSize("FM-index", fm_bytes, ref_fm.size());
  RelativeLCP::lcp_type ref_lcp;
  sdsl::load_from_file(ref_lcp, ref_name + LCP_EXTENSION);
  size_type lcp_bytes = sdsl::size_in_bytes(ref_lcp);
  printSize("LCP array", lcp_bytes, ref_lcp.size());
  printSize("Reference data", fm_bytes + lcp_bytes, ref_fm.size());
  std::cout << std::endl;

  std::cout << std::endl;

  for(int arg = 2; arg < argc; arg++)
  {
    std::string seq_name = argv[arg];
    std::cout << "Sequence: " << seq_name << std::endl;
    std::cout << std::endl;

    {
      RelativeFM<> rfm(ref_fm, seq_name);
      RelativeLCP rlcp(ref_lcp, seq_name);
      RelativeCST<> rcst(rfm, rlcp);
      printSize("Relative CST", rcst.reportSize(), rcst.size());
#ifdef USE_HASH
      traverseHash(rcst, "Relative CST");
#else
      traverse(rcst, "Relative CST");
#endif
      std::cout << std::endl;
    }

    {
      std::string name = "cst_sct3_dac";
      sdsl::cst_sct3<> cst;
      buildCST(cst, seq_name, name);
#ifdef USE_HASH
      traverseHash(cst, name);
#else
      traverse(cst, name);
#endif
      std::cout << std::endl;
    }

    {
      std::string name = "cst_sct3_plcp";
      sdsl::cst_sct3<sdsl::csa_wt<>, sdsl::lcp_support_sada<>> cst;
      buildCST(cst, seq_name, name);
#ifdef USE_HASH
      traverseHash(cst, name);
#else
      traverse(cst, name);
#endif
      std::cout << std::endl;
    }

    {
      std::string name = "cst_sada";
      sdsl::cst_sada<> cst;
      buildCST(cst, seq_name, name);
      traverse(cst, name);
      std::cout << std::endl;
    }

    {
      std::string name = "cst_fully";
      sdsl::cst_fully<> cst;
      buildCST(cst, seq_name, name);
      traverse(cst, name);
      std::cout << std::endl;
    }

    std::cout << std::endl;
  }

  std::cout << "Memory used: " << inMegabytes(memoryUsage()) << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

bool
file_exists(const std::string& name)
{
  std::ifstream in(name.c_str(), std::ios_base::binary);
  if(!in) { return false; }
  in.close();
  return true;
}

template<class CST>
void
buildCST(CST& cst, const std::string& base_name, const std::string& type)
{
  std::string cst_file = base_name + "." + type;
  if(file_exists(cst_file))
  {
    sdsl::load_from_file(cst, cst_file);
  }
  else
  {
    construct(cst, base_name, 1);
    store_to_file(cst, cst_file);
  }
  printSize(type, sdsl::size_in_bytes(cst), cst.size());
}

//------------------------------------------------------------------------------

inline size_type
hashNode(const rcst_node& node, size_type hash)
{
  hash = fnv1a_hash(node.sp, hash);
  return fnv1a_hash(node.ep, hash);
}

inline size_type
hashNode(const sdsl::bp_interval<size_type>& node, size_type hash)
{
  hash = fnv1a_hash(node.i, hash);
  return fnv1a_hash(node.j, hash);
}

template<class CST>
void
traverseHash(CST& cst, const std::string& name, size_type indent)
{
  double start = readTimer();
  size_type nodes = 0, hash = FNV_OFFSET_BASIS;
  for(auto iter = cst.begin(); iter != cst.end(); ++iter)
  {
    if(iter.visit() == 1) { nodes++; }
    hash = hashNode(*iter, hash);
  }
  double seconds = readTimer() - start;

  std::string padding;
  if(name.length() + 1 < indent) { padding = std::string(indent - 1 - name.length(), ' '); }
  std::cout << name << ":" << padding << nodes << " nodes in " << seconds
            << " seconds (hash " << hash << ")" << std::endl;
}

template<class CST>
void
traverse(CST& cst, const std::string& name, size_type indent)
{
  double start = readTimer();
  size_type nodes = 0;
  for(auto iter = cst.begin(); iter != cst.end(); ++iter)
  {
    if(iter.visit() == 1) { nodes++; }
  }
  double seconds = readTimer() - start;

  std::string padding;
  if(name.length() + 1 < indent) { padding = std::string(indent - 1 - name.length(), ' '); }
  std::cout << name << ":" << padding << nodes << " nodes in " << seconds << " seconds" << std::endl;
}

//------------------------------------------------------------------------------
