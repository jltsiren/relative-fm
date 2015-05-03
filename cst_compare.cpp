/*
  Copyright (c) 2015 Genome Research Ltd.

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
#include <sdsl/suffix_tree_algorithm.hpp>

#include "relative_cst.h"

using namespace relative;

//------------------------------------------------------------------------------

//#define USE_HASH

template<class CST>
void buildCST(CST& cst, const std::string& base_name, const std::string& type);

void
buildSelect(RelativeFM<>& fm, const std::string& base_name);

template<class CST>
void matchingStatistics(const CST& cst, const int_vector<8>& seq, const std::string& name, uint64_t indent = 18);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 4)
  {
    std::cerr << "Usage: cst_compare reference target sequence" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  std::cout << "Finding the matching statistics with a CST" << std::endl;
  std::cout << std::endl;

  std::string ref_name = argv[1];
  std::cout << "Reference: " << ref_name << std::endl;

  SimpleFM<> ref_fm(ref_name);
  uint64_t fm_bytes = ref_fm.reportSize();
  printSize("FM-index", fm_bytes, ref_fm.size());
  RelativeLCP::lcp_type ref_lcp;
  load_from_file(ref_lcp, ref_name + LCP_EXTENSION);
  uint64_t lcp_bytes = size_in_bytes(ref_lcp);
  printSize("LCP array", lcp_bytes, ref_lcp.size());
  printSize("Reference data", fm_bytes + lcp_bytes, ref_fm.size());
  std::cout << std::endl;

  std::string target_name = argv[2];
  std::cout << "Target: " << target_name << std::endl;
  std::string seq_name = argv[3];
  std::cout << "Sequence: " << seq_name << std::endl;
  int_vector<8> seq;
  {
    std::ifstream in(seq_name.c_str(), std::ios_base::binary);
    if(!in)
    {
      std::cerr << "cst_compare: Cannot open sequence file " << seq_name << std::endl;
      return 2;
    }
    uint64_t size = util::file_size(seq_name); seq.resize(size);
    in.read((char*)(seq.data()), size); in.close();
  }
  std::cout << std::endl;

  {
    std::string name = "Relative CST";
    RelativeFM<> rfm(ref_fm, target_name);
    buildSelect(rfm, target_name);
    RelativeLCP rlcp(ref_lcp, target_name);
    RelativeCST<> rcst(rfm, rlcp);
    printSize("Relative CST", rcst.reportSize(), rcst.size());
    matchingStatistics(rcst, seq, name);
    std::cout << std::endl;
  }

/*    {
      std::string name = "cst_sct3_dac";
      cst_sct3<> cst;
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
      cst_sct3<csa_wt<>, lcp_support_sada<>> cst;
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
    cst_sada<> cst;
    buildCST(cst, target_name, name);
    maximalMatches(cst, seq, name);
    std::cout << std::endl;
  }*/

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
    load_from_file(cst, cst_file);
  }
  else
  {
    construct(cst, base_name, 1);
    store_to_file(cst, cst_file);
  }
  printSize(type, size_in_bytes(cst), cst.size());
}

void
buildSelect(RelativeFM<>& rfm, const std::string& base_name)
{
  if(rfm.fastSelect()) { return; }
  if(!(rfm.loadSelect(base_name)))
  {
    double start = readTimer();
    rfm.buildSelect();
    double seconds = readTimer() - start;
    std::cout << "Select structures built in " << seconds << " seconds" << std::endl;
    rfm.writeSelect(base_name);
  }
  printSize("Relative select", rfm.selectSize(), rfm.size());
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

template<class CST>
void
maximalMatch(const CST& cst, const int_vector<8>& seq,
  typename CST::node_type& prev, typename CST::node_type& next,
  uint64_t start_offset, typename CST::size_type& depth)
{
  typename CST::node_type old_next = next;
  typename CST::size_type next_depth = cst.depth(next);
  typename CST::size_type bwt_pos = cst.index.Psi(next.sp, depth);

  while(start_offset + depth < seq.size() &&
    cst.Psi(next, next_depth, depth, seq[start_offset + depth], bwt_pos))
  {
    depth++;
    if(next != old_next) { prev = old_next; }
  }
}

template<class CST>
void
matchingStatistics(const CST& cst, const int_vector<8>& seq, const std::string& name, uint64_t indent)
{
  double start = readTimer();
  typename CST::node_type prev = cst.root(), next = cst.root();
  typename CST::size_type depth = 0;

  maximalMatch(cst, seq, prev, next, 0, depth);
  uint64_t total_length = depth;
  for(uint64_t i = 1; i < seq.size(); i++)
  {
    if(depth == 0)
    {
      maximalMatch(cst, seq, prev, next, i, depth);
    }
    else
    {
      next = cst.sl(prev); depth--;
      typename CST::size_type next_depth = cst.depth(next);
      while(next_depth < depth)
      {
        typename CST::size_type bwt_pos;
        next = cst.child(next, seq[i + next_depth], bwt_pos);
        next_depth = cst.depth(next);
        if(next_depth <= depth) { prev = next; }
      }
      maximalMatch(cst, seq, prev, next, i, depth);
    }
    total_length += depth;
    std::cout << range_type(i, depth) << std::endl;
  }
  double seconds = readTimer() - start;

  std::string padding;
  if(name.length() + 1 < indent) { padding = std::string(indent - 1 - name.length(), ' '); }
  std::cout << name << ":" << padding << "Average maximal match: " << (total_length / (double)(cst.size()))
                    << " (" << seconds << " seconds)" << std::endl;
}

//------------------------------------------------------------------------------
