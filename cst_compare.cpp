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

#include <sdsl/suffix_trees.hpp>
#include <sdsl/suffix_tree_algorithm.hpp>

#include "relative_cst.h"

using namespace relative;

//------------------------------------------------------------------------------

//#define VERIFY_RESULTS

template<class CST>
void buildCST(CST& cst, const std::string& base_name, const std::string& type);

void buildSelect(RelativeFM<>& fm, const std::string& base_name);

template<class CST>
void matchingStatistics(const CST& cst, const sdsl::int_vector<8>& seq,
  std::vector<range_type>& ranges, std::vector<size_type>& depths,
  const std::string& name);

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
  size_type fm_bytes = ref_fm.reportSize();
  printSize("FM-index", fm_bytes, ref_fm.size());
  RelativeLCP::lcp_type ref_lcp;
  sdsl::load_from_file(ref_lcp, ref_name + LCP_EXTENSION);
  size_type lcp_bytes = sdsl::size_in_bytes(ref_lcp);
  printSize("LCP array", lcp_bytes, ref_lcp.size());
  printSize("Reference data", fm_bytes + lcp_bytes, ref_fm.size());
  std::cout << std::endl;

  std::string target_name = argv[2];
  std::cout << "Target: " << target_name << std::endl;
  std::string seq_name = argv[3];
  sdsl::int_vector<8> seq;
  {
    std::ifstream in(seq_name.c_str(), std::ios_base::binary);
    if(!in)
    {
      std::cerr << "cst_compare: Cannot open sequence file " << seq_name << std::endl;
      return 2;
    }
    size_type size = sdsl::util::file_size(seq_name); seq.resize(size);
    in.read((char*)(seq.data()), size); in.close();
  }
  std::cout << "Sequence: " << seq_name << " (" << seq.size() << " bytes)" << std::endl;
  std::cout << std::endl;

  std::vector<range_type> rcst_ranges;
  std::vector<size_type>   rcst_depths;
  {
    std::string name = "Relative (slow)";
    RelativeFM<> rfm(ref_fm, target_name);
    RelativeLCP rlcp(ref_lcp, target_name);
    RelativeCST<> rcst(rfm, rlcp);
    printSize(name, rcst.reportSize(), rcst.size());
    matchingStatistics(rcst, seq, rcst_ranges, rcst_depths, name);
    std::cout << std::endl;

    name = "Relative (fast)";
    buildSelect(rfm, target_name);
    printSize(name, rcst.reportSize(), rcst.size());
    matchingStatistics(rcst, seq, rcst_ranges, rcst_depths, name);
    std::cout << std::endl;
  }

  std::vector<range_type> cst_ranges;
  std::vector<size_type>   cst_depths;
  {
    std::string name = "cst_sct3_dac";
    sdsl::cst_sct3<> cst;
    buildCST(cst, target_name, name);
    matchingStatistics(cst, seq, cst_ranges, cst_depths, name);
    std::cout << std::endl;
  }

  {
    std::string name = "cst_sct3_plcp";
    sdsl::cst_sct3<sdsl::csa_wt<>, sdsl::lcp_support_sada<>> cst;
    buildCST(cst, target_name, name);
    matchingStatistics(cst, seq, cst_ranges, cst_depths, name);
    std::cout << std::endl;
  }

  {
    std::string name = "cst_sada";
    sdsl::cst_sada<> cst;
    buildCST(cst, target_name, name);
    matchingStatistics(cst, seq, cst_ranges, cst_depths, name);
    std::cout << std::endl;
  }

  {
    std::string name = "cst_fully";
    sdsl::cst_fully<> cst;
    buildCST(cst, target_name, name);
    // cst_fully is either buggy or orders of magnitude slower than the others,
    // or its interface is incompatible.
    //matchingStatistics(cst, seq, cst_ranges, cst_depths, name);
    std::cout << std::endl;
  }

#ifdef VERIFY_RESULTS
  for(size_type i = 0; i < seq.size(); i++)
  {
    if(rcst_ranges[i] != cst_ranges[i] || rcst_depths[i] != cst_depths[i])
    {
      std::cerr << "cst_compare: Matching statistics for " << seq_name << "[" << i << "]:" << std::endl;
      std::cerr << "  Relative CST: range " << rcst_ranges[i] << ", depth " << rcst_depths[i] << std::endl;
      std::cerr << "  CST:          range " << cst_ranges[i] << ", depth " << cst_depths[i] << std::endl;
      break;
    }
  }
  std::cout << "Matching statistics verified." << std::endl;
  std::cout << std::endl;
#endif

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
    std::cout << std::endl;
    rfm.writeSelect(base_name);
  }
}

//------------------------------------------------------------------------------

template<class CST>
void
maximalMatch(const CST& cst, const sdsl::int_vector<8>& seq,
  typename CST::node_type& prev, typename CST::node_type& next,
  size_type start_offset, typename CST::size_type& depth,
  typename CST::size_type& next_depth)
{
  typename CST::size_type bwt_pos =
    (depth >= next_depth ? 0 : get_char_pos(cst.lb(next), depth - 1, cst.csa));

  while(start_offset + depth < seq.size())
  {
    auto comp = cst.csa.char2comp[seq[start_offset + depth]];
    if(comp == 0 && seq[start_offset + depth] != 0) { break; }
    if(depth >= next_depth) // Next node reached, follow a new edge.
    {
      typename CST::node_type temp = cst.child(next, seq[start_offset + depth], bwt_pos);
      if(temp == cst.root()) { break; }
      next = temp; next_depth = cst.depth(next);
    }
    else  // Continue in the edge.
    {
      bwt_pos = cst.csa.psi[bwt_pos];
      if(bwt_pos < cst.csa.C[comp] || bwt_pos >= cst.csa.C[comp + 1]) { break; }
    }

    depth++;
    if(depth >= next_depth) { prev = next; }
  }
}

template<>
void
maximalMatch(const RelativeCST<>& cst, const sdsl::int_vector<8>& seq,
  RelativeCST<>::node_type& prev, RelativeCST<>::node_type& next,
  size_type start_offset, RelativeCST<>::size_type& depth,
  RelativeCST<>::size_type& next_depth)
{
  RelativeCST<>::size_type bwt_pos =
    (depth >= next_depth ? cst.size() : cst.index.Psi(cst.lb(next), depth - 1));

  while(start_offset + depth < seq.size() &&
    cst.forward_search(next, next_depth, depth, seq[start_offset + depth], bwt_pos))
  {
    depth++;
    if(depth >= next_depth) { prev = next; }
  }
}

template<class CST>
void
matchingStatistics(const CST& cst, const sdsl::int_vector<8>& seq,
  std::vector<range_type>& ranges, std::vector<size_type>& depths,
  const std::string& name)
{
  sdsl::util::clear(ranges); sdsl::util::clear(depths);

  // prev is the last node we have fully matched.
  // If next != prev, we are in the edge from prev to next.
  double start = readTimer();
  typename CST::node_type prev = cst.root(), next = cst.root();
  typename CST::size_type depth = 0, next_depth = 0;
  maximalMatch(cst, seq, prev, next, 0, depth, next_depth);
  ranges.push_back(range_type(cst.lb(next), cst.rb(next)));
  depths.push_back(depth);
  size_type total_length = depth;
  for(size_type i = 1; i < seq.size(); i++)
  {
    if(depth == 0)
    {
      maximalMatch(cst, seq, prev, next, i, depth, next_depth);
    }
    else
    {
      next = prev = cst.sl(prev); depth--;
      next_depth = cst.depth(prev);
      while(next_depth < depth)
      {
        typename CST::size_type bwt_pos = 0;
        next = cst.child(next, seq[i + next_depth], bwt_pos);
        next_depth = cst.depth(next);
        if(next_depth <= depth) { prev = next; }
      }
      maximalMatch(cst, seq, prev, next, i, depth, next_depth);
    }
    ranges.push_back(range_type(cst.lb(next), cst.rb(next)));
    depths.push_back(depth);
    total_length += depth;
  }
  double seconds = readTimer() - start;

  printHeader(name);
  std::cout << "Average maximal match: " << (total_length / (double)(cst.size())) << " (" << seconds << " seconds)" << std::endl;
}

//------------------------------------------------------------------------------
