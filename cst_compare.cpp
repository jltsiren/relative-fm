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

struct Timer
{
  const static size_type  INTERVAL = MILLION; // Check time after INTERVAL positions.
  constexpr static double MAX_TIME = 86400.0; // Stop after MAX_TIME seconds.

  inline static bool check(size_type i, size_type& next_check, double start_time)
  {
    if(i < next_check) { return false; }
    next_check += INTERVAL;
    return (readTimer() - start_time >= MAX_TIME);
  }

  static void print()
  {
    printHeader("Timer");
    std::cout << "interval=" << INTERVAL << ", max_time=" << MAX_TIME << std::endl;
  }
};

//------------------------------------------------------------------------------

struct MaximalMatch
{
  size_type  start, length;
  range_type range;

  MaximalMatch() {}

  MaximalMatch(size_type i, size_type depth, range_type rng) :
    start(i), length(depth), range(rng)
  {
  }

  bool operator!= (const MaximalMatch& another) const
  {
    return (this->start != another.start || this->length != another.length || this->range != another.range);
  }

  bool operator< (const MaximalMatch& another) const
  {
    return (this->start < another.start);
  }
};

std::ostream& operator<<(std::ostream& stream, const MaximalMatch& match)
{
  return stream << "(" << match.start << ", " << match.length << ", " << match.range << ")";
}

//------------------------------------------------------------------------------

template<class CST>
void buildCST(CST& cst, const std::string& base_name, const std::string& type);

void buildSelect(RelativeFM<>& fm, const std::string& base_name);

//------------------------------------------------------------------------------

template<class CST>
void forwardSearch(const CST& cst, const sdsl::int_vector<8>& seq,
  std::vector<MaximalMatch>& results, const std::string& name);

template<>
void forwardSearch(const sdsl::cst_fully<>& cst, const sdsl::int_vector<8>& seq,
  std::vector<MaximalMatch>& results, const std::string& name);

//------------------------------------------------------------------------------

template<class CST>
void backwardSearch(const CST& cst, const sdsl::int_vector<8>& seq,
  std::vector<MaximalMatch>& results, const std::string& name);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 4)
  {
    std::cerr << "Usage: cst_compare reference target sequence" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  std::cout << "Finding maximal matches with a CST" << std::endl;
  std::cout << std::endl;

  Timer::print();
  std::cout << std::endl;

  std::string ref_name = argv[1];
  printHeader("Reference"); std::cout << ref_name << std::endl;
  std::cout << std::endl;

  SimpleFM<> ref_fm(ref_name);
  size_type fm_bytes = ref_fm.reportSize();
  printSize("FM-index", fm_bytes, ref_fm.size());
  NewRelativeLCP::lcp_type ref_lcp;
  sdsl::load_from_file(ref_lcp, ref_name + LCP_EXTENSION);
  size_type lcp_bytes = sdsl::size_in_bytes(ref_lcp);
  printSize("LCP array", lcp_bytes, ref_lcp.size());
  printSize("Reference data", fm_bytes + lcp_bytes, ref_fm.size());
  std::cout << std::endl;
  std::cout << std::endl;

  std::string target_name = argv[2];
  printHeader("Target"); std::cout << target_name << std::endl;
  std::string seq_name = argv[3];
  sdsl::int_vector<8> seq;
  {
    std::ifstream in(seq_name, std::ios_base::binary);
    if(!in)
    {
      std::cerr << "cst_compare: Cannot open query file " << seq_name << std::endl;
      std::exit(EXIT_FAILURE);
    }
    size_type size = sdsl::util::file_size(seq_name); seq.resize(size);
    in.read((char*)(seq.data()), size); in.close();
  }
  printHeader("Query"); std::cout << seq_name << " (" << seq.size() << " bytes)" << std::endl;
  std::cout << std::endl;

  std::vector<MaximalMatch> rcst_results;
  {
    std::string name = "Relative (slow)";
    RelativeFM<> rfm(ref_fm, target_name);
    NewRelativeLCP rlcp(ref_lcp, target_name);
    RelativeCST<> rcst(rfm, rlcp);
    printSize(name, rcst.reportSize(), rcst.size());
    forwardSearch(rcst, seq, rcst_results, name);
    std::cout << std::endl;

    name = "Relative (fast)";
    buildSelect(rfm, target_name);
    printSize(name, rcst.reportSize(), rcst.size());
    forwardSearch(rcst, seq, rcst_results, name);
    std::cout << std::endl;
  }

  std::vector<MaximalMatch> cst_results;
  {
    std::string name = "cst_sct3_dac";
    sdsl::cst_sct3<> cst;
    buildCST(cst, target_name, name);
    forwardSearch(cst, seq, cst_results, name);
    std::cout << std::endl;
  }

  {
    std::string name = "cst_sct3_plcp";
    sdsl::cst_sct3<sdsl::csa_wt<>, sdsl::lcp_support_sada<>> cst;
    buildCST(cst, target_name, name);
    forwardSearch(cst, seq, cst_results, name);
    std::cout << std::endl;
  }

  {
    std::string name = "cst_sada";
    sdsl::cst_sada<> cst;
    buildCST(cst, target_name, name);
    forwardSearch(cst, seq, cst_results, name);
    std::cout << std::endl;
  }

  {
    std::string name = "cst_fully";
    sdsl::cst_fully<> cst;
    buildCST(cst, target_name, name);
    std::vector<MaximalMatch> fcst_results;
    forwardSearch(cst, seq, fcst_results, name);
    std::cout << std::endl;
  }

#ifdef VERIFY_RESULTS
  if(rcst_results.size() != cst_results.size())
  {
    std::cerr << "cst_compare: The number of maximal matches does not match:" << std::endl;
    std::cerr << "  RCST: " << rcst_results.size() << std::endl;
    std::cerr << "  CST:  " << cst_results.size() << std::endl;
  }
  parallelQuickSort(rcst_results.begin(), rcst_results.end());
  parallelQuickSort(cst_results.begin(), cst_results.end());
  for(size_type i = 0; i < rcst_results.size(); i++)
  {
    if(rcst_results[i] != cst_results[i])
    {
      std::cerr << "cst_compare: Maximal match " << i << ":" << std::endl;
      std::cerr << "  RCST: " << rcst_results[i] << std::endl;
      std::cerr << "  CST:  " << cst_results[i] << std::endl;
      break;
    }
  }
  std::cout << "Maximal matches verified." << std::endl;
  std::cout << std::endl;
#endif

  std::cout << "Memory usage: " << inGigabytes(memoryUsage()) << " GB" << std::endl;
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
bool
matchForward(const CST& cst, const sdsl::int_vector<8>& seq,
  typename CST::node_type& prev, typename CST::node_type& next,
  size_type start_offset, typename CST::size_type& depth,
  typename CST::size_type& next_depth)
{
  typename CST::size_type bwt_pos =
    (depth >= next_depth ? 0 : get_char_pos(cst.lb(next), depth - 1, cst.csa));

  bool match_extended = false;
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

    depth++; match_extended = true;
    if(depth >= next_depth) { prev = next; }
  }

  return match_extended;
}

/*
  RCST does not have the same members as SDSL suffix trees.
*/
template<>
bool
matchForward(const RelativeCST<>& cst, const sdsl::int_vector<8>& seq,
  RelativeCST<>::node_type& prev, RelativeCST<>::node_type& next,
  size_type start_offset, RelativeCST<>::size_type& depth,
  RelativeCST<>::size_type& next_depth)
{
  RelativeCST<>::size_type bwt_pos =
    (depth >= next_depth ? cst.size() : cst.index.Psi(cst.lb(next), depth - 1));

  bool match_extended = false;
  while(start_offset + depth < seq.size() &&
    cst.forward_search(next, next_depth, depth, seq[start_offset + depth], bwt_pos))
  {
    depth++; match_extended = true;
    if(depth >= next_depth) { prev = next; }
  }

  return match_extended;
}

/*
  cst_fully uses a different child() function.
*/
template<>
bool
matchForward(const sdsl::cst_fully<>& cst, const sdsl::int_vector<8>& seq,
  sdsl::cst_fully<>::node_type& prev, sdsl::cst_fully<>::node_type& next,
  size_type start_offset, sdsl::cst_fully<>::size_type& depth,
  sdsl::cst_fully<>::size_type& next_depth)
{
  sdsl::cst_fully<>::size_type bwt_pos =
    (depth >= next_depth ? 0 : get_char_pos(cst.lb(next), depth - 1, cst.csa));

  bool match_extended = false;
  while(start_offset + depth < seq.size())
  {
    auto comp = cst.csa.char2comp[seq[start_offset + depth]];
    if(comp == 0 && seq[start_offset + depth] != 0) { break; }
    if(depth >= next_depth) // Next node reached, follow a new edge.
    {
      if(cst.is_leaf(next)) { break; }
      sdsl::cst_fully<>::node_type temp = cst.child(next, seq[start_offset + depth], next_depth);
      if(temp == cst.root()) { break; }
      next = temp; next_depth = cst.depth(next);
      bwt_pos = get_char_pos(cst.lb(next), depth, cst.csa);
    }
    else  // Continue in the edge.
    {
      bwt_pos = cst.csa.psi[bwt_pos];
      if(bwt_pos < cst.csa.C[comp] || bwt_pos >= cst.csa.C[comp + 1]) { break; }
    }

    depth++; match_extended = true;
    if(depth >= next_depth) { prev = next; }
  }

  return match_extended;
}

//------------------------------------------------------------------------------

template<class CST>
void
forwardSearch(const CST& cst, const sdsl::int_vector<8>& seq,
  std::vector<MaximalMatch>& results, const std::string& name)
{
  sdsl::util::clear(results);

  double start = readTimer();

  /*
    'prev' is the last node we have fully matched.
    'next' is the node matching the current substring.
    If next != prev, we are in the edge from 'prev' to 'next'.
  */
  typename CST::node_type prev = cst.root(), next = cst.root();
  typename CST::size_type depth = 0, next_depth = 0;
  size_type total_length = 0;
  bool timeout = false;
  for(size_type i = 0, next_check = Timer::INTERVAL; i < seq.size(); i++)
  {
    if(depth > 0)
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
    }
    if(matchForward(cst, seq, prev, next, i, depth, next_depth))
    {
      results.push_back(MaximalMatch(i, depth, range_type(cst.lb(next), cst.rb(next))));
      total_length += depth;
    }
    if(Timer::check(i + 1, next_check, start)) { timeout = true; break; }
  }

  double seconds = readTimer() - start;
  printHeader(name);
  std::cout << results.size() << " matches of average length "
            <<( total_length / (double)(results.size()))
            << " in " << seconds << " seconds";
  if(timeout) { std::cout << " (timeout)"; }
  std::cout << std::endl;
}

/*
  cst_fully uses a different child() function.
*/
template<>
void
forwardSearch(const sdsl::cst_fully<>& cst, const sdsl::int_vector<8>& seq,
  std::vector<MaximalMatch>& results, const std::string& name)
{
  sdsl::util::clear(results);

  double start = readTimer();

  /*
    'prev' is the last node we have fully matched.
    'next' is the node matching the current substring.
    If next != prev, we are in the edge from 'prev' to 'next'.
  */
  sdsl::cst_fully<>::node_type prev = cst.root(), next = cst.root();
  sdsl::cst_fully<>::size_type depth = 0, next_depth = 0;
  size_type total_length = 0;
  bool timeout = false;
  for(size_type i = 0, next_check = Timer::INTERVAL; i < seq.size(); i++)
  {
    if(depth > 0)
    {
      next = prev = cst.sl(prev); depth--;
      next_depth = cst.depth(prev);
      while(next_depth < depth)
      {
        sdsl::cst_fully<>::size_type bwt_pos = 0;
        next = cst.child(next, seq[i + next_depth], next_depth);
        next_depth = cst.depth(next);
        if(next_depth <= depth) { prev = next; }
      }
    }
    if(matchForward(cst, seq, prev, next, i, depth, next_depth))
    {
      results.push_back(MaximalMatch(i, depth, range_type(cst.lb(next), cst.rb(next))));
      total_length += depth;
    }
    if(Timer::check(i + 1, next_check, start)) { timeout = true; break; }
  }

  double seconds = readTimer() - start;
  printHeader(name);
  std::cout << results.size() << " matches of average length "
            <<( total_length / (double)(results.size()))
            << " in " << seconds << " seconds";
  if(timeout) { std::cout << " (timeout)"; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
/*
template<class CST>
void
backwardSearch(const CST& cst, const sdsl::int_vector<8>& seq,
  std::vector<MaximalMatch>& results, const std::string& name)
{
  sdsl::util::clear(results);

  double start = readTimer();

  size_type total_length = 0;
  typename CST::node_type curr = cst.root();
  bool timeout = false;
  for(size_type i = 1, next_check = Timer::INTERVAL; i <= seq.size(); i++)
  {
    size_type pos = seq.size() - i;
    if(depth > 0)
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
    }
    if(matchForward(cst, seq, prev, next, i, depth, next_depth))
    {
      results.push_back(MaximalMatch(i, depth, range_type(cst.lb(next), cst.rb(next))));
      total_length += depth;
    }
    if(Timer::check(i, next_check, start)) { timeout = true; break; }
  }

  double seconds = readTimer() - start;
  printHeader(name);
  std::cout << results.size() << " matches of average length "
            <<( total_length / (double)(results.size()))
            << " in " << seconds << " seconds";
  if(timeout) { std::cout << " (timeout)"; }
  std::cout << std::endl;
}
*/
//------------------------------------------------------------------------------

