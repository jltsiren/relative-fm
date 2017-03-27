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

#include <cstdlib>
#include <unistd.h>

#include <sdsl/rmq_support.hpp>

#include "relative_fm.h"
#include "new_relative_lcp.h"

using namespace relative;

//------------------------------------------------------------------------------

// This reads patterns from file 'patterns'.
//#define VERIFY_FORWARD_SEARCH

//#define VERIFY_LF
//#define VERIFY_PSI

//#define VERIFY_LCP
#define VERIFY_RMQ
//#define VERIFY_PSV
//#define VERIFY_PSEV
//#define VERIFY_NSV
//#define VERIFY_NSEV

// Verify the queries or just run the speed tests.
#define VERIFY_QUERIES

const size_type LF_QUERIES = 10 * MILLION;
const size_type PSI_QUERIES = 10 * MILLION;

const size_type LCP_QUERIES = 100 * MILLION;
const size_type RMQ_QUERIES = 100 * MILLION;
const size_type RMQ_QUERY_LENGTH = 16;
const size_type PSV_NSV_QUERIES = 100 * MILLION;

//------------------------------------------------------------------------------

void verifyForwardSearch(const SimpleFM<>& fm, const RelativeFM<>& rfm, const NewRelativeLCP& lcp);

template<class Index>
void verifyLF(const Index& index, const std::string& type);

template<class Index>
void verifyPsi(const Index& index, const std::string& type);

void buildSelect(RelativeFM<>& fm, const std::string& base_name);

void verifyLCP(const NewRelativeLCP::lcp_type& lcp, const NewRelativeLCP& rlcp);
void verifyRMQ(const NewRelativeLCP::lcp_type& lcp, const NewRelativeLCP& rlcp);
void verifyPSV(const NewRelativeLCP::lcp_type& lcp, const NewRelativeLCP& rlcp);
void verifyPSEV(const NewRelativeLCP::lcp_type& lcp, const NewRelativeLCP& rlcp);
void verifyNSV(const NewRelativeLCP::lcp_type& lcp, const NewRelativeLCP& rlcp);
void verifyNSEV(const NewRelativeLCP::lcp_type& lcp, const NewRelativeLCP& rlcp);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 3)
  {
    std::cerr << "Usage: verify ref seq1 [seq2 ...]" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  std::cout << "RFM and RLCP verifier" << std::endl;
  std::cout << std::endl;

  std::string ref_name = argv[1];
  std::cout << "Reference: " << ref_name << std::endl;
  std::cout << std::endl;

  SimpleFM<> ref_fm(ref_name);
  printSize("FM-index", ref_fm.reportSize(), ref_fm.size()); std::cout << std::endl;

  NewRelativeLCP::lcp_type ref_lcp;
  sdsl::load_from_file(ref_lcp, ref_name + LCP_EXTENSION);
  printSize("LCP array", sdsl::size_in_bytes(ref_lcp), ref_lcp.size()); std::cout << std::endl;

  std::cout << std::endl;

  for(int arg = 2; arg < argc; arg++)
  {
    std::string seq_name = argv[arg];
    std::cout << "Target: " << seq_name << std::endl;
    std::cout << std::endl;

    SimpleFM<> seq_fm(seq_name);
    seq_fm.reportSize(true); std::cout << std::endl;

    RelativeFM<> rfm(ref_fm, seq_name);
    rfm.reportSize(true); std::cout << std::endl;

    NewRelativeLCP::lcp_type seq_lcp;
    sdsl::load_from_file(seq_lcp, seq_name + LCP_EXTENSION);
    printSize("LCP array", sdsl::size_in_bytes(seq_lcp), seq_lcp.size()); std::cout << std::endl;

    NewRelativeLCP rlcp(ref_lcp, seq_name);
    rlcp.reportSize(true); std::cout << std::endl;

#ifdef VERIFY_FORWARD_SEARCH
    buildSelect(rfm, seq_name);
    verifyForwardSearch(seq_fm, rfm, rlcp);
#endif

#ifdef VERIFY_LF
    verifyLF(seq_fm, "LF (FM)");
    verifyLF(rfm, "LF (RFM)");
#endif

#ifdef VERIFY_PSI
    verifyPsi(seq_fm, "Psi (FM)");
    verifyPsi(rfm, "Psi (RFM, slow)");
    buildSelect(rfm, seq_name);
    verifyPsi(rfm, "Psi (RFM, fast)");
#endif

#ifdef VERIFY_LCP
    verifyLCP(seq_lcp, rlcp);
#endif

#ifdef VERIFY_RMQ
    verifyRMQ(seq_lcp, rlcp);
#endif

#ifdef VERIFY_PSV
    verifyPSV(seq_lcp, rlcp);
#endif

#ifdef VERIFY_PSEV
    verifyPSEV(seq_lcp, rlcp);
#endif

#ifdef VERIFY_NSV
    verifyNSV(seq_lcp, rlcp);
#endif

#ifdef VERIFY_NSEV
    verifyNSEV(seq_lcp, rlcp);
#endif

    std::cout << std::endl;
  }

  std::cout << "Memory used: " << inMegabytes(memoryUsage()) << " MB" << std::endl;
  std::cout << std::endl;
  return 0;
}

//------------------------------------------------------------------------------

/*

template<class Index>
void
findQueries(const Index& index, const std::string& type,
  const std::vector<std::string>& patterns, std::vector<range_type>& results)
{
  double start = readTimer();
  for(size_type i = 0; i < patterns.size(); i++)
  {
    results[i] = index.find(patterns[i].begin(), patterns[i].end());
  }
  double seconds = readTimer() - start;
  printTime(type, patterns.size(), seconds);
}

void
verifyForwardSearch(const SimpleFM<>& fm, const RelativeFM<>& rfm, const NewRelativeLCP& lcp)
{
  std::vector<std::string> patterns;
  size_type total_length = readRows("patterns", patterns, true);
  std::cout << "Read " << patterns.size() << " patterns of total length " << total_length << std::endl;

  std::vector<range_type> fm_results(patterns.size());
  findQueries(fm, "SimpleFM", patterns, fm_results);

  std::vector<range_type> cst_results(patterns.size());
  RelativeCST<> rcst(rfm, lcp);
  findQueries(rcst, "Relative CST", patterns, cst_results);

  {
    bool ok = true;
    for(size_type i = 0; i < patterns.size(); i++)
    {
      if(fm_results[i] != cst_results[i])
      {
        std::cerr << "verify: Query " << i << ", pattern " << patterns[i] << std::endl;
        std::cerr << "  SimpleFM: " << fm_results[i] << std::endl;
        std::cerr << "  Relative CST: " << cst_results[i] << std::endl;
        ok = false; break;
      }
    }
    if(ok) { std::cout << "Results successfully verified" << std::endl; }
  }

  std::cout << std::endl;
}

*/

//------------------------------------------------------------------------------

std::vector<size_type>
randomPositions(size_type n, size_type size)
{
  std::mt19937_64 rng(0xDEADBEEF);
  std::vector<size_type> positions(n);
  for(size_type i = 0; i < positions.size(); i++) { positions[i] = rng() % size; }
  return positions;
}

template<class ArrayType>
size_type
timeQueries(const ArrayType& array, const std::vector<size_type>& queries, const std::string& name)
{
  double start = readTimer();
  size_type sum = 0;
  for(size_type i = 0; i < queries.size(); i++) { sum += array[queries[i]]; }
  double seconds = readTimer() - start;
  printTime(name, queries.size(), seconds);
  return sum;
}

template<class Iterator>
size_type
timeQueries(Iterator begin, Iterator end, const std::string& name)
{
  double start = readTimer();
  size_type sum = 0, size = end - begin;
  while(begin != end) { sum += *begin; ++begin; }
  double seconds = readTimer() - start;
  printTime(name, size, seconds);
  return sum;
}

template<class ArrayType>
size_type
timeQueries(const ArrayType& array, const std::string& name)
{
  double start = readTimer();
  size_type sum = 0;
  for(size_type i = 0; i < array.size(); i++) { sum += array[i]; }
  double seconds = readTimer() - start;
  printTime(name, array.size(), seconds);
  return sum;
}

//------------------------------------------------------------------------------

template<class Index>
void
verifyLF(const Index& index, const std::string& type)
{
  size_type sum = 0;
  std::vector<size_type> positions = randomPositions(LF_QUERIES, index.size());

  {
    double start = readTimer();
    for(size_type i = 0; i < positions.size(); i++)
    {
      sum += index.LF(positions[i]).first;
    }
    double seconds = readTimer() - start;
    printTime(type, positions.size(), seconds);
  }

  if(sum == 0) { std::cout << "This should not happen!" << std::endl; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

template<class Index>
void
verifyPsi(const Index& index, const std::string& type)
{
  size_type sum = 0;
  std::vector<size_type> positions = randomPositions(PSI_QUERIES, index.size());

  {
    double start = readTimer();
    for(size_type i = 0; i < positions.size(); i++)
    {
      sum += index.Psi(positions[i]);
    }
    double seconds = readTimer() - start;
    printTime(type, positions.size(), seconds);
  }

#ifdef VERIFY_QUERIES
  {
    double start = readTimer();
    for(size_type i = index.sequences(); i < index.size(); i++)
    {
      size_type next_pos = index.Psi(i);
      size_type prev_pos = index.LF(next_pos).first;
      if(prev_pos != i)
      {
        std::cerr << "verify: Psi(" << i << ") = " << next_pos
                  << ", LF(Psi(" << i << ")) = " << prev_pos << std::endl;
        std::cerr << "verify: SA[" << i << "] = " << index.locate(i)
                  << ", SA[" << (next_pos - 1) << "] = " << index.locate(next_pos - 1)
                  << ", SA[" << next_pos << "] = " << index.locate(next_pos)
                  << ", SA[" << (next_pos + 1) << "] = " << index.locate(next_pos + 1) << std::endl;
        break;
      }
    }
    double seconds = readTimer() - start;
    std::cout << "Psi verified in " << seconds << " seconds" << std::endl;
  }
#endif

  if(sum == 0) { std::cout << "This should not happen!" << std::endl; }
  std::cout << std::endl;
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

void
verifyLCP(const NewRelativeLCP::lcp_type& lcp, const NewRelativeLCP& rlcp)
{
  size_type sum = 0;
  std::vector<size_type> positions = randomPositions(LCP_QUERIES, lcp.size());

  sum += timeQueries(lcp, positions, "LCP (random)");

  sum += timeQueries(lcp.begin(), lcp.end(), "LCP (seq)");

  sum += timeQueries(rlcp, positions, "RLCP (random)");

  sum += timeQueries(rlcp, "RLCP (seq)");

#ifdef VERIFY_QUERIES
  {
    double start = readTimer();
    NewRelativeLCP::lcp_type::iterator curr = lcp.begin();
    for(size_type i = 0; i < rlcp.size(); i++)
    {
      if(rlcp[i] != *curr)
      {
        std::cerr << "rlcp[" << i << "] = " << rlcp[i] << ", lcp[" << i << "] = " << *curr << std::endl;
        break;
      }
      ++curr;
    }
    double seconds = readTimer() - start;
    std::cout << "Relative LCP array verified in " << seconds << " seconds" << std::endl;
  }
#endif

  if(sum == 0) { std::cout << "This should not happen!" << std::endl; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

std::vector<range_type>
randomRanges(const NewRelativeLCP::lcp_type& lcp, size_type n, size_type max_length)
{
  std::mt19937_64 rng(0xDEADBEEF);
  std::vector<range_type> ranges(n);
  for(size_type i = 0; i < ranges.size(); i++)
  {
    ranges[i].first = rng() % lcp.size();
    size_type len = max_length;
    while(len < lcp.size() && (rng() & 1)) { len *= max_length; }
    ranges[i].second = std::min(lcp.size() - 1, ranges[i].first + rng() % len);
  }
  return ranges;
}

void
verifyRMQ(const NewRelativeLCP::lcp_type& lcp, const NewRelativeLCP& rlcp)
{
  std::vector<range_type> queries = randomRanges(lcp, RMQ_QUERIES, RMQ_QUERY_LENGTH);
  size_type sum = 0;

  {
    double start = readTimer();
    for(size_type i = 0; i < queries.size(); i++)
    {
      sum += rlcp.rmq(queries[i]).first;
    }
    double seconds = readTimer() - start;
    printTime("RMQ", queries.size(), seconds);
  }

#ifdef VERIFY_QUERIES
  {
    sdsl::rmq_succinct_sct<true> rmq(&lcp);
    double start = readTimer();
    for(size_type i = 0; i < queries.size(); i++)
    {
      range_type a = rlcp.rmq(queries[i]);
      size_type b = rmq(queries[i].first, queries[i].second);
      if(a.first != b)
      {
        std::cerr << "Query " << i << ", range " << queries[i] << std::endl;
        std::cerr << "  rlcp: " << a.first << ", value = " << a.second << std::endl;
        std::cerr << "  rmq:  " << b << ", value = " << lcp[b] << std::endl;
        break;
      }
    }
    double seconds = readTimer() - start;
    std::cout << "Verified " << queries.size() << " RMQ queries in " << seconds << " seconds" << std::endl;
  }
#endif

  if(sum == 0) { std::cout << "This should not happen!" << std::endl; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

void
printLCP(const NewRelativeLCP::lcp_type& lcp, size_type pos)
{
  std::cerr << "  lcp[" << pos << "] = " << lcp[pos] << std::endl;
}

void
printError(const NewRelativeLCP::lcp_type& lcp, size_type query, size_type query_pos,
  size_type error_pos, size_type val, bool psv)
{
  std::cerr << "Query " << query << ", position " << query_pos << std::endl;
  printLCP(lcp, query_pos);
  if(psv) { std::cerr << "  psv("; } else { std::cerr << "  nsv("; }
  std::cerr << query_pos << ") = ";
  if(val >= lcp.size()) { std::cerr << "null"; } else { std::cerr << val; }
  std::cerr << std::endl;
  printLCP(lcp, error_pos);
}

//------------------------------------------------------------------------------

/*

void
verifyPSV(const NewRelativeLCP::lcp_type& lcp, const NewRelativeLCP& rlcp)
{
  std::vector<size_type> queries = randomPositions(PSV_NSV_QUERIES, lcp.size());
  size_type sum = 0;

  {
    double start = readTimer();
    for(size_type i = 0; i < queries.size(); i++)
    {
      sum += rlcp.psv(queries[i]).first;
    }
    double seconds = readTimer() - start;
    printTime("PSV", queries.size(), seconds);
  }

#ifdef VERIFY_QUERIES
  {
    double start = readTimer();
    bool ok = true;
    for(size_type i = 0; i < queries.size() && ok; i++)
    {
      size_type rank = lcp.initBackward(queries[i]);
      size_type res = rlcp.psv(queries[i]).first, comp = lcp.accessBackward(queries[i], rank);
      if(res >= rlcp.size())
      {
        for(size_type j = queries[i]; j > 0; j--)
        {
          if(lcp.accessBackward(j - 1, rank) < comp)
          {
            printError(lcp, i, queries[i], j - 1, res, true);
            ok = false; break;
          }
        }
      }
      else
      {
        for(size_type j = queries[i] - 1; j > res; j--)
        {
          if(lcp.accessBackward(j, rank) < comp)
          {
            printError(lcp, i, queries[i], j, res, true);
            ok = false; break;
          }
        }
        if(lcp.accessBackward(res, rank) >= comp)
        {
          printError(lcp, i, queries[i], res, res, true);
          ok = false; break;
        }
      }
    }
    double seconds = readTimer() - start;
    std::cout << "Verified " << queries.size() << " PSV queries in " << seconds << " seconds" << std::endl;
  }
#endif

  if(sum == 0) { std::cout << "This should not happen!" << std::endl; }
  std::cout << std::endl;
}

*/

//------------------------------------------------------------------------------

/*

void
verifyPSEV(const NewRelativeLCP::lcp_type& lcp, const NewRelativeLCP& rlcp)
{
  std::vector<size_type> queries = randomPositions(PSV_NSV_QUERIES, lcp.size());
  size_type sum = 0;

  {
    double start = readTimer();
    for(size_type i = 0; i < queries.size(); i++)
    {
      sum += rlcp.psev(queries[i]).first;
    }
    double seconds = readTimer() - start;
    printTime("PSEV", queries.size(), seconds);
  }

#ifdef VERIFY_QUERIES
  {
    double start = readTimer();
    bool ok = true;
    for(size_type i = 0; i < queries.size() && ok; i++)
    {
      size_type rank = lcp.initBackward(queries[i]);
      size_type res = rlcp.psev(queries[i]).first, comp = lcp.accessBackward(queries[i], rank);
      if(res >= rlcp.size())
      {
        for(size_type j = queries[i]; j > 0; j--)
        {
          if(lcp.accessBackward(j - 1, rank) <= comp)
          {
            printError(lcp, i, queries[i], j - 1, res, true);
            ok = false; break;
          }
        }
      }
      else
      {
        for(size_type j = queries[i] - 1; j > res; j--)
        {
          if(lcp.accessBackward(j, rank) <= comp)
          {
            printError(lcp, i, queries[i], j, res, true);
            ok = false; break;
          }
        }
        if(lcp.accessBackward(res, rank) > comp)
        {
          printError(lcp, i, queries[i], res, res, true);
          ok = false; break;
        }
      }
    }
    double seconds = readTimer() - start;
    std::cout << "Verified " << queries.size() << " PSEV queries in " << seconds << " seconds" << std::endl;
  }
#endif

  if(sum == 0) { std::cout << "This should not happen!" << std::endl; }
  std::cout << std::endl;
}

*/

//------------------------------------------------------------------------------

/*

void
verifyNSV(const NewRelativeLCP::lcp_type& lcp, const NewRelativeLCP& rlcp)
{
  std::vector<size_type> queries = randomPositions(PSV_NSV_QUERIES, lcp.size());
  size_type sum = 0;

  {
    double start = readTimer();
    for(size_type i = 0; i < queries.size(); i++)
    {
      sum += rlcp.nsv(queries[i]).first;
    }
    double seconds = readTimer() - start;
    printTime("NSV", queries.size(), seconds);
  }

#ifdef VERIFY_QUERIES
  {
    double start = readTimer();
    bool ok = true;
    for(size_type i = 0; i < queries.size() && ok; i++)
    {
      size_type rank = lcp.initForward(queries[i]);
      size_type res = rlcp.nsv(queries[i]).first, comp = lcp.accessForward(queries[i], rank);
      if(res >= rlcp.size())
      {
        for(size_type j = queries[i] + 1; j < lcp.size(); j++)
        {
          if(lcp.accessForward(j, rank) < comp)
          {
            printError(lcp, i, queries[i], j, res, false);
            ok = false; break;
          }
        }
      }
      else
      {
        if(lcp[res] >= comp)
        {
          printError(lcp, i, queries[i], res, res, false);
          ok = false; break;
        }
        for(size_type j = queries[i] + 1; j < res; j++)
        {
          if(lcp.accessForward(j, rank) < comp)
          {
            printError(lcp, i, queries[i], j, res, false);
            ok = false; break;
          }
        }
      }
    }
    double seconds = readTimer() - start;
    std::cout << "Verified " << queries.size() << " NSV queries in " << seconds << " seconds" << std::endl;
  }
#endif

  if(sum == 0) { std::cout << "This should not happen!" << std::endl; }
  std::cout << std::endl;
}

*/

//------------------------------------------------------------------------------

/*

void
verifyNSEV(const NewRelativeLCP::lcp_type& lcp, const NewRelativeLCP& rlcp)
{
  std::vector<size_type> queries = randomPositions(PSV_NSV_QUERIES, lcp.size());
  size_type sum = 0;

  {
    double start = readTimer();
    for(size_type i = 0; i < queries.size(); i++)
    {
      sum += rlcp.nsev(queries[i]).first;
    }
    double seconds = readTimer() - start;
    printTime("NSEV", queries.size(), seconds);
  }

#ifdef VERIFY_QUERIES
  {
    double start = readTimer();
    bool ok = true;
    for(size_type i = 0; i < queries.size() && ok; i++)
    {
      size_type rank = lcp.initForward(queries[i]);
      size_type res = rlcp.nsev(queries[i]).first, comp = lcp.accessForward(queries[i], rank);
      if(res >= rlcp.size())
      {
        for(size_type j = queries[i] + 1; j < lcp.size(); j++)
        {
          if(lcp.accessForward(j, rank) <= comp)
          {
            printError(lcp, i, queries[i], j, res, false);
            ok = false; break;
          }
        }
      }
      else
      {
        if(lcp[res] > comp)
        {
          printError(lcp, i, queries[i], res, res, false);
          ok = false; break;
        }
        for(size_type j = queries[i] + 1; j < res; j++)
        {
          if(lcp.accessForward(j, rank) <= comp)
          {
            printError(lcp, i, queries[i], j, res, false);
            ok = false; break;
          }
        }
      }
    }
    double seconds = readTimer() - start;
    std::cout << "Verified " << queries.size() << " NSEV queries in " << seconds << " seconds" << std::endl;
  }
#endif

  if(sum == 0) { std::cout << "This should not happen!" << std::endl; }
  std::cout << std::endl;
}

*/

//------------------------------------------------------------------------------
