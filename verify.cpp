#include <cstdlib>
#include <unistd.h>

#include <sdsl/rmq_support.hpp>

#include "relative_cst.h"

using namespace relative;

//------------------------------------------------------------------------------

// This reads patterns from file 'patterns'.
#define VERIFY_FORWARD_SEARCH

//#define VERIFY_LF
//#define VERIFY_PSI

//#define VERIFY_LCP
//#define VERIFY_RMQ
//#define VERIFY_PSV
//#define VERIFY_PSEV
//#define VERIFY_NSV
//#define VERIFY_NSEV

// Verify the queries or just run the speed tests.
//#define VERIFY_QUERIES

const uint64_t MILLION = 1000000;

const uint64_t LF_QUERIES = 10 * MILLION;
const uint64_t PSI_QUERIES = 10 * MILLION;

const uint64_t LCP_QUERIES = 100 * MILLION;
const uint64_t RMQ_QUERIES = 100 * MILLION;
const uint64_t RMQ_QUERY_LENGTH = 16;
const uint64_t PSV_NSV_QUERIES = 100 * MILLION;

//------------------------------------------------------------------------------

void verifyForwardSearch(const SimpleFM<>& fm, const RelativeFM<>& rfm, const RelativeLCP& lcp);

template<class Index>
void verifyLF(const Index& index, const std::string& type);

template<class Index>
void verifyPsi(const Index& index, const std::string& type);

void buildSelect(RelativeFM<>& fm, const std::string& base_name);

void verifyLCP(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp);
void verifyRMQ(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp);
void verifyPSV(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp);
void verifyPSEV(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp);
void verifyNSV(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp);
void verifyNSEV(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp);

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

  RelativeLCP::lcp_type ref_lcp;
  load_from_file(ref_lcp, ref_name + LCP_EXTENSION);
  printSize("LCP array", size_in_bytes(ref_lcp), ref_lcp.size()); std::cout << std::endl;

  std::cout << std::endl;

  for(int arg = 2; arg < argc; arg++)
  {
    std::string seq_name = argv[arg];
    std::cout << "Target: " << seq_name << std::endl;
    std::cout << std::endl;

    SimpleFM<> seq_fm(seq_name);
    seq_fm.reportSize(true);

    RelativeFM<> rfm(ref_fm, seq_name);
    rfm.reportSize(true);

    RelativeLCP::lcp_type seq_lcp;
    load_from_file(seq_lcp, seq_name + LCP_EXTENSION);
    printSize("LCP array", size_in_bytes(seq_lcp), seq_lcp.size()); std::cout << std::endl;

    RelativeLCP rlcp(ref_lcp, seq_name);
    rlcp.reportSize(true);

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

template<class Index>
void
findQueries(const Index& index, const std::string& type,
  const std::vector<std::string>& patterns, std::vector<range_type>& results)
{
  double start = readTimer();
  for(uint64_t i = 0; i < patterns.size(); i++)
  {
    results[i] = index.find(patterns[i].begin(), patterns[i].end());
  }
  double seconds = readTimer() - start;
  printTime(type, patterns.size(), seconds);
}

void
verifyForwardSearch(const SimpleFM<>& fm, const RelativeFM<>& rfm, const RelativeLCP& lcp)
{
  std::vector<std::string> patterns;
  uint64_t total_length = readRows("patterns", patterns, true);
  std::cout << "Read " << patterns.size() << " patterns of total length " << total_length << std::endl;

  std::vector<range_type> fm_results(patterns.size());
  findQueries(fm, "SimpleFM", patterns, fm_results);

  std::vector<range_type> cst_results(patterns.size());
  RelativeCST<> rcst(rfm, lcp);
  findQueries(rcst, "Relative CST", patterns, cst_results);

  {
    bool ok = true;
    for(uint64_t i = 0; i < patterns.size(); i++)
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

//------------------------------------------------------------------------------

std::vector<uint64_t>
randomPositions(uint64_t n, uint64_t size)
{
  std::mt19937_64 rng(0xDEADBEEF);
  std::vector<uint64_t> positions(n);
  for(uint64_t i = 0; i < positions.size(); i++) { positions[i] = rng() % size; }
  return positions;
}

template<class ArrayType>
uint64_t
timeQueries(const ArrayType& array, const std::vector<uint64_t>& queries, std::string name)
{
  double start = readTimer();
  uint64_t sum = 0;
  for(uint64_t i = 0; i < queries.size(); i++) { sum += array[queries[i]]; }
  double seconds = readTimer() - start;
  printTime(name, queries.size(), seconds);
  return sum;
}

template<class ArrayType>
uint64_t
timeQueries(const ArrayType& array, std::string name)
{
  double start = readTimer();
  uint64_t sum = 0;
  for(uint64_t i = 0; i < array.size(); i++) { sum += array[i]; }
  double seconds = readTimer() - start;
  printTime(name, array.size(), seconds);
  return sum;
}

//------------------------------------------------------------------------------

template<class Index>
void
verifyLF(const Index& index, const std::string& type)
{
  uint64_t sum = 0;
  std::vector<uint64_t> positions = randomPositions(LF_QUERIES, index.size());

  {
    double start = readTimer();
    for(uint64_t i = 0; i < positions.size(); i++)
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
  uint64_t sum = 0;
  std::vector<uint64_t> positions = randomPositions(PSI_QUERIES, index.size());

  {
    double start = readTimer();
    for(uint64_t i = 0; i < positions.size(); i++)
    {
      sum += index.Psi(positions[i]);
    }
    double seconds = readTimer() - start;
    printTime(type, positions.size(), seconds);
  }

#ifdef VERIFY_QUERIES
  {
    double start = readTimer();
    for(uint64_t i = index.sequences(); i < index.size(); i++)
    {
      uint64_t next_pos = index.Psi(i);
      uint64_t prev_pos = index.LF(next_pos).first;
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
verifyLCP(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp)
{
  uint64_t sum = 0;
  std::vector<uint64_t> positions = randomPositions(LCP_QUERIES, lcp.size());

  sum += timeQueries(lcp, positions, "LCP (random)");

  {
    double start = readTimer();
    uint64_t rank = lcp.initForward(0);
    for(uint64_t i = 0; i < lcp.size(); i++)
    {
      sum += lcp.accessForward(i, rank);
    }
    double seconds = readTimer() - start;
    printTime("LCP (seq)", lcp.size(), seconds);
  }

  sum += timeQueries(rlcp, positions, "RLCP (random)");

  sum += timeQueries(rlcp, "RLCP (seq)");

#ifdef VERIFY_QUERIES
  {
    double start = readTimer();
    uint64_t rank = lcp.initForward(0);
    for(uint64_t i = 0; i < rlcp.size(); i++)
    {
      if(rlcp[i] != lcp.accessForward(i, rank))
      {
        std::cerr << "rlcp[" << i << "] = " << rlcp[i] << ", lcp[" << i << "] = " << lcp[i] << std::endl;
        break;
      }
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
randomRanges(const RelativeLCP::lcp_type& lcp, uint64_t n, uint64_t max_length)
{
  std::mt19937_64 rng(0xDEADBEEF);
  std::vector<range_type> ranges(n);
  for(uint64_t i = 0; i < ranges.size(); i++)
  {
    ranges[i].first = rng() % lcp.size();
    uint64_t len = max_length;
    while(len < lcp.size() && (rng() & 1)) { len *= max_length; }
    ranges[i].second = std::min(lcp.size() - 1, ranges[i].first + rng() % len);
  }
  return ranges;
}

void
verifyRMQ(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp)
{
  std::vector<range_type> queries = randomRanges(lcp, RMQ_QUERIES, RMQ_QUERY_LENGTH);
  uint64_t sum = 0;

  {
    double start = readTimer();
    for(uint64_t i = 0; i < queries.size(); i++)
    {
      sum += rlcp.rmq(queries[i]).first;
    }
    double seconds = readTimer() - start;
    printTime("RMQ", queries.size(), seconds);
  }

#ifdef VERIFY_QUERIES
  {
    rmq_succinct_sct<true> rmq(&lcp);
    double start = readTimer();
    for(uint64_t i = 0; i < queries.size(); i++)
    {
      range_type a = rlcp.rmq(queries[i]);
      uint64_t b = rmq(queries[i].first, queries[i].second);
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
printLCP(const RelativeLCP::lcp_type& lcp, uint64_t pos)
{
  std::cerr << "  lcp[" << pos << "] = " << lcp[pos] << std::endl;
}

void
printError(const RelativeLCP::lcp_type& lcp, uint64_t query, uint64_t query_pos,
  uint64_t error_pos, uint64_t val, bool psv)
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

void
verifyPSV(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp)
{
  std::vector<uint64_t> queries = randomPositions(PSV_NSV_QUERIES, lcp.size());
  uint64_t sum = 0;

  {
    double start = readTimer();
    for(uint64_t i = 0; i < queries.size(); i++)
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
    for(uint64_t i = 0; i < queries.size() && ok; i++)
    {
      uint64_t rank = lcp.initBackward(queries[i]);
      uint64_t res = rlcp.psv(queries[i]).first, comp = lcp.accessBackward(queries[i], rank);
      if(res >= rlcp.size())
      {
        for(uint64_t j = queries[i]; j > 0; j--)
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
        for(uint64_t j = queries[i] - 1; j > res; j--)
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

//------------------------------------------------------------------------------

void
verifyPSEV(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp)
{
  std::vector<uint64_t> queries = randomPositions(PSV_NSV_QUERIES, lcp.size());
  uint64_t sum = 0;

  {
    double start = readTimer();
    for(uint64_t i = 0; i < queries.size(); i++)
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
    for(uint64_t i = 0; i < queries.size() && ok; i++)
    {
      uint64_t rank = lcp.initBackward(queries[i]);
      uint64_t res = rlcp.psev(queries[i]).first, comp = lcp.accessBackward(queries[i], rank);
      if(res >= rlcp.size())
      {
        for(uint64_t j = queries[i]; j > 0; j--)
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
        for(uint64_t j = queries[i] - 1; j > res; j--)
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

//------------------------------------------------------------------------------

void
verifyNSV(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp)
{
  std::vector<uint64_t> queries = randomPositions(PSV_NSV_QUERIES, lcp.size());
  uint64_t sum = 0;

  {
    double start = readTimer();
    for(uint64_t i = 0; i < queries.size(); i++)
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
    for(uint64_t i = 0; i < queries.size() && ok; i++)
    {
      uint64_t rank = lcp.initForward(queries[i]);
      uint64_t res = rlcp.nsv(queries[i]).first, comp = lcp.accessForward(queries[i], rank);
      if(res >= rlcp.size())
      {
        for(uint64_t j = queries[i] + 1; j < lcp.size(); j++)
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
        for(uint64_t j = queries[i] + 1; j < res; j++)
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

//------------------------------------------------------------------------------

void
verifyNSEV(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp)
{
  std::vector<uint64_t> queries = randomPositions(PSV_NSV_QUERIES, lcp.size());
  uint64_t sum = 0;

  {
    double start = readTimer();
    for(uint64_t i = 0; i < queries.size(); i++)
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
    for(uint64_t i = 0; i < queries.size() && ok; i++)
    {
      uint64_t rank = lcp.initForward(queries[i]);
      uint64_t res = rlcp.nsev(queries[i]).first, comp = lcp.accessForward(queries[i], rank);
      if(res >= rlcp.size())
      {
        for(uint64_t j = queries[i] + 1; j < lcp.size(); j++)
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
        for(uint64_t j = queries[i] + 1; j < res; j++)
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

//------------------------------------------------------------------------------
