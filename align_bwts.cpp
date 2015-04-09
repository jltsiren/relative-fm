#include <cstdlib>
#include <unistd.h>

#include <sdsl/rmq_support.hpp>

#include "relative_fm.h"
#include "relative_lcp.h"
#include "sequence.h"

using namespace relative;

//------------------------------------------------------------------------------

#define VERIFY_LCP
#define VERIFY_RMQ
#define VERIFY_PSV
#define VERIFY_NSV

const uint64_t LCP_QUERIES = 1000000;
const uint64_t RMQ_QUERIES = 1000000;
const uint64_t RMQ_QUERY_LENGTH = 1024;
const uint64_t PSV_NSV_QUERIES = 1000000;

//------------------------------------------------------------------------------

template<class BWTType>
void mainLoop(int argc, char** argv, const align_parameters& parameters, LoadMode mode, bool lcp);

void verifyLCP(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp);
void verifyRMQ(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp);
void verifyPSV(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp);
void verifyNSV(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  // FIXME add support for mode_ropebwt
  if(argc < 3)
  {
    std::cerr << "Usage: align_bwts [parameters] ref seq1 [seq2 ...]" << std::endl;

    std::cerr << "  -b N  Set BWT block size to N (default "
              << align_parameters::BLOCK_SIZE << ")" << std::endl;
    std::cerr << "  -d N  Set maximum diagonal in LCS computation to N (default "
              << align_parameters::MAX_D << ")" << std::endl;
    std::cerr << "  -l N  Partition by patterns of length up to N (default "
              << align_parameters::MAX_LENGTH << ")" << std::endl;
    std::cerr << "  -p    Preallocate buffers for LCS computation" << std::endl;
    std::cerr << "  -r    BWTs were built with ropebwt2" << std::endl;

    std::cerr << "  -i    Find a BWT-invariant subsequence that supports SA/ISA samples" << std::endl;
    std::cerr << "  -L    Build also the relative LCP array" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  LoadMode mode = mode_plain;
  align_parameters parameters;
  bool lcp = false;
  int c = 0;
  while((c = getopt(argc, argv, "b:d:l:priL")) != -1)
  {
    switch(c)
    {
    case 'b':
      parameters.block_size = atol(optarg); break;
    case 'd':
      parameters.max_d = atol(optarg); break;
    case 'l':
      parameters.max_length = atol(optarg); break;
    case 'p':
      parameters.preallocate = true; break;
    case 'r':
      mode = mode_ropebwt2; parameters.sorted_alphabet = false; break;
    case 'i':
      parameters.invariant = true;
      if(parameters.sa_sample_rate == align_parameters::SA_SAMPLE_RATE)
      {
        parameters.sa_sample_rate = align_parameters::SECONDARY_SA_SAMPLE_RATE;
      }
      if(parameters.isa_sample_rate == align_parameters::ISA_SAMPLE_RATE)
      {
        parameters.isa_sample_rate = align_parameters::SECONDARY_ISA_SAMPLE_RATE;
      }
      break;
    case 'L':
      lcp = true; break;
    case '?':
      return 2;
    default:
      return 3;
    }
  }

  if(lcp)
  {
    std::cout << "Relative FM-index and LCP array builder" << std::endl;
  }
  else
  {
    std::cout << "Relative FM-index builder" << std::endl;
  }
  std::cout << "Using OpenMP with " << omp_get_max_threads() << " threads" << std::endl;
  std::cout << std::endl;
  std::cout << "Algorithm: " << (parameters.invariant ? "invariant" : "partitioning") << std::endl;
  std::cout << "Input format: " << (mode == mode_ropebwt2 ? "ropebwt2" : "plain") << std::endl;
  if(parameters.sa_sample_rate != 0)
  {
    std::cout << "SA sample rate: " << parameters.sa_sample_rate << std::endl;
  }
  if(parameters.isa_sample_rate != 0)
  {
    std::cout << "ISA sample rate: " << parameters.isa_sample_rate << std::endl;
  }
  if(!(parameters.invariant))
  {
    std::cout << "Block size: " << parameters.block_size << std::endl;
    std::cout << "Maximum diagonal: " << parameters.max_d << std::endl;
    std::cout << "Maximum length: " << parameters.max_length << std::endl;
    std::cout << "Buffers: " << (parameters.preallocate ? "preallocated" : "on demand") << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Reference: " << argv[optind] << std::endl;
  std::cout << std::endl;

  if(mode == mode_ropebwt2) { mainLoop<RLSequence>(argc - optind, argv + optind, parameters, mode, lcp); }
  else { mainLoop<bwt_type>(argc - optind, argv + optind, parameters, mode, lcp); }

  double memory = inMegabytes(memoryUsage());
  std::cout << "Memory usage: " << memory << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

template<class BWTType>
void
mainLoop(int argc, char** argv, const align_parameters& parameters, LoadMode mode, bool lcp)
{
  std::string ref_name = argv[0];
  SimpleFM<BWTType> ref(ref_name, mode);
  if(mode == mode_ropebwt2) { ref.alpha.assign(ROPEBWT_ALPHABET); }
  ref.reportSize(true);
  RelativeLCP::lcp_type ref_lcp;
  RelativeLCP::index_type ref_index;
  if(lcp)
  {
    load_from_file(ref_lcp, ref_name + LCP_EXTENSION);
    load_from_file(ref_index, ref_name + DLCP_INDEX_EXTENSION);
    printSize("LCP array", size_in_bytes(ref_lcp), ref.size()); std::cout << std::endl;
  }
  std::cout << std::endl;

  for(int arg = 1; arg < argc; arg++)
  {
    std::string seq_name = argv[arg];
    std::cout << "Target: " << seq_name << std::endl;
    SimpleFM<BWTType> seq(seq_name, mode);
    if(mode == mode_ropebwt2) { seq.alpha.assign(ROPEBWT_ALPHABET); }
    double start = readTimer();
    RelativeFM<BWTType> rel(ref, seq, parameters, true);
    double seconds = readTimer() - start;
    std::cout << "Index built in " << seconds << " seconds" << std::endl;
    std::cout << std::endl;

    rel.writeTo(seq_name);
    seq.reportSize(true);
    rel.reportSize(true);

    if(lcp)
    {
      RelativeLCP::lcp_type seq_lcp;
      load_from_file(seq_lcp, seq_name + LCP_EXTENSION);
      start = readTimer();
      RelativeLCP rlcp(ref_lcp, seq_lcp, ref_index, true);
      seconds = readTimer() - start;
      std::cout << "Relative LCP array built in " << seconds << " seconds" << std::endl;
      std::cout << std::endl;

      printSize("LCP array", size_in_bytes(seq_lcp), seq.size()); std::cout << std::endl;
      rlcp.writeTo(seq_name);
      rlcp.reportSize(true);

#ifdef VERIFY_LCP
      verifyLCP(seq_lcp, rlcp);
#endif

#ifdef VERIFY_RMQ
      verifyRMQ(seq_lcp, rlcp);
#endif

#ifdef VERIFY_PSV
      verifyPSV(seq_lcp, rlcp);
#endif

#ifdef VERIFY_NSV
      verifyNSV(seq_lcp, rlcp);
#endif
    }

    std::cout << std::endl;
  }
}

//------------------------------------------------------------------------------

std::vector<uint64_t>
randomPositions(const RelativeLCP::lcp_type& lcp, uint64_t n)
{
  std::mt19937_64 rng(0xDEADBEEF);
  std::vector<uint64_t> positions(n);
  for(uint64_t i = 0; i < positions.size(); i++) { positions[i] = rng() % lcp.size(); }
  return positions;
}

void
verifyLCP(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp)
{
  uint64_t sum = 0;

  {
    std::vector<uint64_t> positions = randomPositions(lcp, LCP_QUERIES);
    double start = readTimer();
    for(uint64_t i = 0; i < positions.size(); i++)
    {
      sum += lcp[positions[i]];
    }
    double seconds = readTimer() - start;
    printTime("LCP (random)", positions.size(), seconds);
  }

  {
    double start = readTimer();
    uint64_t rank = lcp.size();
    for(uint64_t i = 0; i < lcp.size(); i++)
    {
      sum += lcp.accessForward(i, rank);
    }
    double seconds = readTimer() - start;
    printTime("LCP (seq)", lcp.size(), seconds);
  }

  {
    std::vector<uint64_t> positions = randomPositions(lcp, LCP_QUERIES);
    double start = readTimer();
    for(uint64_t i = 0; i < positions.size(); i++)
    {
      sum += rlcp[positions[i]];
    }
    double seconds = readTimer() - start;
    printTime("RLCP (random)", positions.size(), seconds);
  }

  {
    double start = readTimer();
    for(uint64_t i = 0; i < rlcp.size(); i++)
    {
      sum += rlcp[i];
    }
    double seconds = readTimer() - start;
    printTime("RLCP (seq)", rlcp.size(), seconds);
  }

  {
    double start = readTimer();
    uint64_t rank = lcp.size();
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
    ranges[i].second = std::min(lcp.size() - 1, ranges[i].first + rng() % max_length);
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
      sum += rlcp.rmq(queries[i]);
    }
    double seconds = readTimer() - start;
    printTime("RMQ", queries.size(), seconds);
  }

  {
    rmq_succinct_sct<true> rmq(&lcp);
    double start = readTimer();
    for(uint64_t i = 0; i < queries.size(); i++)
    {
      uint64_t a = rlcp.rmq(queries[i]);
      uint64_t b = rmq(queries[i].first, queries[i].second);
      if(a != b)
      {
        std::cerr << "Query " << i << ", range " << queries[i] << std::endl;
        std::cerr << "  rlcp: " << a << ", value = " << lcp[a] << std::endl;
        std::cerr << "  rmq:  " << b << ", value = " << lcp[b] << std::endl;
        break;
      }
    }
    double seconds = readTimer() - start;
    std::cout << "Verified " << queries.size() << " RMQ queries in " << seconds << " seconds" << std::endl;
  }

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

void
verifyPSV(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp)
{
  std::vector<uint64_t> queries = randomPositions(lcp, PSV_NSV_QUERIES);
  uint64_t sum = 0;

  {
    double start = readTimer();
    for(uint64_t i = 0; i < queries.size(); i++)
    {
      sum += rlcp.psv(queries[i]);
    }
    double seconds = readTimer() - start;
    printTime("PSV", queries.size(), seconds);
  }

  {
    double start = readTimer();
    bool ok = true;
    for(uint64_t i = 0; i < queries.size() && ok; i++)
    {
      uint64_t res = rlcp.psv(queries[i]), comp = lcp[queries[i]];
      if(res >= rlcp.size())
      {
        for(uint64_t j = 0; j < queries[i]; j++)
        {
          if(lcp[j] < comp)
          {
            printError(lcp, i, queries[i], j, res, true);
            ok = false; break;
          }
        }
      }
      else
      {
        if(lcp[res] >= comp)
        {
          printError(lcp, i, queries[i], res, res, true);
          ok = false; break;
        }
        for(uint64_t j = res + 1; j < queries[i]; j++)
        {
          if(lcp[j] < comp)
          {
            printError(lcp, i, queries[i], j, res, true);
            ok = false; break;
          }
        }
      }
    }
    double seconds = readTimer() - start;
    std::cout << "Verified " << queries.size() << " PSV queries in " << seconds << " seconds" << std::endl;
  }

  if(sum == 0) { std::cout << "This should not happen!" << std::endl; }
  std::cout << std::endl;
}

void
verifyNSV(const RelativeLCP::lcp_type& lcp, const RelativeLCP& rlcp)
{
  std::vector<uint64_t> queries = randomPositions(lcp, PSV_NSV_QUERIES);
  uint64_t sum = 0;

  {
    double start = readTimer();
    for(uint64_t i = 0; i < queries.size(); i++)
    {
      sum += rlcp.nsv(queries[i]);
    }
    double seconds = readTimer() - start;
    printTime("NSV", queries.size(), seconds);
  }

  {
    double start = readTimer();
    bool ok = true;
    for(uint64_t i = 0; i < queries.size() && ok; i++)
    {
      uint64_t res = rlcp.nsv(queries[i]), comp = lcp[queries[i]];
      if(res >= rlcp.size())
      {
        for(uint64_t j = queries[i] + 1; j < lcp.size(); j++)
        {
          if(lcp[j] < comp)
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
          if(lcp[j] < comp)
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

  if(sum == 0) { std::cout << "This should not happen!" << std::endl; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
