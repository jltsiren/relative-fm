/*
  Copyright (c) 2015, 2016, 2017 Genome Research Ltd.
  Copyright (c) 2014 Jouni Siren

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

#include "relative_fm.h"
#include "relative_lcp.h"

using namespace relative;

//------------------------------------------------------------------------------

void mainLoop(int argc, char** argv, const align_parameters& parameters, bool lcp);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  // FIXME add support for mode_ropebwt
  if(argc < 3)
  {
    std::cerr << "Usage: align_bwts [parameters] reference target1 [target2 ...]" << std::endl;

    std::cerr << "  -b N  Set BWT block size to N (default "
              << align_parameters::BLOCK_SIZE << ")" << std::endl;
    std::cerr << "  -d N  Set maximum diagonal in LCS computation to N (default "
              << align_parameters::MAX_D << ")" << std::endl;
    std::cerr << "  -l N  Partition by patterns of length up to N (default "
              << align_parameters::MAX_LENGTH << ")" << std::endl;
    std::cerr << "  -p    Preallocate buffers for LCS computation" << std::endl;

    std::cerr << "  -i    Find a BWT-invariant subsequence that supports SA/ISA samples" << std::endl;
    std::cerr << "  -L    Build also the relative LCP array" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  align_parameters parameters;
  bool lcp = false;
  int c = 0;
  while((c = getopt(argc, argv, "b:d:l:piL")) != -1)
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

  mainLoop(argc - optind, argv + optind, parameters, lcp);

  double memory = inMegabytes(memoryUsage());
  std::cout << "Memory usage: " << memory << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

void
mainLoop(int argc, char** argv, const align_parameters& parameters, bool lcp)
{
  std::string ref_name = argv[0];
  SimpleFM<bwt_type> ref(ref_name);
  ref.reportSize(true);
  RelativeLCP::lcp_type ref_lcp;
  RelativeLCP::index_type ref_index;
  if(lcp)
  {
    sdsl::load_from_file(ref_lcp, ref_name + LCP_EXTENSION);
    sdsl::load_from_file(ref_index, ref_name + DLCP_INDEX_EXTENSION);
    printSize("LCP array", sdsl::size_in_bytes(ref_lcp), ref.size()); std::cout << std::endl;
  }
  std::cout << std::endl;

  for(int arg = 1; arg < argc; arg++)
  {
    std::string seq_name = argv[arg];
    std::cout << "Target: " << seq_name << std::endl;
    SimpleFM<bwt_type> seq(seq_name);
    double start = readTimer();
    RelativeFM<bwt_type> rel(ref, seq, parameters, true);
    double seconds = readTimer() - start;
    std::cout << "Index built in " << seconds << " seconds" << std::endl;
    std::cout << std::endl;

    rel.writeTo(seq_name);
    seq.reportSize(true);
    rel.reportSize(true);

    if(lcp)
    {
      RelativeLCP::lcp_type seq_lcp;
      sdsl::load_from_file(seq_lcp, seq_name + LCP_EXTENSION);
      start = readTimer();
      RelativeLCP rlcp(ref_lcp, seq_lcp, ref_index, true);
      seconds = readTimer() - start;
      std::cout << "Relative LCP array built in " << seconds << " seconds" << std::endl;
      std::cout << std::endl;

      printSize("LCP array", sdsl::size_in_bytes(seq_lcp), seq.size()); std::cout << std::endl;
      rlcp.writeTo(seq_name);
      rlcp.reportSize(true);
    }

    std::cout << std::endl;
  }
}

//------------------------------------------------------------------------------
