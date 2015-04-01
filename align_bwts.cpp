#include <cstdlib>
#include <unistd.h>

#include "relative_fm.h"
#include "relative_lcp.h"
#include "sequence.h"

using namespace relative;

template<class BWTType>
void mainLoop(int argc, char** argv, const align_parameters& parameters, LoadMode mode, bool lcp);


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

/*      // Verify the correctness of the relative LCP array.
      start = readTimer();
      for(uint64_t i = 0; i < seq_lcp.size(); i++)
      {
        if(rlcp[i] != seq_lcp[i])
        {
          std::cerr << "rlcp[" << i << "] = " << rlcp[i] << ", seq_lcp[" << i << "] = " << seq_lcp[i] << std::endl;
          break;
        }
      }
      seconds = readTimer() - start;
      std::cout << "Relative LCP array verified in " << seconds << " seconds" << std::endl;*/
    }

    std::cout << std::endl;
  }
}
