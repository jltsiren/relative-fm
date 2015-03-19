#include <cstdlib>
#include <unistd.h>

#include "relative_fm.h"
#include "sequence.h"

using namespace relative;

template<class BWTType>
void mainLoop(int argc, char** argv, const align_parameters& parameters, LoadMode mode);


int
main(int argc, char** argv)
{
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
    std::cerr << std::endl;
    return 1;
  }

  LoadMode mode = mode_plain;
  align_parameters parameters;
  int c = 0;
  while((c = getopt(argc, argv, "b:d:l:pri")) != -1)
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
    case '?':
      return 2;
    default:
      return 3;
    }
  }

  std::cout << "Relative FM-index builder" << std::endl;
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

  if(mode == mode_ropebwt2) { mainLoop<RLSequence>(argc - optind, argv + optind, parameters, mode); }
  else { mainLoop<bwt_type>(argc - optind, argv + optind, parameters, mode); }

  double memory = inMegabytes(memoryUsage());
  std::cout << "Memory usage: " << memory << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
}


template<class BWTType>
void
mainLoop(int argc, char** argv, const align_parameters& parameters, LoadMode mode)
{
  SimpleFM<BWTType> ref(argv[0], mode);
  if(mode == mode_ropebwt2) { ref.alpha.assign(ROPEBWT2_ALPHABET); }
  ref.reportSize(true);
  std::cout << std::endl;

  for(int arg = 1; arg < argc; arg++)
  {
    std::cout << "Target: " << argv[arg] << std::endl;
    SimpleFM<BWTType> seq(argv[arg], mode);
    if(mode == mode_ropebwt2) { seq.alpha.assign(ROPEBWT2_ALPHABET); }
    double start = readTimer();
    RelativeFM<BWTType> rel(ref, seq, parameters, true);
    double seconds = readTimer() - start;
    std::cout << "Index built in " << seconds << " seconds" << std::endl;
    std::cout << std::endl;

    rel.writeTo(argv[arg]);
    seq.reportSize(true);
    rel.reportSize(true);
    std::cout << std::endl;
  }
}
