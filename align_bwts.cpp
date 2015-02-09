#include "relative_fm.h"
#include "sequence.h"


template<class BWTType>
void mainLoop(int argc, char** argv, bool ropebwt2);

int
main(int argc, char** argv)
{
  if(argc < 3)
  {
    std::cerr << "Usage: align_bwts [-r] ref seq1 [seq2 ...]" << std::endl;
    std::cerr << "  -r  BWTs were built by ropebwt2." << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  bool ropebwt2 = false;
  int ref_arg = 1;
  if(argc >= 4 && argv[1][0] == '-' && argv[1][1] == 'r')
  {
    ropebwt2 = true; ref_arg++;
  }

  std::cout << "Relative FM-index builder" << std::endl;
  std::cout << std::endl;
  std::cout << "Input format: " << (ropebwt2 ? "ropebwt2" : "plain") << std::endl;
  std::cout << "Reference: " << argv[ref_arg] << std::endl;
  std::cout << std::endl;

  if(ropebwt2) { mainLoop<bwt_type>(argc - ref_arg, argv + ref_arg, ropebwt2); }
  else         { mainLoop<RLSequence>(argc - ref_arg, argv + ref_arg, ropebwt2); }

  double memory = inMegabytes(memoryUsage());
  std::cout << "Memory usage: " << memory << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
}


template<class BWTType>
void
mainLoop(int argc, char** argv, bool ropebwt2)
{
  SimpleFM<BWTType> ref(argv[0], ropebwt2);
  if(ropebwt2) { ref.alpha.assign(ROPEBWT2_ALPHABET); }
  ref.reportSize(true);
  std::cout << std::endl;

  for(int arg = 1; arg < argc; arg++)
  {
    std::cout << "Target: " << argv[arg] << std::endl;
    SimpleFM<BWTType> seq(argv[arg]);
    if(ropebwt2) { seq.alpha.assign(ROPEBWT2_ALPHABET); }
    double start = readTimer();
    RelativeFM<SimpleFM<BWTType> > rel(ref, seq, ropebwt2, true);
    double seconds = readTimer() - start;
    std::cout << "Index built in " << seconds << " seconds" << std::endl;
    std::cout << std::endl;

    rel.writeTo(argv[arg]);
    seq.reportSize(true);
    rel.reportSize(true);
    std::cout << std::endl;
  }
}
