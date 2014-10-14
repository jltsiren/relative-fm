#include "relative_fm.h"


int
main(int argc, char** argv)
{
  if(argc < 3)
  {
    std::cerr << "Usage: align_bwts ref seq1 [seq2 ...]" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  std::cout << "Relative FM-index builder" << std::endl;
  std::cout << std::endl;
  std::cout << "Reference: " << argv[1] << std::endl;
  std::cout << std::endl;

  SimpleFM<> ref(argv[1]);
  ref.reportSize(true);
  std::cout << std::endl;

  for(int arg = 2; arg < argc; arg++)
  {
    std::cout << "Target: " << argv[arg] << std::endl;
    SimpleFM<> seq(argv[arg]);
    double start = readTimer();
    RelativeFM rel(ref, seq, true);
    double seconds = readTimer() - start;
    std::cout << "Index built in " << seconds << " seconds" << std::endl;
    std::cout << std::endl;

    rel.writeTo(argv[arg]);
    seq.reportSize(true);
    rel.reportSize(true);
    std::cout << std::endl;
  }

  double memory = inMegabytes(memoryUsage());
  std::cout << "Memory usage: " << memory << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
}
