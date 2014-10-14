#include "rlz_fm.h"


int
main(int argc, char** argv)
{
  if(argc < 3)
  {
    std::cerr << "Usage: build_rlzfm ref seq1 [seq2 ...]" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  std::cout << "Relative Lempel-Ziv FM-index builder" << std::endl;
  std::cout << std::endl;
  std::cout << "Reference: " << argv[1] << std::endl;
  std::cout << std::endl;

  SimpleFM<> ref(argv[1]);
  ref.reportSize(true);
  std::cout << std::endl;
  csa_wt<> csa;
  {
    std::string bwt_name = std::string(argv[1]) + BWT_EXTENSION;
    std::ifstream in(bwt_name.c_str(), std::ios_base::binary);
    if(!in)
    {
      std::cerr << "build_rlzfm: Cannot open input file " << bwt_name.c_str() << std::endl;
      return 2;
    }
    int_vector<8> buffer; buffer.load(in);
    in.close();
    reverseIndex(buffer, csa);
  }

  for(int arg = 2; arg < argc; arg++)
  {
    std::cout << "Target: " << argv[arg] << std::endl;
    SimpleFM<> seq(argv[arg]);
    double start = readTimer();
    RLZFM rel(ref, seq, &csa);
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
