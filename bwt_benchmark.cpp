#include <fstream>

#include <sdsl/construct.hpp>
#include <sdsl/csa_wt.hpp>

#include "utils.h"


int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: bwt_benchmark input" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  int_vector<8> text;
  std::cout << "File: " << argv[1] << std::endl;
  load_vector_from_file(text, argv[1], 1);
  std::cout << "Text size: " << text.size() << std::endl;

  double start = readTimer();
  csa_wt<bwt_type, 4096, 4096> csa;
  construct_im(csa, text);
  double seconds = readTimer() - start;

  std::cout << "Index built in " << seconds << " seconds" << std::endl;
  double memory = inMegabytes(memoryUsage());
  std::cout << "Memory usage: " << memory << " MB" << std::endl;

  return 0;
}
