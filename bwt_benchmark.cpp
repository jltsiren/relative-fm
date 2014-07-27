#include <fstream>

#include <sdsl/construct.hpp>

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

  std::cout << "File: " << argv[1] << std::endl;
  uint64_t file_size = sdsl::util::file_size(argv[1]); 
  std::cout << "File size: " << file_size << std::endl;
  int_vector<8> text(file_size + 1, 0); // append 0 byte
  uint64_t n = text.size();
  std::ifstream in(argv[1], std::ios_base::binary);
  in.read((char*)text.data(), file_size); 
  std::cout << "Text size: " << text.size() << std::endl;

  // construct SA & BWT
  double start = readTimer();
  int_vector<64> sa(n);
  divsufsort64((const unsigned char*)text.data(), (int64_t*)sa.data(), n);
  uint8_t* bwt = (uint8_t*)sa.data(); // override SA with BWT values
  uint64_t to_add[2] = {(uint64_t)-1,n-1};
  for (uint64_t i=0; i < n; ++i) { bwt[i] = text[sa[i] + to_add[sa[i] == 0]]; }
  double seconds = readTimer() - start;

  std::cout << "BWT built in " << seconds << " seconds" << std::endl;
  double memory = inMegabytes(memoryUsage());
  std::cout << "Memory usage: " << memory << " MB" << std::endl;

  return 0;
}
