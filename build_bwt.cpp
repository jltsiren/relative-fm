#include <fstream>

#include <sdsl/construct.hpp>

#include "utils.h"


int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: build_bwt input1 [input2 ...]" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  for(int i = 1; i < argc; i++)
  {
    std::string base_name = argv[i];
    std::cout << "File: " << base_name << std::endl;
    uint64_t size = util::file_size(base_name);
    std::cout << "Text size: " << size << std::endl;

    // Read text.
    int_vector<8> text(size + 1, 0);  // Append an endmarker.
    std::ifstream in(base_name.c_str(), std::ios_base::binary);
    if(!in)
    {
      std::cerr << "build_bwt: Cannot open input file " << base_name << std::endl;
      std::cout << std::endl;
      continue;
    }
    in.read((char*)(text.data()), size);

    // Build BWT.
    double start = readTimer();
    divbwt64((const unsigned char*)(text.data()), (unsigned char*)(text.data()), 0, size + 1);
    double seconds = readTimer() - start;
    std::cout << "BWT built in " << seconds << " seconds (" << (inMegabytes(size) / seconds) << " MB/s)" << std::endl;

    // Write BWT.
    std::string bwt_name = base_name + BWT_EXTENSION;
    std::ofstream out(bwt_name.c_str(), std::ios_base::binary);
    if(!out)
    {
      std::cerr << "build_bwt: Cannot open output file " << bwt_name << std::endl;
      std::cout << std::endl;
      continue;
    }
    text.serialize(out);
    out.close();
    std::cout << "BWT written to " << bwt_name << std::endl;
    
    std::cout << std::endl;
  }

  return 0;
}
