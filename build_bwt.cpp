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
    int_vector<8> text;
    uint64_t size = 0;

    // Read text.
    {
      std::ifstream in(base_name.c_str(), std::ios_base::binary);
      if(!in)
      {
        std::cerr << "build_bwt: Cannot open input file " << base_name << std::endl;
        std::cout << std::endl;
        continue;
      }
      size = util::file_size(base_name); text.resize(size + 1);
      std::cout << "Text size: " << size << std::endl;
      in.read((char*)(text.data()), size); text[size] = 0; in.close();
    }

    // Build BWT.
    {
      double start = readTimer();
      int_vector<64> sa(size + 1);
      divsufsort64((const unsigned char*)(text.data()), (int64_t*)(sa.data()), size + 1);
      uint8_t* bwt = (uint8_t*)(sa.data()); // Overwrite SA with BWT.
      uint64_t to_add[2] = { (uint64_t)-1, size };
      for(uint64_t i = 0; i <= size; i++) { bwt[i] = text[sa[i] + to_add[sa[i] == 0]]; }
      for(uint64_t i = 0; i <= size; i++) { text[i] = bwt[i]; }
      util::clear(sa);
      double seconds = readTimer() - start;
      std::cout << "BWT built in " << seconds << " seconds (" << (inMegabytes(size) / seconds) << " MB/s)" << std::endl;
    }

    // Write alphabet.
    {
      std::string ramfile = ram_file_name(base_name + BWT_EXTENSION);
      store_to_file(text, ramfile);
      alphabet_type alpha;
      {
        int_vector_buffer<8> ram_buffer(ramfile);
        alphabet_type temp(ram_buffer, size + 1);
        alpha.swap(temp);
      }
      ram_fs::remove(ramfile);

      std::string filename = base_name + ALPHA_EXTENSION;
      std::ofstream out(filename.c_str(), std::ios_base::binary);
      if(!out)
      {
        std::cerr << "build_bwt: Cannot open alphabet file " << filename << std::endl;
        std::cout << std::endl;
        continue;
      }
      alpha.serialize(out); out.close();
      std::cout << "Alphabet written to " << filename << std::endl;
      for(uint64_t i = 0; i <= size; i++) { text[i] = alpha.char2comp[text[i]]; }
    }

    // Write BWT.
    {
      std::string filename = base_name + BWT_EXTENSION;
      std::ofstream out(filename.c_str(), std::ios_base::binary);
      if(!out)
      {
        std::cerr << "build_bwt: Cannot open BWT file " << filename << std::endl;
        std::cout << std::endl;
        continue;
      }
      text.serialize(out); out.close();
      std::cout << "BWT written to " << filename << std::endl;
    }

    std::cout << std::endl;
  }

  return 0;
}
