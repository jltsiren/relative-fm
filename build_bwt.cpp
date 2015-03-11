#include <cstdlib>
#include <fstream>
#include <unistd.h>

#include <sdsl/construct.hpp>

#include "utils.h"

using namespace relative;


int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: build_bwt [options] input1 [input2 ...]" << std::endl;
    std::cerr << "  -a    Write the alphabet file." << std::endl;
    std::cerr << "  -s N  Sample one out of N SA values (default 0)." << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  bool write_alphabet = false, options = false;
  uint64_t sample_rate = 0;
  int c = 0;
  while((c = getopt(argc, argv, "as:")) != -1)
  {
    switch(c)
    {
    case 'a':
      write_alphabet = true; options = true; break;
    case 's':
      sample_rate = atol(optarg); options = true; break;
    case '?':
      return 2;
    default:
      return 3;
    }
  }

  std::cout << "BWT construction" << std::endl;
  if(options)
  {
    std::cout << "Options:";
    if(write_alphabet) { std::cout << " alphabet"; }
    if(sample_rate > 0) { std::cout << " sample_rate=" << sample_rate; }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  for(int i = optind; i < argc; i++)
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

    // Build BWT and sample SA.
    int_vector<0> samples;
    {
      double start = readTimer();
      int_vector<64> sa(size + 1);
      divsufsort64((const unsigned char*)(text.data()), (int64_t*)(sa.data()), size + 1);
      if(sample_rate > 0)
      {
        util::assign(samples, int_vector<0>(size / sample_rate + 1, 0, bitlength(size)));
        for(uint64_t i = 0; i <= size; i += sample_rate) { samples[i / sample_rate] = sa[i]; }
      }
      uint8_t* bwt = (uint8_t*)(sa.data()); // Overwrite SA with BWT.
      uint64_t to_add[2] = { (uint64_t)-1, size };
      for(uint64_t i = 0; i <= size; i++) { bwt[i] = text[sa[i] + to_add[sa[i] == 0]]; }
      for(uint64_t i = 0; i <= size; i++) { text[i] = bwt[i]; }
      util::clear(sa);
      double seconds = readTimer() - start;
      std::cout << "BWT built in " << seconds << " seconds (" << (inMegabytes(size) / seconds) << " MB/s)" << std::endl;
    }

    // Compact the alphabet and write it if necessary.
    {
      Alphabet alpha(text);
      for(uint64_t i = 0; i <= size; i++) { text[i] = alpha.char2comp[text[i]]; }
      if(write_alphabet)
      {
        std::string filename = base_name + ALPHA_EXTENSION;
        std::ofstream out(filename.c_str(), std::ios_base::binary);
        if(!out)
        {
          std::cerr << "build_bwt: Cannot open alphabet file " << filename << std::endl;
        }
        else
        {
          alpha.serialize(out); out.close();
          std::cout << "Alphabet written to " << filename << std::endl;
        }
      }
    }

    // Write BWT.
    {
      std::string filename = base_name + BWT_EXTENSION;
      std::ofstream out(filename.c_str(), std::ios_base::binary);
      if(!out)
      {
        std::cerr << "build_bwt: Cannot open BWT file " << filename << std::endl;
      }
      else
      {
        text.serialize(out); out.close();
        std::cout << "BWT written to " << filename << std::endl;
      }
    }

    // Write SA samples.
    if(sample_rate > 0)
    {
      std::string filename = base_name + SAMPLE_EXTENSION;
      std::ofstream out(filename.c_str(), std::ios_base::binary);
      if(!out)
      {
        std::cerr << "build_bwt: Cannot open SA sample file " << filename << std::endl;
      }
      else
      {
        write_member(sample_rate, out); samples.serialize(out); out.close();
        std::cout << "SA samples written to " << filename << std::endl;
      }
    }

    std::cout << std::endl;
  }

  return 0;
}
