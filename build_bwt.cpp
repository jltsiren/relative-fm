#include <cstdlib>
#include <unistd.h>

#include <sdsl/construct.hpp>

#include "support.h"

using namespace relative;


int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: build_bwt [options] input1 [input2 ...]" << std::endl;
    std::cerr << "  -a    Write the alphabet file." << std::endl;
    std::cerr << "  -i N  Sample one out of N ISA values (default 0)." << std::endl;
    std::cerr << "  -s N  Sample one out of N SA values (default 0)." << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  bool write_alphabet = false, options = false;
  uint64_t sa_sample_rate = 0, isa_sample_rate = 0;
  int c = 0;
  while((c = getopt(argc, argv, "ai:s:")) != -1)
  {
    switch(c)
    {
    case 'a':
      write_alphabet = true; options = true; break;
    case 'i':
      isa_sample_rate = atol(optarg); options = true; break;
    case 's':
      sa_sample_rate = atol(optarg); options = true; break;
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
    if(sa_sample_rate > 0) { std::cout << " sa_sample_rate=" << sa_sample_rate; }
    if(isa_sample_rate > 0) { std::cout << " isa_sample_rate=" << isa_sample_rate; }
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
    int_vector<0> sa_samples, isa_samples;
    {
      double start = readTimer();
      int_vector<64> sa(size + 1);
      divsufsort64((const unsigned char*)(text.data()), (int64_t*)(sa.data()), size + 1);
      if(sa_sample_rate > 0)
      {
        util::assign(sa_samples, int_vector<0>(size / sa_sample_rate + 1, 0, bitlength(size)));
        for(uint64_t i = 0; i <= size; i += sa_sample_rate) { sa_samples[i / sa_sample_rate] = sa[i]; }
      }
      if(isa_sample_rate > 0)
      {
        util::assign(isa_samples, int_vector<0>(size / isa_sample_rate + 1, 0, bitlength(size)));
        for(uint64_t i = 0; i <= size; i++)
        {
          if(sa[i] % isa_sample_rate == 0) { isa_samples[sa[i] / isa_sample_rate] = i; }
        }
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

    // Write SA/ISA samples.
    if(sa_sample_rate > 0 || isa_sample_rate > 0)
    {
      std::string filename = base_name + SAMPLE_EXTENSION;
      std::ofstream out(filename.c_str(), std::ios_base::binary);
      if(!out)
      {
        std::cerr << "build_bwt: Cannot open the sample file " << filename << std::endl;
      }
      else
      {
        write_member(sa_sample_rate, out); sa_samples.serialize(out);
        write_member(isa_sample_rate, out); isa_samples.serialize(out);
        out.close();
        std::cout << "Samples written to " << filename << std::endl;
      }
    }

    std::cout << std::endl;
  }

  return 0;
}
