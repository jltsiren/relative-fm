/*
  Copyright (c) 2015, 2016, 2017 Genome Research Ltd.
  Copyright (c) 2014 Jouni Siren and Simon Gog

  Author: Jouni Siren <jouni.siren@iki.fi>
  Author: Simon Gog

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#include <cstdlib>
#include <unistd.h>

#include <sdsl/construct.hpp>
#include <sdsl/lcp.hpp>

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
    std::cerr << "  -l    Build the LCP array." << std::endl;
    std::cerr << "  -s N  Sample one out of N SA values (default 0)." << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  bool write_alphabet = false, build_lcp = false, options = false;
  size_type sa_sample_rate = 0, isa_sample_rate = 0;
  int c = 0;
  while((c = getopt(argc, argv, "ai:ls:")) != -1)
  {
    switch(c)
    {
    case 'a':
      write_alphabet = true; options = true; break;
    case 'i':
      isa_sample_rate = atol(optarg); options = true; break;
    case 'l':
      build_lcp = true; options = true; break;
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
    if(isa_sample_rate > 0) { std::cout << " isa_sample_rate=" << isa_sample_rate; }
    if(build_lcp) { std::cout << " lcp"; }
    if(sa_sample_rate > 0) { std::cout << " sa_sample_rate=" << sa_sample_rate; }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  for(int i = optind; i < argc; i++)
  {
    std::string base_name = argv[i];
    std::cout << "File: " << base_name << std::endl;
    sdsl::int_vector<8> text;
    sdsl::cache_config config;
    size_type size = 0;

    // Read text.
    {
      std::ifstream in(base_name.c_str(), std::ios_base::binary);
      if(!in)
      {
        std::cerr << "build_bwt: Cannot open input file " << base_name << std::endl;
        std::cout << std::endl;
        continue;
      }
      size = sdsl::util::file_size(base_name); text.resize(size + 1);
      std::cout << "Text size: " << size << std::endl;
      in.read((char*)(text.data()), size); text[size] = 0; in.close();
    }

    // Build BWT and sample SA.
    sdsl::int_vector<0> sa_samples, isa_samples;
    {
      double start = readTimer();
      sdsl::int_vector<64> sa(size + 1);
      divsufsort64((const unsigned char*)(text.data()), (int64_t*)(sa.data()), size + 1);
      if(sa_sample_rate > 0)
      {
        sdsl::util::assign(sa_samples, sdsl::int_vector<0>(size / sa_sample_rate + 1, 0, bit_length(size)));
        for(size_type i = 0; i <= size; i += sa_sample_rate) { sa_samples[i / sa_sample_rate] = sa[i]; }
      }
      if(isa_sample_rate > 0)
      {
        sdsl::util::assign(isa_samples, sdsl::int_vector<0>(size / isa_sample_rate + 1, 0, bit_length(size)));
        for(size_type i = 0; i <= size; i++)
        {
          if(sa[i] % isa_sample_rate == 0) { isa_samples[sa[i] / isa_sample_rate] = i; }
        }
      }
      char_type* bwt = (char_type*)(sa.data()); // Overwrite SA with BWT.
      size_type to_add[2] = { (size_type)-1, size };
      for(size_type i = 0; i <= size; i++) { bwt[i] = text[sa[i] + to_add[sa[i] == 0]]; }
      for(size_type i = 0; i <= size; i++) { text[i] = bwt[i]; }
      sdsl::util::clear(sa);
      double seconds = readTimer() - start;
      std::cout << "BWT built in " << seconds << " seconds (" << (inMegabytes(size) / seconds) << " MB/s)" << std::endl;
    }

    // Compact the alphabet and write it if necessary.
    {
      Alphabet alpha(text);
      for(size_type i = 0; i <= size; i++) { text[i] = alpha.char2comp[text[i]]; }
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
        if(build_lcp)
        {
          config.file_map[sdsl::conf::KEY_BWT] = filename;
        }
      }
      sdsl::util::clear(text);
    }

    // Write SA/ISA samples.
    if(sa_sample_rate > 0 || isa_sample_rate > 0)
    {
      std::string filename = base_name + SAMPLE_EXTENSION;
      std::ofstream out(filename.c_str(), std::ios_base::binary);
      if(!out)
      {
        std::cerr << "build_bwt: Cannot open sample file " << filename << std::endl;
      }
      else
      {
        sdsl::write_member(sa_sample_rate, out); sa_samples.serialize(out);
        sdsl::write_member(isa_sample_rate, out); isa_samples.serialize(out);
        out.close();
        std::cout << "Samples written to " << filename << std::endl;
      }
      sdsl::util::clear(sa_samples); sdsl::util::clear(isa_samples);
    }

    // Build and write LCP.
    if(build_lcp)
    {
      double start = readTimer();
      construct_lcp_bwt_based(config);
      sdsl::int_vector_buffer<0> lcp_buffer(sdsl::cache_file_name(sdsl::conf::KEY_LCP, config));
      SLArray lcp(lcp_buffer);
      double seconds = readTimer() - start;
      std::cout << "LCP array built in " << seconds << " seconds (" << (inMegabytes(size) / seconds) << " MB/s)" << std::endl;
      std::string filename = base_name + LCP_EXTENSION;
      sdsl::store_to_file(lcp, filename);
    }

    std::cout << std::endl;
  }

  return 0;
}
