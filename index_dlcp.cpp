/*
  Copyright (c) 2015, 2016 Genome Research Ltd.

  Author: Jouni Siren <jouni.siren@iki.fi>

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

#include <sdsl/construct_sa.hpp>

#include "relative_lcp.h"

using namespace relative;

//------------------------------------------------------------------------------

//#define STORE_PLAIN_DLCP

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: index_dlcp input1 [input2 ...]" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  std::cout << "Indexing the differential LCP array" << std::endl;
  std::cout << std::endl;

  for(int i = 1; i < argc; i++)
  {
    double start = readTimer();
    std::string base_name = argv[i];
    std::cout << "File: " << base_name << std::endl;
    std::string lcp_file = base_name + LCP_EXTENSION;
    std::ifstream input(lcp_file.c_str(), std::ios_base::binary);
    if(!input)
    {
      std::cerr << "index_dlcp: Cannot open LCP file " << lcp_file << std::endl;
      std::cout << std::endl;
      continue;
    }
    SLArray lcp; lcp.load(input); input.close();
    std::cout << "LCP size: " << lcp.size() << std::endl;
    double megabytes = lcp.size() / MEGABYTE_DOUBLE;

#ifdef STORE_PLAIN_DLCP
    {
      std::vector<int32_t> dlcp_plain(lcp.size(), 0);
      int32_t prev = 0;
      for(size_type i = 1; i < dlcp_plain.size(); i++)
      {
        int32_t curr = lcp[i];
        dlcp_plain[i] = curr - prev;
        prev = curr;
      }
      std::string dlcp_name = base_name + DLCP_EXTENSION;
      std::ofstream out(dlcp_name.c_str(), std::ios_base::binary);
      if(out)
      {
        out.write((char*)(dlcp_plain.data()), dlcp_plain.size() * sizeof(int32_t));
        out.close();
      }
      else
      {
        std::cerr << "index_dlcp: Cannot open DLCP file " << dlcp_name << std::endl;
      }
    }
#endif

    sdsl::int_vector<0> dlcp; dlcp.width(lcp.large.width() + 2);
    differentialArray<SLArray, DiffEncoderNZ>(lcp, dlcp, true);
    sdsl::util::clear(lcp);
    std::cout << "DLCP width: " << (size_type)(dlcp.width()) << std::endl;

    sdsl::int_vector<0> sa;
    sdsl::qsufsort::construct_sa(sa, dlcp); sdsl::util::bit_compress(sa);
    std::cout << "SA width: " << (size_type)(sa.width()) << std::endl;
    store_to_file(sa, base_name + DLCP_INDEX_EXTENSION);

    double seconds = readTimer() - start;
    std::cout << "Index built in " << seconds << " seconds (" << (megabytes / seconds) << " MB/s)" << std::endl;
    std::cout << std::endl;
  }

  return 0;
}

//------------------------------------------------------------------------------
