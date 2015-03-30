#include <sdsl/construct_sa.hpp>

#include "relative_lcp.h"

using namespace relative;

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

    int_vector<0> dlcp; RelativeLCP::differentialArray(lcp, dlcp, true);
    util::clear(lcp);
    std::cout << "DLCP width: " << (uint64_t)(dlcp.width()) << std::endl;

    int_vector<0> sa;
    qsufsort::construct_sa(sa, dlcp); util::bit_compress(sa);
    std::cout << "SA width: " << (uint64_t)(sa.width()) << std::endl;
    store_to_file(sa, base_name + DLCP_INDEX_EXTENSION);

    double seconds = readTimer() - start;
    std::cout << "Index built in " << seconds << " seconds (" << (megabytes / seconds) << " MB/s)" << std::endl;
    std::cout << std::endl;
  }

  return 0;
}

//------------------------------------------------------------------------------
