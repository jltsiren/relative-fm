#include <cstdlib>
#include <fstream>

#include <sys/resource.h>

#include "utils.h"

namespace relative
{

//------------------------------------------------------------------------------

const std::string BWT_EXTENSION = ".bwt";
const std::string NATIVE_BWT_EXTENSION = ".cbwt";
const std::string ALPHA_EXTENSION = ".alpha";
const std::string SAMPLE_EXTENSION = ".samples";
const std::string LCP_EXTENSION = ".lcp";
const std::string SIMPLE_FM_DEFAULT_ALPHABET("\0ACGNT", 6);
const std::string ROPEBWT_ALPHABET("\0ACGTN", 6);

//------------------------------------------------------------------------------

void
printSize(const std::string& header, uint64_t bytes, uint64_t data_size, uint64_t indent)
{
  std::string padding;
  if(header.length() + 1 < indent) { padding = std::string(indent - 1 - header.length(), ' '); }

  std::cout << header << ':' << padding << inMegabytes(bytes) << " MB (" << inBPC(bytes, data_size) << " bpc)" << std::endl;
}

void
printTime(const std::string& header, uint64_t found, uint64_t matches, uint64_t bytes, double seconds, bool occs, uint64_t indent)
{
  std::string padding;
  if(header.length() + 1 < indent) { padding = std::string(indent - 1 - header.length(), ' '); }

  std::cout << header << ':' << padding << "Found " << found << " patterns with " << matches << " occ in "
    << seconds << " seconds (";
  if(occs) { std::cout << (matches / seconds) << " occ/s)" << std::endl; }
  else     { std::cout << (inMegabytes(bytes) / seconds) << " MB/s)" << std::endl; }
}

//------------------------------------------------------------------------------

double
readTimer()
{
  return omp_get_wtime();
}

uint64_t
memoryUsage()
{
  rusage usage;
  getrusage(RUSAGE_SELF, &usage);
#ifdef RUSAGE_IN_BYTES
  return usage.ru_maxrss;
#else
  return ((uint64_t)1024) * usage.ru_maxrss;
#endif
}

//------------------------------------------------------------------------------

uint64_t
readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows)
{
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "readRows(): Cannot open input file " << filename << std::endl;
    return 0;
  }

  uint64_t chars = 0;
  while(input)
  {
    std::string buf;
    std::getline(input, buf);
    if(skip_empty_rows && buf.length() == 0) { continue; }
    rows.push_back(buf);
    chars += buf.length();
  }

  input.close();
  return chars;
}

//------------------------------------------------------------------------------

} // namespace relative
