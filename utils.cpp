/*
  Copyright (c) 2015, 2016, 2017 Genome Research Ltd.
  Copyright (c) 2014 Jouni Siren

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
const std::string DLCP_EXTENSION = ".dlcp";
const std::string DLCP_INDEX_EXTENSION = ".dlcp_index";
const std::string SIMPLE_FM_DEFAULT_ALPHABET("\0ACGNT", 6);

//------------------------------------------------------------------------------

void
printHeader(const std::string& header, size_type indent)
{
  std::string padding;
  if(header.length() + 1 < indent) { padding = std::string(indent - 1 - header.length(), ' '); }
  std::cout << header << ":" << padding;
}

void
printSize(const std::string& header, size_type bytes, size_type data_size, size_type indent)
{
  printHeader(header, indent);
  std::cout << inMegabytes(bytes) << " MB (" << inBPC(bytes, data_size) << " bpc)" << std::endl;
}

void
printTime(const std::string& header, size_type found, size_type matches, size_type bytes, double seconds, bool occs, size_type indent)
{
  printHeader(header, indent);

  std::cout << "Found " << found << " patterns with " << matches << " occ in " << seconds << " seconds (";
  if(occs) { std::cout << (matches / seconds) << " occ/s)" << std::endl; }
  else     { std::cout << (inMegabytes(bytes) / seconds) << " MB/s)" << std::endl; }
}

void
printTime(const std::string& header, size_type queries, double seconds, size_type indent)
{
  printHeader(header, indent);
  std::cout << queries << " queries in " << seconds << " seconds ("
            << inMicroseconds(seconds / queries) << " µs/query)" << std::endl;
}

//------------------------------------------------------------------------------

double
readTimer()
{
  return omp_get_wtime();
}

size_type
memoryUsage()
{
  rusage usage;
  getrusage(RUSAGE_SELF, &usage);
#ifdef RUSAGE_IN_BYTES
  return usage.ru_maxrss;
#else
  return KILOBYTE * usage.ru_maxrss;
#endif
}

//------------------------------------------------------------------------------

size_type
readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows)
{
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "readRows(): Cannot open input file " << filename << std::endl;
    return 0;
  }

  size_type chars = 0;
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
