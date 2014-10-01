#include <cstdlib>
#include <fstream>

#include <sys/resource.h>

#include "utils.h"

//------------------------------------------------------------------------------

void
printSize(const std::string& header, uint64_t bytes, uint64_t data_size, uint indent)
{
  std::string padding;
  if(header.length() + 1 < indent) { padding = std::string(indent - 1 - header.length(), ' '); }

  std::cout << header << ':' << padding << inMegabytes(bytes) << " MB (" << inBPC(bytes, data_size) << " bpc)" << std::endl;
}

void
printTime(const std::string& header, uint64_t found, uint64_t matches, uint64_t bytes, double seconds, uint indent)
{
  std::string padding;
  if(header.length() + 1 < indent) { padding = std::string(indent - 1 - header.length(), ' '); }

  std::cout << header << ':' << padding << "Found " << found << " patterns with " << matches << " occ in "
    << seconds << " seconds (" << (inMegabytes(bytes) / seconds) << " MB / s)" << std::endl;
}

//------------------------------------------------------------------------------

double
readTimer()
{
  return clock() / (double)CLOCKS_PER_SEC;
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

void relativeLZ(const bit_vector& text, const bit_vector& reference, std::vector<range_type>& phrases)
{
  // Handle the reference.
  uint64_t ones = util::cnt_one_bits(reference);
  uint64_t zeros = reference.size() - ones;
  unsigned char* buffer = new unsigned char[reference.size() + 1];
  for(uint64_t i = 0; i < reference.size(); i++)
  {
    buffer[i] = reference[i] + 1;
  }
  buffer[reference.size()] = 0;
  if(ones == 0 || zeros == 0)
  {
    std::cerr << "relativeLZ(): Reference must contain both 0-bits and 1-bits!" << std::endl;
    return;
  }

  // Build SA.
  int_vector<> sa(reference.size(), 0, bits::hi(reference.size()) + 1);
  algorithm::calculate_sa(buffer, reference.size(), sa);
  delete[] buffer; buffer = 0;  

  // Parse the text.
  uint64_t pos = 0;
  phrases.clear();
  while(pos < text.size())
  {
    range_type range = (text[pos] ? range_type(zeros, reference.size() - 1) : range_type(0, zeros - 1));
    uint64_t len = 1;
    while(pos + len < text.size())
    {
      uint64_t low = range.first, high = range.second, last_high = range.second;
      while(low < high) // Lower bound for pattern text[pos, pos + len].
      {
        uint64_t mid = low + (high - low) / 2;
        if(sa[mid] + len >= reference.size() || reference[sa[mid] + len] < text[pos + len]) { low = mid + 1; }
        else if(reference[sa[mid] + len] == text[pos + len]) { high = mid; }
        else { last_high = high = std::max(mid, (uint64_t)1) - 1; }
      }
      if(sa[low] + len >= reference.size() || reference[sa[low] + len] != text[pos + len]) { break; }
      range.first = low; high = last_high;
      while(low < high) // Upper bound for pattern text[pos, pos + len].
      {
        uint64_t mid = low + (high - low + 1) / 2;
        if(reference[sa[mid] + len] == text[pos + len]) { low = mid; }
        else { high = mid - 1; }
      }
      range.second = high; len++;
    }
    phrases.push_back(range_type(sa[range.first], len));
    pos += len;
  }
}

//------------------------------------------------------------------------------
