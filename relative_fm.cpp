#include "relative_fm.h"

namespace relative
{

//------------------------------------------------------------------------------

uint64_t
mostFrequentChar(std::vector<uint8_t>& ref_buffer, std::vector<uint8_t>& seq_buffer,
  bit_vector& ref_lcs, bit_vector& seq_lcs,
  range_type ref_range, range_type seq_range)
{
  std::vector<range_type> freqs(256, range_type(0, 0));
  for(auto c : ref_buffer) { freqs[c].first++; }
  for(auto c : seq_buffer) { freqs[c].second++; }
  for(auto& f : freqs) { f.first = std::min(f.first, f.second); }
  auto max_pos = std::max_element(freqs.begin(), freqs.end());
  uint64_t c = max_pos - freqs.begin(), max_f = max_pos->first;

  for(uint64_t i = 0, found = 0; i < ref_buffer.size() && found < max_f; i++)
  {
    if(ref_buffer[i] == c) { ref_lcs[ref_range.first + i] = 1; found++; }
  }
  for(uint64_t i = 0, found = 0; i < seq_buffer.size() && found < max_f; i++)
  {
    if(seq_buffer[i] == c) { seq_lcs[seq_range.first + i] = 1; found++; }
  }

#ifdef VERBOSE_STATUS_INFO
  std::cout << "  Matched " << max_f << " copies of character " << c << " (" << ((char)c) << ")";
  std::cout << " from sequences of length " << range_type(ref_buffer.size(), seq_buffer.size()) << std::endl;
#endif
  return max_f;
}

//------------------------------------------------------------------------------

} // namespace relative
