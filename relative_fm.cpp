#include "relative_fm.h"

//------------------------------------------------------------------------------

void
verifyRanges(std::vector<record_type>& ranges, uint64_t ref_len, uint64_t seq_len)
{
  range_type expect(0, 0);
  for(auto curr : ranges)
  {
    if(curr.left.first != expect.first ||  curr.right.first != expect.second)
    {
      std::cerr << "verifyRanges(): Expected " << expect << ", got " << curr << std::endl;
    }
    expect = range_type(curr.left.second + 1, curr.right.second + 1);
  }
  if(expect.first != ref_len || expect.second != seq_len)
  {
    std::cerr << "verifyRanges(): Final expect was " << expect << std::endl;
  }
}

std::vector<std::pair<int, int> >
mostFrequentChar(std::vector<uint8_t>& ref_buffer, std::vector<uint8_t>& seq_buffer)
{
  std::vector<std::pair<int, int> > freqs(256, std::make_pair(0, 0));
  for(auto c : ref_buffer) { freqs[c].first++; }
  for(auto c : seq_buffer) { freqs[c].second++; }
  for(auto& f : freqs) { f.first = std::min(f.first, f.second); }
  auto max_pos = std::max_element(freqs.begin(), freqs.end());
  uint max_c = max_pos - freqs.begin(), max_f = max_pos->first;

  std::vector<std::pair<int, int> > res; res.reserve(max_f);
  auto ref_p = ref_buffer.begin(), seq_p = seq_buffer.begin();
  while(res.size() < max_f)
  {
    while(*ref_p != max_c) { ++ref_p; }
    while(*seq_p != max_c) { ++seq_p; }
    res.push_back(std::make_pair(ref_p - ref_buffer.begin(), seq_p - seq_buffer.begin()));
    ++ref_p; ++seq_p;
  }

#ifdef VERBOSE_STATUS_INFO
  std::cout << "  Matched " << max_f << " copies of character " << max_c << " (" << ((char)max_c) << ")";
  std::cout << " from sequences of length " << range_type(ref_buffer.size(), seq_buffer.size()) << std::endl;
#endif

  return res;
}

//------------------------------------------------------------------------------
