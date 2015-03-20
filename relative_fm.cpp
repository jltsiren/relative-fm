#include <map>

#include "relative_fm.h"

namespace relative
{

//------------------------------------------------------------------------------

range_type
mostFrequentChar(const std::vector<uint8_t>& ref_extract, const std::vector<uint8_t>& seq_extract,
  bit_vector& ref_lcs, bit_vector& seq_lcs,
  range_type ref_range, range_type seq_range)
{
  std::vector<range_type> freqs(256, range_type(0, 0));
  for(uint64_t i = 0; i < ref_extract.size(); i++) { freqs[ref_extract[i]].first++; }
  for(uint64_t i = 0; i < seq_extract.size(); i++) { freqs[seq_extract[i]].second++; }
  uint64_t max_c = 0, max_f = 0;
  for(uint64_t c = 0; c < freqs.size(); c++)
  {
    uint64_t temp = std::min(freqs[c].first, freqs[c].second);
    if(temp > max_f) { max_c = c; max_f = temp; }
  }

  uint64_t ref_padding = ref_range.first % 64, seq_padding = seq_range.first % 64;
  for(uint64_t i = 0, found = 0; i < length(ref_range) && found < max_f; i++)
  {
    if(ref_extract[i] == max_c) { ref_lcs[ref_padding + i] = 1; found++; }
  }
  for(uint64_t i = 0, found = 0; i < length(seq_range) && found < max_f; i++)
  {
    if(seq_extract[i] == max_c) { seq_lcs[seq_padding + i] = 1; found++; }
  }

#ifdef VERBOSE_STATUS_INFO
  #pragma omp critical (stderr)
  {
    std::cerr << "Finding the most frequent character in " << std::make_pair(ref_range, seq_range) << std::endl;
    if(max_f == 0)
    {
      for(uint64_t c = 0; c < freqs.size(); c++)
      {
        if(freqs[c].first != 0 || freqs[c].second != 0)
        {
          std::cerr << "  Character " << c << " (" << ((char)c) << "): " << freqs[c] << std::endl;
        }
      }
    }
    std::cerr << "  Matched " << max_f << " copies of character " << max_c << " (" << ((char)max_c) << ")";
    std::cerr << " from sequences of length " << range_type(ref_extract.size(), seq_extract.size()) << std::endl;
  }
#endif
  return range_type(max_f, std::min(length(ref_range), length(seq_range)) - max_f);
}

range_type
matchRuns(const std::vector<uint8_t>& ref_extract, const std::vector<uint8_t>& seq_extract,
  bit_vector& ref_lcs, bit_vector& seq_lcs,
  range_type ref_range, range_type seq_range,
  const uint8_t* alphabet, uint64_t sigma)
{
  uint64_t ref_pos = 0, seq_pos = 0, found = 0;
  range_type ref_run(0, 0), seq_run(0, 0);
  uint64_t ref_padding = ref_range.first % 64, seq_padding = seq_range.first % 64;
  while(ref_pos < length(ref_range) && seq_pos < length(seq_range))
  {
    if(ref_run.second <= ref_pos) // Find the next ref_run if needed.
    {
      ref_run.first = ref_extract[ref_pos]; ref_run.second = ref_pos + 1;
      while(ref_run.second < length(ref_range) && ref_extract[ref_run.second] == ref_run.first) { ref_run.second++; }
    }
    if(seq_run.second <= seq_pos) // Find the next seq_run if needed.
    {
      seq_run.first = seq_extract[seq_pos]; seq_run.second = seq_pos + 1;
      while(seq_run.second < length(seq_range) && seq_extract[seq_run.second] == seq_run.first) { seq_run.second++; }
    }
    if(ref_run.first == seq_run.first)  // Match the current runs if the characters match.
    {
      uint64_t run_length = std::min(ref_run.second - ref_pos, seq_run.second - seq_pos);
      for(uint64_t i = 0; i < run_length; i++)
      {
        ref_lcs[ref_padding + ref_pos + i] = 1; seq_lcs[seq_padding + seq_pos + i] = 1;
      }
      ref_pos = ref_run.second; seq_pos = seq_run.second; found += run_length;
    }
    else  // Otherwise advance the run with the smaller character comp value.
    {
      uint64_t ref_c = 0, seq_c = 0;
      for(uint64_t c = 0; c < sigma; c++)
      {
        if(alphabet[c] == ref_run.first) { ref_c = c; }
        if(alphabet[c] == seq_run.first) { seq_c = c; }
      }
      if(ref_c < seq_c) { ref_pos = ref_run.second; }
      else              { seq_pos = seq_run.second; }
    }
  }

#ifdef VERBOSE_STATUS_INFO
  #pragma omp critical (stderr)
  {
    std::cerr << "Matching the runs in " << std::make_pair(ref_range, seq_range) << std::endl;
    std::cerr << "  Matched " << found << " characters from sequences of length "
              << range_type(ref_extract.size(), seq_extract.size()) << std::endl;
  }
#endif
  return range_type(found, std::min(length(ref_range), length(seq_range)) - found);
}

// Diagonal i has 2i+3 elements: -(i+1) to (i+1).
inline int offsetFor(int pos, int d) { return (d + 3) * d + pos + 1; }

range_type
greedyLCS(const std::vector<uint8_t>& ref_extract, const std::vector<uint8_t>& seq_extract,
  bit_vector& ref_lcs, bit_vector& seq_lcs,
  short_record_type& record,
  const uint8_t* alphabet, uint64_t sigma,
  const align_parameters& parameters,
  std::vector<int>* buf)
{
  range_type ref_range = record.first(), seq_range = record.second();
  bool onlyNs = record.onlyNs(), endmarker = record.endmarker();
  if(isEmpty(ref_range) || isEmpty(seq_range)) { return range_type(0, 0); }
  int ref_len = length(ref_range), seq_len = length(seq_range);

  if(endmarker)
  {
    return matchRuns(ref_extract, seq_extract, ref_lcs, seq_lcs, ref_range, seq_range, alphabet, sigma);
  }
  else if(onlyNs || abs(ref_len - seq_len) > parameters.max_d)
  {
    return mostFrequentChar(ref_extract, seq_extract, ref_lcs, seq_lcs, ref_range, seq_range);
  }

  // buffer[offsetFor(k, d)] stores how many characters of ref have been processed on diagonal k.
  bool found = false;
  if(!(parameters.preallocate)) { buf = new std::vector<int>; buf->reserve(align_parameters::BUFFER_SIZE); }
  std::vector<int>& buffer = *buf; buffer.resize(3, 0);
  int d = 0;  // Current diagonal; becomes the number of diagonals.
  for(d = 0; !found && d <= parameters.max_d; d++)
  {
    for(int k = -d; k <= d; k += 2)
    {
      int x = 0;
      if(k == -d || (k != d && buffer[offsetFor(k - 1, d)] < buffer[offsetFor(k + 1, d)]))
      {
        x = buffer[offsetFor(k + 1, d)];
      }
      else
      {
        x = buffer[offsetFor(k - 1, d)] + 1;
      }
      int y = x - k;
      while(x < ref_len && y < seq_len && ref_extract[x] == seq_extract[y]) { x++; y++; }
      buffer[offsetFor(k, d)] = x;
      if(x >= ref_len && y >= seq_len) { found = true; }
    }
    buffer.resize(buffer.size() + 2 * d + 5, 0);  // Add space for diagonal d+1.
    for(int i = 0, curr_pos = offsetFor(-d, d), next_pos = offsetFor(-d, d + 1); i < 2 * d + 3; i++)
    {
      buffer[next_pos + i] = buffer[curr_pos + i];  // Initialize the next diagonal.
    }

    // Too much memory required, match just the most frequent characters.
    if(d >= parameters.max_d && !found)
    {
#ifdef VERBOSE_STATUS_INFO
      #pragma omp critical (stderr)
      std::cerr << "MAX_D exceeded on ranges " << std::make_pair(ref_range, seq_range) << std::endl;
#endif
      if(!(parameters.preallocate)) { delete buf; buf = 0; }
      return mostFrequentChar(ref_extract, seq_extract, ref_lcs, seq_lcs, ref_range, seq_range);
    }
  }

  // Extract the LCS.
  uint64_t lcs = 0;
  d--;  // The last diagonal.
  uint64_t ref_padding = ref_range.first % 64, seq_padding = seq_range.first % 64;
  for(int k = ref_len - seq_len; d >= 0; d--)
  {
    int x_lim = 0, next_k = 0;
    if(d == 0) { x_lim = 0; }
    else if(k == -d || (k != d && buffer[offsetFor(k - 1, d - 1)] < buffer[offsetFor(k + 1, d - 1)]))
    {
      x_lim = buffer[offsetFor(k + 1, d - 1)]; next_k = k + 1;
    }
    else { x_lim = buffer[offsetFor(k - 1, d - 1)] + 1; next_k = k - 1; }
    for(int x = buffer[offsetFor(k, d)]; x > x_lim; x--)
    {
      ref_lcs[ref_padding + x - 1] = 1;
      seq_lcs[seq_padding + x - 1 - k] = 1;
      lcs++;
    }
    k = next_k;
  }

  if(!(parameters.preallocate)) { delete buf; buf = 0; }
  return range_type(lcs, 0);
}

//------------------------------------------------------------------------------

void
copyBits(const bit_vector& source, bit_vector& target, range_type range)
{
  // Copy the first bits.
  uint64_t offset = source.size() - length(range);
  uint64_t bits = std::min((uint64_t)64, source.size()) - offset;
  uint64_t val = source.get_int(offset, bits);
  target.set_int(range.first, val, bits);

  // Copy the bulk of the bits.
  for(uint64_t source_pos = 1, target_pos = (range.first + bits) / 64;
      source_pos < source.size() / 64; source_pos++, target_pos++)
  {
    target.data()[target_pos] = source.data()[source_pos];
  }

  // Copy the last bits.
  if(source.size() > 64 && (bits = source.size() % 64) > 0)
  {
    val = source.get_int(source.size() - bits, bits);
    target.set_int(range.second + 1 - bits, val, bits);
  }
}

//------------------------------------------------------------------------------

struct range_pair_comparator
{
  bool operator() (const range_pair& a, const range_pair& b) const
  {
    return (a.ref_range < b.ref_range);
  }
};

void
processRange(uint64_t pos, uint64_t val, uint64_t len,
  int_vector<0>& smallest, int_vector<0>& previous, uint64_t& lcs)
{
  // Find the longest increasing subsequence with last value < val.
  uint64_t low = 1, high = lcs;
  while(low <= high)
  {
    uint64_t mid = low + (high - low) / 2;
    if(smallest[2 * mid + 1] < val) { low = mid + 1; }
    else { high = mid - 1; }
  }

  // low is now the length of the found subsequence + 1.
  for(uint64_t i = 0; i < len; i++)
  {
    previous[2 * (pos + i)] = smallest[2 * (low + i) - 2];
    previous[2 * (pos + i) + 1] = smallest[2 * (low + i) - 1];
    smallest[2 * (low + i)] = pos + i;
    smallest[2 * (low + i) + 1] = val + i;
  }
  lcs = std::max(lcs, low + len - 1);
}

void
sortMatches(std::vector<range_pair>& matches, uint64_t ref_len, uint64_t seq_len)
{
  range_pair_comparator comp;
  parallelMergeSort(matches.begin(), matches.end(), comp);
  matches.push_back(range_pair(range_type(ref_len, ref_len + 1), range_type(seq_len, seq_len + 1)));
}

uint64_t
countCoverage(std::vector<range_pair>& left_matches, std::vector<range_pair>& right_matches, uint64_t ref_len)
{
  uint64_t total_matches = 0, pos = 0;
  std::vector<range_pair>::iterator l = left_matches.begin(), r = right_matches.begin();
  while(l->ref_range.first < ref_len || r->ref_range.first < ref_len)
  {
    if(l->ref_range.first <= r->ref_range.first)
    {
      if(l->ref_range.second > pos)
      {
        total_matches += l->ref_range.second - std::max(l->ref_range.first, pos);
        pos = l->ref_range.second;
      }
      ++l;
    }
    else
    {
      if(r->ref_range.second > pos)
      {
        total_matches += r->ref_range.second - std::max(r->ref_range.first, pos);
        pos = r->ref_range.second;
      }
      ++r;
    }
  }
  return total_matches;
}

uint64_t
increasingSubsequence(std::vector<range_pair>& left_matches, std::vector<range_pair>& right_matches,
  bit_vector& ref_lcs, bit_vector& seq_lcs, uint64_t ref_len, uint64_t seq_len)
{
  if(left_matches.size() == 0 || right_matches.size() == 0)
  {
    util::assign(ref_lcs, bit_vector(ref_len)); util::assign(seq_lcs, bit_vector(seq_len));
    return 0;
  }

  /*
    Let smallest(i) = smallest[2 * i, 2 * i + 1] and previous(i) = previous[2 * i, 2 * i + 1].
    smallest(i) = (pos, val) is the ending position of the increasing subsequence of
    length i with the smallest val. previous(pos) = (pos', val') is the predecessor of
    position pos in the longest increasing subsequence pos is a part of. Note that if pos
    is a right match, the actual position is pos - ref_len.
  */
  int_vector<0> smallest(2 * (ref_len + 1), 2 * ref_len, bitlength(2 * ref_len));
  int_vector<0> previous(4 * ref_len, 2 * ref_len, bitlength(2 * ref_len));
  uint64_t lcs = 0;
  {
    std::vector<range_pair>::iterator l = left_matches.begin(), r = right_matches.begin();
    while(l->ref_range.first < ref_len || r->ref_range.first < ref_len)
    {
      if(l->ref_range.first < r->ref_range.first)
      {
        uint64_t len = std::min(l->ref_range.second, r->ref_range.first) - l->ref_range.first;
        processRange(l->ref_range.first, l->seq_range.first, len, smallest, previous, lcs);
        l->ref_range.first += len; l->seq_range.first += len;
      }
      else if(l->ref_range.first > r->ref_range.first)
      {
        uint64_t len = std::min(r->ref_range.second, l->ref_range.first) - r->ref_range.first;
        processRange(r->ref_range.first + ref_len, r->seq_range.first, len, smallest, previous, lcs);
        r->ref_range.first += len; r->seq_range.first += len;
      }
      else
      {
        // Process the larger values first, because they cannot interfere with the smaller values.
        // This works because the ranges of left and right matches cannot overlap.
        uint64_t len = std::min(length(l->ref_range), length(r->ref_range)) - 1;
        if(l->seq_range.first > r->seq_range.first)
        {
          processRange(l->ref_range.first, l->seq_range.first, len, smallest, previous, lcs);
          processRange(r->ref_range.first + ref_len, r->seq_range.first, len, smallest, previous, lcs);
        }
        else
        {
          processRange(r->ref_range.first + ref_len, r->seq_range.first, len, smallest, previous, lcs);
          processRange(l->ref_range.first, l->seq_range.first, len, smallest, previous, lcs);
        }
        l->ref_range.first += len; l->seq_range.first += len;
        r->ref_range.first += len; r->seq_range.first += len;
      }
      if(l->ref_range.first >= l->ref_range.second) { ++l; }
      if(r->ref_range.first >= r->ref_range.second) { ++r; }
    }
    util::clear(left_matches); util::clear(right_matches);
  }

  // Fill the bitvectors.
  util::assign(ref_lcs, bit_vector(ref_len, 0)); util::assign(seq_lcs, bit_vector(seq_len, 0));
  uint64_t pos = smallest[2 * lcs], val = smallest[2 * lcs + 1];
  for(uint64_t i = lcs; i > 0; i--)
  {
    ref_lcs[pos % ref_len] = 1; seq_lcs[val] = 1;
    val = previous[2 * pos + 1]; pos = previous[2 * pos];
  }

  return lcs;
}

//------------------------------------------------------------------------------

} // namespace relative
