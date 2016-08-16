/*
  Copyright (c) 2015, 2016 Genome Research Ltd.
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

#include <map>

#include "relative_fm.h"

namespace relative
{

//------------------------------------------------------------------------------

range_type
mostFrequentChar(const std::vector<uint8_t>& ref_extract, const std::vector<uint8_t>& seq_extract,
  sdsl::bit_vector& ref_lcs, sdsl::bit_vector& seq_lcs,
  range_type ref_range, range_type seq_range)
{
  std::vector<range_type> freqs(256, range_type(0, 0));
  for(size_type i = 0; i < ref_extract.size(); i++) { freqs[ref_extract[i]].first++; }
  for(size_type i = 0; i < seq_extract.size(); i++) { freqs[seq_extract[i]].second++; }
  size_type max_c = 0, max_f = 0;
  for(size_type c = 0; c < freqs.size(); c++)
  {
    size_type temp = std::min(freqs[c].first, freqs[c].second);
    if(temp > max_f) { max_c = c; max_f = temp; }
  }

  size_type ref_padding = ref_range.first % 64, seq_padding = seq_range.first % 64;
  for(size_type i = 0, found = 0; i < Range::length(ref_range) && found < max_f; i++)
  {
    if(ref_extract[i] == max_c) { ref_lcs[ref_padding + i] = 1; found++; }
  }
  for(size_type i = 0, found = 0; i < Range::length(seq_range) && found < max_f; i++)
  {
    if(seq_extract[i] == max_c) { seq_lcs[seq_padding + i] = 1; found++; }
  }

#ifdef VERBOSE_STATUS_INFO
  #pragma omp critical (stderr)
  {
    std::cerr << "Finding the most frequent character in " << std::make_pair(ref_range, seq_range) << std::endl;
    if(max_f == 0)
    {
      for(size_type c = 0; c < freqs.size(); c++)
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
  return range_type(max_f, std::min(Range::length(ref_range), Range::length(seq_range)) - max_f);
}

range_type
matchRuns(const std::vector<uint8_t>& ref_extract, const std::vector<uint8_t>& seq_extract,
  sdsl::bit_vector& ref_lcs, sdsl::bit_vector& seq_lcs,
  range_type ref_range, range_type seq_range,
  const uint8_t* alphabet, size_type sigma)
{
  size_type ref_pos = 0, seq_pos = 0, found = 0;
  range_type ref_run(0, 0), seq_run(0, 0);
  size_type ref_padding = ref_range.first % 64, seq_padding = seq_range.first % 64;
  while(ref_pos < Range::length(ref_range) && seq_pos < Range::length(seq_range))
  {
    if(ref_run.second <= ref_pos) // Find the next ref_run if needed.
    {
      ref_run.first = ref_extract[ref_pos]; ref_run.second = ref_pos + 1;
      while(ref_run.second < Range::length(ref_range) && ref_extract[ref_run.second] == ref_run.first) { ref_run.second++; }
    }
    if(seq_run.second <= seq_pos) // Find the next seq_run if needed.
    {
      seq_run.first = seq_extract[seq_pos]; seq_run.second = seq_pos + 1;
      while(seq_run.second < Range::length(seq_range) && seq_extract[seq_run.second] == seq_run.first) { seq_run.second++; }
    }
    if(ref_run.first == seq_run.first)  // Match the current runs if the characters match.
    {
      size_type run_length = std::min(ref_run.second - ref_pos, seq_run.second - seq_pos);
      for(size_type i = 0; i < run_length; i++)
      {
        ref_lcs[ref_padding + ref_pos + i] = 1; seq_lcs[seq_padding + seq_pos + i] = 1;
      }
      ref_pos = ref_run.second; seq_pos = seq_run.second; found += run_length;
    }
    else  // Otherwise advance the run with the smaller character comp value.
    {
      size_type ref_c = 0, seq_c = 0;
      for(size_type c = 0; c < sigma; c++)
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
  return range_type(found, std::min(Range::length(ref_range), Range::length(seq_range)) - found);
}

// Diagonal i has 2i+3 elements: -(i+1) to (i+1).
inline int offsetFor(int pos, int d) { return (d + 3) * d + pos + 1; }

range_type
greedyLCS(const std::vector<uint8_t>& ref_extract, const std::vector<uint8_t>& seq_extract,
  sdsl::bit_vector& ref_lcs, sdsl::bit_vector& seq_lcs,
  short_record_type& record,
  const uint8_t* alphabet, size_type sigma,
  const align_parameters& parameters,
  std::vector<int>* buf)
{
  range_type ref_range = record.first(), seq_range = record.second();
  bool onlyNs = record.onlyNs(), endmarker = record.endmarker();
  if(Range::empty(ref_range) || Range::empty(seq_range)) { return range_type(0, 0); }
  int ref_len = Range::length(ref_range), seq_len = Range::length(seq_range);

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
  size_type lcs = 0;
  d--;  // The last diagonal.
  size_type ref_padding = ref_range.first % 64, seq_padding = seq_range.first % 64;
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
copyBits(const sdsl::bit_vector& source, sdsl::bit_vector& target, range_type range)
{
  // Copy the first bits.
  size_type offset = source.size() - Range::length(range);
  size_type bits = std::min((size_type)64, source.size()) - offset;
  size_type val = source.get_int(offset, bits);
  target.set_int(range.first, val, bits);

  // Copy the bulk of the bits.
  for(size_type source_pos = 1, target_pos = (range.first + bits) / 64;
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
processRange(size_type pos, size_type val, size_type len,
  sdsl::int_vector<0>& smallest, sdsl::int_vector<0>& previous, size_type& lcs)
{
  // Find the longest increasing subsequence with last value < val.
  size_type low = 1, high = lcs;
  while(low <= high)
  {
    size_type mid = low + (high - low) / 2;
    if(smallest[2 * mid + 1] < val) { low = mid + 1; }
    else { high = mid - 1; }
  }

  // low is now the length of the found subsequence + 1.
  for(size_type i = 0; i < len; i++)
  {
    previous[2 * (pos + i)] = smallest[2 * (low + i) - 2];
    previous[2 * (pos + i) + 1] = smallest[2 * (low + i) - 1];
    smallest[2 * (low + i)] = pos + i;
    smallest[2 * (low + i) + 1] = val + i;
  }
  lcs = std::max(lcs, low + len - 1);
}

void
sortMatches(std::vector<range_pair>& matches, size_type ref_len, size_type seq_len)
{
  range_pair_comparator comp;
  parallelMergeSort(matches.begin(), matches.end(), comp);
  matches.push_back(range_pair(range_type(ref_len, ref_len + 1), range_type(seq_len, seq_len + 1)));
}

size_type
countCoverage(std::vector<range_pair>& left_matches, std::vector<range_pair>& right_matches, size_type ref_len)
{
  size_type total_matches = 0, pos = 0;
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

size_type
increasingSubsequence(std::vector<range_pair>& left_matches, std::vector<range_pair>& right_matches,
  sdsl::bit_vector& ref_lcs, sdsl::bit_vector& seq_lcs, size_type ref_len, size_type seq_len)
{
  if(left_matches.size() == 0 || right_matches.size() == 0)
  {
    sdsl::util::assign(ref_lcs, sdsl::bit_vector(ref_len)); sdsl::util::assign(seq_lcs, sdsl::bit_vector(seq_len));
    return 0;
  }

  /*
    Let smallest(i) = smallest[2 * i, 2 * i + 1] and previous(i) = previous[2 * i, 2 * i + 1].
    smallest(i) = (pos, val) is the ending position of the increasing subsequence of
    length i with the smallest val. previous(pos) = (pos', val') is the predecessor of
    position pos in the longest increasing subsequence pos is a part of. Note that if pos
    is a right match, the actual position is pos - ref_len.
  */
  sdsl::int_vector<0> smallest(2 * (ref_len + 1), 2 * ref_len, bit_length(2 * ref_len));
  sdsl::int_vector<0> previous(4 * ref_len, 2 * ref_len, bit_length(2 * ref_len));
  size_type lcs = 0;
  {
    std::vector<range_pair>::iterator l = left_matches.begin(), r = right_matches.begin();
    while(l->ref_range.first < ref_len || r->ref_range.first < ref_len)
    {
      if(l->ref_range.first < r->ref_range.first)
      {
        size_type len = std::min(l->ref_range.second, r->ref_range.first) - l->ref_range.first;
        processRange(l->ref_range.first, l->seq_range.first, len, smallest, previous, lcs);
        l->ref_range.first += len; l->seq_range.first += len;
      }
      else if(l->ref_range.first > r->ref_range.first)
      {
        size_type len = std::min(r->ref_range.second, l->ref_range.first) - r->ref_range.first;
        processRange(r->ref_range.first + ref_len, r->seq_range.first, len, smallest, previous, lcs);
        r->ref_range.first += len; r->seq_range.first += len;
      }
      else
      {
        // Process the larger values first, because they cannot interfere with the smaller values.
        // This works because the ranges of left and right matches cannot overlap.
        size_type len = std::min(Range::length(l->ref_range), Range::length(r->ref_range)) - 1;
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
    sdsl::util::clear(left_matches); sdsl::util::clear(right_matches);
  }

  // Fill the bitvectors.
  sdsl::util::assign(ref_lcs, sdsl::bit_vector(ref_len, 0)); sdsl::util::assign(seq_lcs, sdsl::bit_vector(seq_len, 0));
  size_type pos = smallest[2 * lcs], val = smallest[2 * lcs + 1];
  for(size_type i = lcs; i > 0; i--)
  {
    ref_lcs[pos % ref_len] = 1; seq_lcs[val] = 1;
    val = previous[2 * pos + 1]; pos = previous[2 * pos];
  }

  return lcs;
}

//------------------------------------------------------------------------------

} // namespace relative
