#include <cstdlib>
#include <fstream>

#include <sdsl/qsufsort.hpp>

#include "rlz.h"
#include "utils.h"

//------------------------------------------------------------------------------

void relativeLZ(const bit_vector& text, const bit_vector& reference,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, bit_vector& mismatches)
{
  // Handle the reference.
  uint64_t ones = util::cnt_one_bits(reference);
  uint64_t zeros = reference.size() - ones;
  if(ones == 0 || zeros == 0)
  {
    std::cerr << "relativeLZ(): Reference must contain both 0-bits and 1-bits!" << std::endl;
    return;
  }
  unsigned char* buffer = new unsigned char[reference.size() + 1];
  for(uint64_t i = 0; i < reference.size(); i++)
  {
    buffer[i] = reference[i] + 1;
  }
  buffer[reference.size()] = 0;

  // Build SA.
  int_vector<> sa(reference.size(), 0, bitlength(reference.size()));
  algorithm::calculate_sa(buffer, reference.size(), sa);
  delete[] buffer; buffer = 0;  

  // Parse the text.
  uint64_t pos = 0;
  starts.clear(); lengths.clear();
  std::vector<bool> char_buffer;
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
      if(sa[low] + len >= reference.size()) { break; }  // We use the last matched char as the mismatch.
      if(reference[sa[low] + len] != text[pos + len])
      {
        len++; break; // We use the mismatch.
      }
      range.first = low; high = last_high;
      while(low < high) // Upper bound for pattern text[pos, pos + len].
      {
        uint64_t mid = low + (high - low + 1) / 2;
        if(reference[sa[mid] + len] == text[pos + len]) { low = mid; }
        else { high = mid - 1; }
      }
      range.second = high; len++;
    }
    starts.push_back(sa[range.first]); lengths.push_back(len);
    pos += len;
    char_buffer.push_back(text[pos - 1]);
  }
  mismatches.resize(char_buffer.size());
  for(uint64_t i = 0; i < char_buffer.size(); i++) { mismatches[i] = char_buffer[i]; }
}

//------------------------------------------------------------------------------

void
relativeLZSuccinct(const bit_vector& text, const bit_vector& reference,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, bit_vector& mismatches)
{
  if(text.size() == 0) { return; }

  // Handle the reference.
  uint64_t ones = util::cnt_one_bits(reference);
  uint64_t zeros = reference.size() - ones;
  if(ones == 0 || zeros == 0)
  {
    std::cerr << "relativeLZSuccinct(): Reference must contain both 0-bits and 1-bits!" << std::endl;
    return;
  }
  bit_vector bwt(reference.size() + 1);
  for(uint64_t i = 0; i < reference.size(); i++) { bwt[i] = reference[reference.size() - i - 1]; }
  bwt[reference.size()] = 0;

  // Build BWT for the reverse reference.
  uint64_t endmarker = incrementalBWT(bwt);
  bit_vector::rank_1_type bwt_rank(&bwt);

  // Sample the SA.
  int_vector<0> sa_samples;
  uint64_t sample_rate = RLZ_SA_SAMPLE_RATE;
  sampleSA(bwt, bwt_rank, endmarker, sa_samples, sample_rate);

  // Use backward searching on the BWT of the reverse reference to parse the text.
  uint64_t text_pos = 0;
  starts.clear(); lengths.clear();
  std::vector<bool> char_buffer;
  while(text_pos < text.size())
  {
    range_type range = (text[text_pos] ? range_type(zeros + 1, reference.size()) : range_type(1, zeros));
    uint64_t len = 1; // We have matched len bits in the reverse reference.
    while(text_pos + len < text.size())
    {
      range_type new_range;
      if(text[text_pos + len])
      {
        new_range.first = _LF1(bwt_rank, range.first, zeros);
        new_range.second = _LF1(bwt_rank, range.second + 1, zeros) - 1;
      }
      else
      {
        new_range.first = _LF0(bwt_rank, range.first, endmarker);
        new_range.second = _LF0(bwt_rank, range.second + 1, endmarker) - 1;
      }
      if(isEmpty(new_range)) { break; }
      else { range = new_range; len++; }
    }
    uint64_t reverse_pos = 0; // Position of the reverse pattern in reverse reference.
    while(true)
    {
      if(range.first == endmarker) { break; }
      if(range.first % sample_rate == 1) { reverse_pos += sa_samples[range.first / sample_rate]; break; }
      range.first = (bwt[range.first] ? _LF1(bwt_rank, range.first, zeros) : _LF0(bwt_rank, range.first, endmarker));
      reverse_pos++;
    }
    starts.push_back(reference.size() - reverse_pos - len); // Convert into actual pattern position.
    if(text_pos + len < text.size()) { len++; } // Add the mismatching character.
    lengths.push_back(len);
    char_buffer.push_back(text[text_pos + len - 1]);
    text_pos += len;
  }
  mismatches.resize(char_buffer.size());
  for(uint64_t i = 0; i < char_buffer.size(); i++) { mismatches[i] = char_buffer[i]; }
}

//------------------------------------------------------------------------------

/*
  Input:    Bit sequence bwt[0, n-1], with bwt[n] = 0 as an endmarker.
  Output:   bwt[0, offset - 1] remains intact.
            bwt[offset, n] contains BWT(bwt[offset, n-1] + $).
            The endmarker is encoded with a 0-bit.
  Returns:  Position in bwt that contains the endmarker.
*/
uint64_t
lastBWTBlock(bit_vector& bwt, uint64_t offset)
{
#ifdef VERBOSE_OUTPUT
  std::cout << "Offset: " << offset << std::endl;
#endif
  if(offset + 1 >= bwt.size()) { return offset; }

  // Prepare the text.
  uint64_t block_size = bwt.size() - offset;
  unsigned char* buffer = new unsigned char[block_size];
  for(uint64_t i = offset; i < bwt.size() - 1; i++)
  {
    buffer[i - offset] = bwt[i] + 1;
  }
  buffer[block_size - 1] = 0;

  // Build SA.
  int_vector<> sa(block_size - 1, 0, bitlength(block_size - 1));
  algorithm::calculate_sa(buffer, block_size - 1, sa);

  // Build BWT.
  uint64_t endmarker = offset;
  bwt[offset] = buffer[block_size - 2] - 1;
  for(uint64_t i = 0, j = offset + 1; i < sa.size(); i++, j++)
  {
    if(sa[i] == 0) { endmarker = j; bwt[j] = 0; }
    else { bwt[j] = buffer[sa[i] - 1] - 1; }
  }
  delete[] buffer; buffer = 0;

  return endmarker;
}

/*
  Input:    bwt[0, offset + block_size - 1] contains a prefix of a binary sequence.
            bwt[offset + block_size, n] contains the BWT of the suffix.
            The endmarker is stored at bwt[endmarker].
            No sanity checks are done for the input (call incrementalBWT() instead).
  Output:   Prefix bwt[0, offset - 1] remains intact.
            bwt[offset, n] contains the BWT of the suffix.
            The endmarker is encoded by a 0-bit.
  Returns:  Position in bwt that contains the endmarker.
*/
uint64_t
prevBWTBlock(bit_vector& bwt, uint64_t offset, uint64_t block_size, uint64_t endmarker)
{
#ifdef VERBOSE_OUTPUT
  std::cout << "Offset: " << offset << std::endl;
#endif

  // Prepare the BWT.
  bit_vector::rank_1_type bwt_rank(&bwt);
  uint64_t bwt_start = offset + block_size, bwt_size = bwt.size() - offset - block_size;
  uint64_t false_ones = bwt_rank(bwt_start);      // 1-bits before the BWT.
  uint64_t false_zeros = bwt_start - false_ones;  // 0-bits before the BWT.
  uint64_t ones = bwt_rank(bwt.size()) - false_ones;
  uint64_t zeros = bwt_size - ones - 1;

  // Build the rank array.
  int_vector<64> ra(block_size + 2, 0);
  uint64_t text_pos = bwt_start, bwt_pos = endmarker;
  while(text_pos > offset)
  {
    text_pos--;
    if(bwt[text_pos]) { bwt_pos = _LF1(bwt_rank, bwt_pos + 1, zeros) - 1 + bwt_start - false_ones; }
    else              { bwt_pos = _LF0(bwt_rank, bwt_pos + 1, endmarker) - 1 + bwt_start - false_zeros; }
    ra[text_pos - offset] = bwt_pos - bwt_start;
  }

  // Build SA using the rank array.
  for(uint64_t i = 0; i < block_size; i++)
  {
    // Add the character to get the correct sorting order, 1 to handle the endmarker suffix,
    // and 0 or 1 to handle the endmarker in BWT.
    ra[i] += bwt[i + offset] + 1 + (ra[i] >= endmarker ? 1 : 0);
  }
  ra[block_size] = endmarker - bwt_start + (endmarker >= bwt_start + zeros + 1 ? 1 : 0) + 1;
  ra[block_size + 1] = 0;
  int_vector<64> sa(ra.size(), 0);
  qsufsort::sorter<int_vector<64> > sorter;
  sorter.do_sort(sa, ra); // The rank array gets overwritten.
  util::clear(ra);

  // Build BWT and determine the new endmarker.
  bit_vector increment(block_size);
  uint64_t sa_offset = 1, new_endmarker = 0;
  for(uint64_t i = 1; i < sa.size(); i++)
  {
    if(sa[i] == 0) { new_endmarker = i - sa_offset; increment[i - sa_offset] = 0; }
    else if(sa[i] == block_size) { sa_offset++; } // Skip the next suffix.
    else { increment[i - sa_offset] = bwt[offset + sa[i] - 1]; }
  }
  util::clear(sa);

  // Rebuild the rank array that was overwritten during SA construction.
  ra.resize(block_size);
  text_pos = bwt_start; bwt_pos = endmarker;
  while(text_pos > offset)
  {
    text_pos--;
    if(bwt[text_pos]) { bwt_pos = _LF1(bwt_rank, bwt_pos + 1, zeros) - 1 + bwt_start - false_ones; }
    else              { bwt_pos = _LF0(bwt_rank, bwt_pos + 1, endmarker) - 1 + bwt_start - false_zeros; }
    ra[text_pos - offset] = bwt_pos;  // No need to correct by offset anymore.
  }

  // Merge the BWTs.
  bwt[endmarker] = bwt[offset + block_size - 1];  // No longer an endmarker.
  std::sort(ra.begin(), ra.end());
  uint64_t inc_pos = 0, old_pos = bwt_start, next_pos = offset;
  while(inc_pos < block_size)
  {
    while(old_pos <= ra[inc_pos])
    {
      uint64_t bits = std::min((uint64_t)64, ra[inc_pos] + 1 - old_pos);
      uint64_t temp = bwt.get_int(old_pos, bits); old_pos += bits;
      bwt.set_int(next_pos, temp, bits); next_pos += bits;
//        bwt[next_pos++] = bwt[old_pos++];
    }
    if(inc_pos == new_endmarker) { endmarker = next_pos; }
    bwt[next_pos++] = increment[inc_pos++];
  }
  // Note that the remaining bits of the old BWT are already in place.

  return endmarker;
}

uint64_t
incrementalBWT(bit_vector& bwt, uint64_t block_size)
{
  if(block_size < 2)
  {
    std::cerr << "incrementalBWT(): Block size must be at least 2!" << std::endl;
    return 0;
  }

  // Handle small inputs.
  if(bwt.size() <= block_size) { return lastBWTBlock(bwt, 0); }

  // Handle the last block.
  uint64_t offset = (bwt.size() / block_size) * block_size;
  if(offset >= bwt.size() - 1) { offset -= block_size; }
  uint64_t endmarker = lastBWTBlock(bwt, offset);

  // Handle the rest of the blocks.
  while(offset > 0)
  {
    offset -= block_size;
    endmarker = prevBWTBlock(bwt, offset, block_size, endmarker);
  }

  return endmarker;
}

void
sampleSA(bit_vector& bwt, bit_vector::rank_1_type& bwt_rank, uint64_t endmarker,
  int_vector<0>& sa_samples, uint64_t sample_rate)
{
  if(bwt.size() <= 1 || sample_rate == 0) { return; }
#ifdef VERBOSE_OUTPUT
  std::cout << "Sampling SA... "; std::cout.flush();
#endif

  uint64_t zeros = bwt.size() - util::cnt_one_bits(bwt) - 1;
  sa_samples.width(bitlength(bwt.size() - 1));
  sa_samples.resize((bwt.size() + sample_rate - 2) / sample_rate);
  uint64_t bwt_pos = 0, text_pos = bwt.size() - 1;
  while(text_pos > 0)
  {
    text_pos--;
    if(bwt[bwt_pos]) { bwt_pos = _LF1(bwt_rank, bwt_pos, zeros); }
    else             { bwt_pos = _LF0(bwt_rank, bwt_pos, endmarker); }
    if(bwt_pos % sample_rate == 1) { sa_samples[bwt_pos / sample_rate] = text_pos; }
  }
#ifdef VERBOSE_OUTPUT
  std::cout << "done!" << std::endl;
#endif
}

//------------------------------------------------------------------------------
