#include <cstdlib>
#include <fstream>

#include <sys/resource.h>

#include <sdsl/qsufsort.hpp>

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

const uint64_t SA_SAMPLE_RATE     = 128;
const uint64_t DEFAULT_BLOCK_SIZE = 64 * 1048576;

inline uint64_t
_LF1(bit_vector::rank_1_type& rank, uint64_t pos, uint64_t zeros)
{
  return zeros + 1 + rank(pos);
}

// The endmarker is encoded as an 0-bit, which should be ignored in rank_0.
inline uint64_t
_LF0(bit_vector::rank_1_type& rank, uint64_t pos, uint64_t endmarker)
{
  return 1 + pos - rank(pos) - (pos > endmarker ? 1 : 0);
}

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
  uint64_t first_pos = offset;
  bwt[offset] = buffer[block_size - 2] - 1;
  for(uint64_t i = 0; i < sa.size(); i++)
  {
    if(sa[i] == 0) { first_pos = offset + i + 1; bwt[offset + i + 1] = 0; }
    else { bwt[offset + i + 1] = buffer[sa[i] - 1] - 1; }
  }
  delete[] buffer; buffer = 0;

  return first_pos;
}

/*
  Input:    bwt[0, offset + block_size - 1] contains a prefix of a binary sequence.
            bwt[offset + block_size, n] contains the BWT of the suffix.
            The endmarker is stored at bwt[first_pos].
            No sanity checks are done for the input (call incrementalBWT() instead).
  Output:   Prefix bwt[0, offset - 1] remains intact.
            bwt[offset, n] contains the BWT of the suffix.
            The endmarker is encoded by a 0-bit.
  Returns:  Position in bwt that contains the endmarker.
*/
uint64_t
prevBWTBlock(bit_vector& bwt, uint64_t offset, uint64_t block_size, uint64_t first_pos)
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
  uint64_t text_pos = bwt_start, bwt_pos = first_pos;
  while(text_pos > offset)
  {
    text_pos--;
    if(bwt[text_pos]) { bwt_pos = _LF1(bwt_rank, bwt_pos + 1, zeros) - 1 + bwt_start - false_ones; }
    else              { bwt_pos = _LF0(bwt_rank, bwt_pos + 1, first_pos) - 1 + bwt_start - false_zeros; }
    ra[text_pos - offset] = bwt_pos - bwt_start;
  }

  // Build SA using the rank array.
  for(uint64_t i = 0; i < block_size; i++)
  {
    // Add the character to get the correct sorting order, 1 to handle the endmarker,
    // and 0 or 1 to handle first_pos.
    ra[i] += bwt[i + offset] + 1 + (ra[i] >= first_pos ? 1 : 0);
  }
  ra[block_size] = first_pos - bwt_start + (first_pos >= bwt_start + zeros + 1 ? 1 : 0) + 1;
  ra[block_size + 1] = 0;
  int_vector<64> sa(ra.size(), 0);
  qsufsort::sorter<int_vector<64> > sorter;
  sorter.do_sort(sa, ra); // The rank array gets overwritten.
  util::clear(ra);

  // Build BWT and determine the new first_pos.
  bit_vector increment(block_size);
  uint64_t sa_offset = 1, new_first_pos = 0;
  for(uint64_t i = 1; i < sa.size(); i++)
  {
    if(sa[i] == 0) { new_first_pos = i - sa_offset; increment[i - sa_offset] = 0; }
    else if(sa[i] == block_size) { sa_offset++; } // Skip the next suffix.
    else { increment[i - sa_offset] = bwt[offset + sa[i] - 1]; }
  }
  util::clear(sa);

  // Rebuild the rank array.
  ra.resize(block_size);
  text_pos = bwt_start; bwt_pos = first_pos;
  while(text_pos > offset)
  {
    text_pos--;
    if(bwt[text_pos]) { bwt_pos = _LF1(bwt_rank, bwt_pos + 1, zeros) - 1 + bwt_start - false_ones; }
    else              { bwt_pos = _LF0(bwt_rank, bwt_pos + 1, first_pos) - 1 + bwt_start - false_zeros; }
    ra[text_pos - offset] = bwt_pos - bwt_start;
  }

  // Merge the BWTs.
  bwt[first_pos] = bwt[offset + block_size - 1];  // No longer an endmarker.
  std::sort(ra.begin(), ra.begin() + block_size);
  uint64_t inc_pos = 0, old_pos = 0;
  while(inc_pos < block_size && old_pos < bwt_size)
  {
    while(old_pos <= ra[inc_pos])
    {
      bwt[offset + inc_pos + old_pos] = bwt[bwt_start + old_pos]; old_pos++;
    }
    if(inc_pos == new_first_pos) { first_pos = offset + inc_pos + old_pos; }
    bwt[offset + inc_pos + old_pos] = increment[inc_pos]; inc_pos++;
  }
  while(inc_pos < block_size)
  {
    if(inc_pos == new_first_pos) { first_pos = offset + inc_pos + old_pos; }
    bwt[offset + inc_pos + old_pos] = increment[inc_pos]; inc_pos++;
  }
  // Note that the remaining bits of the old BWT are already in place.

  return first_pos;
}

/*
  Input:    Bit sequence bwt[0, n-1], with bwt[n] = 0 as an endmarker.
  Output:   BWT(bwt + $) stored in the same bit sequence. The endmarker is encoded with a 0-bit.
  Returns:  Position in bwt that contains the endmarker.
*/
uint64_t
incrementalBWT(bit_vector& bwt, uint64_t block_size = DEFAULT_BLOCK_SIZE)
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
  uint64_t first_pos = lastBWTBlock(bwt, offset);

  // Handle the rest of the blocks.
  while(offset > 0)
  {
    offset -= block_size;
    first_pos = prevBWTBlock(bwt, offset, block_size, first_pos);
  }

  return first_pos;
}

/*
  Input:    BWT, the position of the endmarker, and the number of real 0-bits in the BWT.
  Output:   sa_samples contains SA[i * sample_rate + 1] for all i.
*/
void
sampleSA(bit_vector& bwt, bit_vector::rank_1_type& bwt_rank, uint64_t first_pos, uint64_t zeros,
  int_vector<0>& sa_samples, uint64_t sample_rate)
{
#ifdef VERBOSE_OUTPUT
  std::cout << "Sampling SA... "; std::cout.flush();
#endif
  if(bwt.size() <= 1 || sample_rate == 0) { return; }

  sa_samples.width(bitlength(bwt.size() - 1));
  sa_samples.resize((bwt.size() + sample_rate - 2) / sample_rate);
  uint64_t bwt_pos = 0, text_pos = bwt.size() - 1;
  while(text_pos > 0)
  {
    text_pos--;
    if(bwt[bwt_pos]) { bwt_pos = _LF1(bwt_rank, bwt_pos, zeros); }
    else             { bwt_pos = _LF0(bwt_rank, bwt_pos, first_pos); }
    if(bwt_pos % sample_rate == 1) { sa_samples[bwt_pos / sample_rate] = text_pos; }
  }
#ifdef VERBOSE_OUTPUT
  std::cout << "done!" << std::endl;
#endif
}

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
  uint64_t first_pos = incrementalBWT(bwt);
  bit_vector::rank_1_type bwt_rank(&bwt);

  // Sample the SA.
  int_vector<0> sa_samples;
  uint64_t sample_rate = SA_SAMPLE_RATE;
  sampleSA(bwt, bwt_rank, first_pos, zeros, sa_samples, sample_rate);

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
        new_range.first = _LF0(bwt_rank, range.first, first_pos);
        new_range.second = _LF0(bwt_rank, range.second + 1, first_pos) - 1;
      }
      if(isEmpty(new_range)) { break; }
      else { range = new_range; len++; }
    }
    uint64_t reverse_pos = 0; // Position of the reverse pattern in reverse refefence.
    while(true)
    {
      if(range.first == first_pos) { break; }
      if(range.first % sample_rate == 1) { reverse_pos += sa_samples[range.first / sample_rate]; break; }
      range.first = (bwt[range.first] ? _LF1(bwt_rank, range.first, zeros) : _LF0(bwt_rank, range.first, first_pos));
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
