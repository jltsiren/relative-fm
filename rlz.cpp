#include <cstdlib>
#include <fstream>

#include <sdsl/construct_sa.hpp>

#include "rlz.h"
#include "utils.h"

namespace relative
{

//------------------------------------------------------------------------------

void
relativeLZSuccinct(const bit_vector& text, const bit_vector& reference,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, bit_vector& mismatches)
{
  if(text.size() == 0) { return; }

#ifdef VERBOSE_STATUS_INFO
  std::cout << "RLZ parsing: text length " << text.size() << ", reference length " << reference.size() << "." << std::endl;
#endif

  bv_fmi fmi(reference);
  relativeLZSuccinct(text, fmi, starts, lengths, mismatches);
}

void
relativeLZSuccinct(const bit_vector& text, const bv_fmi& reference,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, bit_vector& mismatches)
{
  if(text.size() == 0) { return; }

  // Use backward searching on the BWT of the reverse reference to parse the text.
  uint64_t text_pos = 0;
  starts.clear(); lengths.clear();
  std::vector<bool> char_buffer;
  while(text_pos < text.size())
  {
    bv_fmi::range_type range = reference.bitRange(text[text_pos]);
    uint64_t len = 1; // We have matched len bits in the reverse reference.
    while(text_pos + len < text.size())
    {
      bv_fmi::range_type new_range;
      if(text[text_pos + len])
      {
        new_range.first = reference.LF1(range.first);
        new_range.second = reference.LF1(range.second + 1) - 1;
      }
      else
      {
        new_range.first = reference.LF0(range.first);
        new_range.second = reference.LF0(range.second + 1) - 1;
      }
      if(new_range.first > new_range.second) { break; }
      else { range = new_range; len++; }
    }
    uint64_t reverse_pos = 0; // Position of the reverse pattern in reverse reference.
    while(true)
    {
      if(range.first == reference.endmarker)
      {
        break;
      }
      if(range.first % reference.sample_rate == 1)
      {
        reverse_pos += reference.sampleAt(range.first);
        break;
      }
      range.first = (reference.bwt[range.first] ? reference.LF1(range.first) : reference.LF0(range.first));
      reverse_pos++;
    }
    starts.push_back(reference.bwt.size() - reverse_pos - len - 1); // Convert into actual pattern position.
    if(text_pos + len < text.size()) { len++; } // Add the mismatching character.
    lengths.push_back(len);
    char_buffer.push_back(text[text_pos + len - 1]);
    text_pos += len;
  }
  mismatches.resize(char_buffer.size());
  for(uint64_t i = 0; i < char_buffer.size(); i++) { mismatches[i] = char_buffer[i]; }

#ifdef VERBOSE_STATUS_INFO
  std::cout << "Parsed the text as " << starts.size() << " phrases." << std::endl;
#endif
}

//------------------------------------------------------------------------------

void
relativeLZ(const int_vector<8>& text, const int_vector<8>& reference,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, int_vector<8>& mismatches)
{
  if(text.size() == 0) { return; }

#ifdef VERBOSE_STATUS_INFO
  std::cout << "RLZ parsing: text length " << text.size() << ", reference length " << reference.size() << "." << std::endl;
#endif

  csa_wt<> csa;
  reverseIndex(reference, csa);
  relativeLZ(text, csa, starts, lengths, mismatches);
}

void
relativeLZ(const int_vector<0>& text, const int_vector<0>& reference,
  std::vector<uint64_t>& starts, std::vector<uint64_t>& lengths, int_vector<0>& mismatches)
{
  if(text.size() == 0) { return; }

#ifdef VERBOSE_STATUS_INFO
  std::cout << "RLZ parsing: text length " << text.size() << ", reference length " << reference.size() << "." << std::endl;
#endif

  csa_wt<wm_int<> > csa;
  reverseIndex(reference, csa);
  relativeLZ(text, csa, starts, lengths, mismatches);
}

//------------------------------------------------------------------------------

bv_fmi::bv_fmi(const bit_vector& source, uint64_t block_size, uint64_t _sample_rate)
{
  if(block_size < 2)
  {
    std::cerr << "bv_fmi::bv_fmi(): Block size must be at least 2" << std::endl;
    return;
  }

  // Copy and reverse the source.
  this->zeros = source.size() - util::cnt_one_bits(source);
  if(this->zeros == 0 || this->zeros == source.size())
  {
    std::cerr << "bv_fmi::bv_fmi(): Source must contain both 0-bits and 1-bits" << std::endl;
    return;
  }
  this->bwt.resize(source.size() + 1);
  for(uint64_t i = 0; i < source.size(); i++) { this->bwt[i] = source[source.size() - i - 1]; }
  this->bwt[source.size()] = 0;

  // Build BWT for the reverse source.
  this->incrementalBWT(block_size);
  util::clear(this->rank); util::init_support(this->rank, &(this->bwt));

  // Sample the SA.
  this->sample_rate = _sample_rate;
  this->sampleSA();
}

bv_fmi::bv_fmi(std::istream& in)
{
  this->bwt.load(in);
  this->rank.load(in, &(this->bwt));
  read_member(this->zeros, in);
  read_member(this->endmarker, in);
  this->sa_samples.load(in);
  read_member(this->sample_rate, in);
}

uint64_t
bv_fmi::serialize(std::ostream& out)
{
  uint64_t bytes = 0;
  bytes += this->bwt.serialize(out);
  bytes += this->rank.serialize(out);
  bytes += write_member(this->zeros, out);
  bytes += write_member(this->endmarker, out);
  bytes += this->sa_samples.serialize(out);
  bytes += write_member(this->sample_rate, out);
  return bytes;
}

void
bv_fmi::incrementalBWT(uint64_t block_size)
{
  // Handle small inputs.
  if(this->bwt.size() <= block_size)
  {
    this->lastBWTBlock(0);
  }

  // Handle the last block.
  uint64_t offset = (this->bwt.size() / block_size) * block_size;
  if(offset + 2 >= this->bwt.size()) { offset -= block_size; }
  lastBWTBlock(offset);

  // Handle the rest of the blocks.
  while(offset > 0)
  {
    offset -= block_size;
    prevBWTBlock(offset, block_size);
  }
}

/*
  Input:    Bit sequence bwt[0, n-1], with bwt[n] = 0 as an endmarker.
  Output:   bwt[0, offset - 1] remains intact.
            bwt[offset, n] contains BWT(bwt[offset, n-1] + $).
            The endmarker is encoded with a 0-bit and its position is stored.
*/
void
bv_fmi::lastBWTBlock(uint64_t offset)
{
#ifdef VERBOSE_STATUS_INFO
  std::cout << "Offset: " << offset << std::endl;
#endif
  if(offset + 1 >= this->bwt.size()) { this->endmarker = offset; return; }

  // Prepare the text.
  uint64_t block_size = this->bwt.size() - offset;
  unsigned char* buffer = new unsigned char[block_size];
  for(uint64_t i = offset; i < bwt.size() - 1; i++)
  {
    buffer[i - offset] = this->bwt[i] + 1;
  }
  buffer[block_size - 1] = 0;

  // Build SA.
  int_vector<> sa(block_size - 1, 0, bits::hi(block_size - 1) + 1);
  algorithm::calculate_sa(buffer, block_size - 1, sa);

  // Build BWT.
  bwt[offset] = buffer[block_size - 2] - 1;
  for(uint64_t i = 0, j = offset + 1; i < sa.size(); i++, j++)
  {
    if(sa[i] == 0) { this->endmarker = j; this->bwt[j] = 0; }
    else { this->bwt[j] = buffer[sa[i] - 1] - 1; }
  }
  delete[] buffer; buffer = 0;
}

/*
  Input:    bwt[0, offset + block_size - 1] contains a prefix of a binary sequence.
            bwt[offset + block_size, n] contains the BWT of the suffix.
            The endmarker is stored at bwt[endmarker].
  Output:   Prefix bwt[0, offset - 1] remains intact.
            bwt[offset, n] contains the BWT of the suffix.
            The endmarker is updated and encoded by a 0-bit.
*/
void
bv_fmi::prevBWTBlock(uint64_t offset, uint64_t block_size)
{
#ifdef VERBOSE_STATUS_INFO
  std::cout << "Offset: " << offset << std::endl;
#endif

  // Prepare the BWT.
  util::clear(this->rank); util::init_support(this->rank, &(this->bwt));
  uint64_t bwt_start = offset + block_size;
  uint64_t false_ones = this->rank(bwt_start);    // 1-bits before the BWT.
  uint64_t false_zeros = bwt_start - false_ones;  // 0-bits before the BWT.

  // Build the rank array.
  int_vector<64> ra(block_size + 2, 0);
  uint64_t text_pos = bwt_start, bwt_pos = this->endmarker;
  while(text_pos > offset)
  {
    text_pos--;
    if(this->bwt[text_pos]) { bwt_pos = this->LF1(bwt_pos + 1) - 1; }
    else                    { bwt_pos = this->LF0(bwt_pos + 1) - 1 + false_ones; }
    ra[text_pos - offset] = bwt_pos - bwt_start;
  }

  // Build SA using the rank array.
  for(uint64_t i = 0; i < block_size; i++)
  {
    // Add the character to get the correct sorting order, 1 to handle the endmarker suffix,
    // and 0 or 1 to handle the endmarker in BWT.
    ra[i] += this->bwt[i + offset] + 1 + (ra[i] >= this->endmarker - bwt_start ? 1 : 0);
  }
  ra[block_size] = this->endmarker - bwt_start +
    (this->endmarker >= bwt_start + this->zeros - false_zeros + 1 ? 1 : 0) + 1;
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
    else { increment[i - sa_offset] = this->bwt[offset + sa[i] - 1]; }
  }
  util::clear(sa);

  // Rebuild the rank array that was overwritten during SA construction.
  ra.resize(block_size);
  text_pos = bwt_start; bwt_pos = this->endmarker;
  while(text_pos > offset)
  {
    text_pos--;
    if(this->bwt[text_pos]) { bwt_pos = this->LF1(bwt_pos + 1) - 1; }
    else                    { bwt_pos = this->LF0(bwt_pos + 1) - 1 + false_ones; }
    ra[text_pos - offset] = bwt_pos;  // No need to correct by offset anymore.
  }

  // Merge the BWTs.
  this->bwt[this->endmarker] = this->bwt[offset + block_size - 1];  // No longer an endmarker.
  sequentialSort(ra.begin(), ra.end());
  uint64_t inc_pos = 0, old_pos = bwt_start, next_pos = offset;
  while(inc_pos < block_size)
  {
    while(old_pos <= ra[inc_pos])
    {
      uint64_t bits = std::min((uint64_t)64, ra[inc_pos] + 1 - old_pos);
      uint64_t temp = this->bwt.get_int(old_pos, bits); old_pos += bits;
      this->bwt.set_int(next_pos, temp, bits); next_pos += bits;
//        bwt[next_pos++] = bwt[old_pos++];
    }
    if(inc_pos == new_endmarker) { this->endmarker = next_pos; }
    this->bwt[next_pos++] = increment[inc_pos++];
  }
  // Note that the remaining bits of the old BWT are already in place.
}

void
bv_fmi::sampleSA()
{
  if(this->bwt.size() <= 1 || this->sample_rate == 0) { return; }
#ifdef VERBOSE_STATUS_INFO
  std::cout << "Sampling SA... "; std::cout.flush();
#endif

  this->sa_samples.width(bits::hi(this->bwt.size() - 1) + 1);
  this->sa_samples.resize((this->bwt.size() + sample_rate - 2) / sample_rate);
  uint64_t bwt_pos = 0, text_pos = this->bwt.size() - 1;
  while(text_pos > 0)
  {
    text_pos--;
    if(bwt[bwt_pos]) { bwt_pos = this->LF1(bwt_pos); }
    else             { bwt_pos = this->LF0(bwt_pos); }
    if(bwt_pos % this->sample_rate == 1) { this->sa_samples[bwt_pos / this->sample_rate] = text_pos; }
  }
#ifdef VERBOSE_STATUS_INFO
  std::cout << "done." << std::endl;
#endif
}

//------------------------------------------------------------------------------

relative_encoder::relative_encoder()
{
}

relative_encoder::relative_encoder(const relative_encoder& r)
{
  this->copy(r);
}

relative_encoder::relative_encoder(relative_encoder&& r)
{
  *this = std::move(r);
}

relative_encoder&
relative_encoder::operator=(const relative_encoder& r)
{
  if(this != &r) { this->copy(r); }
  return *this;
}

relative_encoder&
relative_encoder::operator=(relative_encoder&& r)
{
  if(this != &r)
  {
    this->values = std::move(r.values);
    this->rle = std::move(r.rle);
    this->rank = std::move(r.rank); this->rank.set_vector(&(this->rle));
  }
  return *this;
}

void
relative_encoder::copy(const relative_encoder& r)
{
  this->values = r.values;
  this->rle = r.rle;
  this->rank = r.rank; this->rank.set_vector(&(this->rle));
}

uint64_t
relative_encoder::reportSize() const
{
  return size_in_bytes(this->values) + size_in_bytes(this->rle) + size_in_bytes(this->rank);
}

void
relative_encoder::load(std::istream& input)
{
  this->values.load(input);
  this->rle.load(input);
  this->rank.load(input, &(this->rle));
}

uint64_t
relative_encoder::serialize(std::ostream& output) const
{
  uint64_t bytes = 0;
  bytes += this->values.serialize(output);
  bytes += this->rle.serialize(output);
  bytes += this->rank.serialize(output);
  return bytes;
}

//------------------------------------------------------------------------------

} // namespace relative
