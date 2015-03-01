#include "relative_fm.h"

namespace relative
{

//------------------------------------------------------------------------------

range_type
mostFrequentChar(std::vector<uint8_t>& ref_buffer, std::vector<uint8_t>& seq_buffer,
  bit_vector& ref_lcs, bit_vector& seq_lcs,
  range_type ref_range, range_type seq_range)
{
  std::vector<range_type> freqs(256, range_type(0, 0));
  for(uint64_t i = 0; i < ref_buffer.size(); i++) { freqs[ref_buffer[i]].first++; }
  for(uint64_t i = 0; i < seq_buffer.size(); i++) { freqs[seq_buffer[i]].second++; }
  uint64_t max_c = 0, max_f = 0;
  for(uint64_t c = 0; c < freqs.size(); c++)
  {
    uint64_t temp = std::min(freqs[c].first, freqs[c].second);
    if(temp > max_f) { max_c = c; max_f = temp; }
  }

#ifdef _OPENMP
  uint64_t ref_padding = ref_range.first % 64, seq_padding = seq_range.first % 64;
#else
  uint64_t ref_padding = ref_range.first, seq_padding = seq_range.first;
#endif
  for(uint64_t i = 0, found = 0; i < length(ref_range) && found < max_f; i++)
  {
    if(ref_buffer[i] == max_c) { ref_lcs[ref_padding + i] = 1; found++; }
  }
  for(uint64_t i = 0, found = 0; i < length(seq_range) && found < max_f; i++)
  {
    if(seq_buffer[i] == max_c) { seq_lcs[seq_padding + i] = 1; found++; }
  }

#ifdef VERBOSE_STATUS_INFO
  #ifdef _OPENMP
  #pragma omp critical (stderr)
  #endif
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
    std::cerr << " from sequences of length " << range_type(ref_buffer.size(), seq_buffer.size()) << std::endl;
  }
#endif
  return range_type(max_f, std::min(length(ref_range), length(seq_range)) - max_f);
}

range_type
matchRuns(std::vector<uint8_t>& ref_buffer, std::vector<uint8_t>& seq_buffer,
  bit_vector& ref_lcs, bit_vector& seq_lcs,
  range_type ref_range, range_type seq_range,
  const uint8_t* alphabet, uint64_t sigma)
{
  uint64_t ref_pos = 0, seq_pos = 0, found = 0;
  range_type ref_run(0, 0), seq_run(0, 0);
#ifdef _OPENMP
  uint64_t ref_padding = ref_range.first % 64, seq_padding = seq_range.first % 64;
#else
  uint64_t ref_padding = ref_range.first, seq_padding = seq_range.first;
#endif
  while(ref_pos < length(ref_range) && seq_pos < length(seq_range))
  {
    if(ref_run.second <= ref_pos) // Find the next ref_run if needed.
    {
      ref_run.first = ref_buffer[ref_pos]; ref_run.second = ref_pos + 1;
      while(ref_run.second < length(ref_range) && ref_buffer[ref_run.second] == ref_run.first) { ref_run.second++; }
    }
    if(seq_run.second <= seq_pos) // Find the next seq_run if needed.
    {
      seq_run.first = seq_buffer[seq_pos]; seq_run.second = seq_pos + 1;
      while(seq_run.second < length(seq_range) && seq_buffer[seq_run.second] == seq_run.first) { seq_run.second++; }
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
  #ifdef _OPENMP
  #pragma omp critical (stderr)
  #endif
  {
    std::cerr << "Matching the runs in " << std::make_pair(ref_range, seq_range) << std::endl;
    std::cerr << "  Matched " << found << " characters from sequences of length "
              << range_type(ref_buffer.size(), seq_buffer.size()) << std::endl;
  }
#endif
  return range_type(found, std::min(length(ref_range), length(seq_range)) - found);
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

LCS::LCS()
{
  this->lcs_size = 0;
}

LCS::LCS(const bit_vector& a, const bit_vector& b, LCS::size_type _lcs_size)
{
  this->ref = a; this->seq = b;
  this->lcs_size = _lcs_size;

  util::init_support(this->ref_rank, &(this->ref));
#ifndef USE_HYBRID_BITVECTORS
  util::init_support(this->ref_select, &(this->ref));
#endif

  util::init_support(this->seq_rank, &(this->seq));
#ifndef USE_HYBRID_BITVECTORS
  util::init_support(this->seq_select, &(this->seq));
#endif
}

LCS::LCS(const LCS& l)
{
  this->copy(l);
}

LCS::LCS(LCS&& l)
{
  *this = std::move(l);
}

LCS::~LCS()
{
}

void
LCS::copy(const LCS& l)
{
  this->ref = l.ref;
  this->ref_rank = l.ref_rank;
#ifndef USE_HYBRID_BITVECTORS
  this->ref_select = l.ref_select;
#endif

  this->seq = l.seq;
  this->seq_rank = l.seq_rank;
#ifndef USE_HYBRID_BITVECTORS
  this->seq_select = l.seq_select;
#endif

  this->lcs_size = l.lcs_size;

  this->set_vectors();
}

void
LCS::swap(LCS& l)
{
  if(this != &l)
  {
    this->ref.swap(l.ref);
    util::swap_support(this->ref_rank, l.ref_rank, &(this->ref), &(l.ref));
#ifndef USE_HYBRID_BITVECTORS
    util::swap_support(this->ref_select, l.ref_select, &(this->ref), &(l.ref));
#endif

    this->seq.swap(l.seq);
    util::swap_support(this->seq_rank, l.seq_rank, &(this->seq), &(l.seq));
#ifndef USE_HYBRID_BITVECTORS
    util::swap_support(this->seq_select, l.seq_select, &(this->seq), &(l.seq));
#endif

    std::swap(this->lcs_size, l.lcs_size);
  }
}

LCS&
LCS::operator=(const LCS& l)
{
  if(this != &l) { this->copy(l); }
  return *this;
}

LCS&
LCS::operator=(LCS&& l)
{
  if(this != &l)
  {
    this->ref = std::move(l.ref);
    this->ref_rank = std::move(l.ref_rank);
#ifndef USE_HYBRID_BITVECTORS
    this->ref_select = std::move(l.ref_select);
#endif

    this->seq = std::move(l.seq);
    this->seq_rank = std::move(l.seq_rank);
#ifndef USE_HYBRID_BITVECTORS
    this->seq_select = std::move(l.seq_select);
#endif

    this->lcs_size = l.lcs_size;

    this->set_vectors();
  }
  return *this;
}

uint64_t
LCS::serialize(std::ostream& out, structure_tree_node* s, std::string name) const
{
  structure_tree_node* child = structure_tree::add_child(s, name, util::class_name(*this));
  uint64_t written_bytes = 0;

  written_bytes += this->ref.serialize(out, child, "ref");
  written_bytes += this->ref_rank.serialize(out, child, "ref_rank");
#ifndef USE_HYBRID_BITVECTORS
  written_bytes += this->ref_select.serialize(out, child, "ref_select");
#endif

  written_bytes += this->seq.serialize(out, child, "seq");
  written_bytes += this->seq_rank.serialize(out, child, "seq_rank");
#ifndef USE_HYBRID_BITVECTORS
  written_bytes += this->seq_select.serialize(out, child, "seq_select");
#endif

  written_bytes += write_member(this->lcs_size, out, child, "lcs_size");

  structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
LCS::load(std::istream& in)
{
  this->ref.load(in);
  this->ref_rank.load(in, &(this->ref));
#ifndef USE_HYBRID_BITVECTORS
  this->ref_select.load(in, &(this->ref));
#endif

  this->seq.load(in);
  this->seq_rank.load(in, &(this->seq));
#ifndef USE_HYBRID_BITVECTORS
  this->seq_select.load(in, &(this->seq));
#endif

  read_member(this->lcs_size, in);
}

void
LCS::set_vectors()
{
  this->ref_rank.set_vector(&(this->ref));
#ifndef USE_HYBRID_BITVECTORS
  this->ref_select.set_vector(&(this->ref));
#endif

  this->seq_rank.set_vector(&(this->seq));
#ifndef USE_HYBRID_BITVECTORS
  this->seq_select.set_vector(&(this->seq));
#endif
}

#ifdef USE_HYBRID_BITVECTORS
uint64_t
LCS::select(const LCS::vector_type& vec, const LCS::vector_type::rank_1_type& rank, uint64_t i) const
{
  if(i == 0) { return 0; }

  // Find the last position, where rank < i.
  uint64_t low = 0, high = vec.size();
  while(low < high)
  {
    uint64_t mid = low + (high - low + 1) / 2;
    if(rank(mid) >= i) { high = mid - 1; }
    else { low = mid; }
  }

  return low;
}
#endif

//------------------------------------------------------------------------------

} // namespace relative
