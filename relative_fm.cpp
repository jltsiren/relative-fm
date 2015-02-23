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

  for(uint64_t i = 0, found = 0; i < length(ref_range) && found < max_f; i++)
  {
#ifdef _OPENMP
    if(ref_buffer[i] == c) { ref_lcs[i] = 1; found++; }
#else
    if(ref_buffer[i] == c) { ref_lcs[ref_range.first + i] = 1; found++; }
#endif
  }
  for(uint64_t i = 0, found = 0; i < length(seq_range) && found < max_f; i++)
  {
#ifdef _OPENMP
    if(seq_buffer[i] == c) { seq_lcs[i] = 1; found++; }
#else
    if(seq_buffer[i] == c) { seq_lcs[seq_range.first + i] = 1; found++; }
#endif
  }

#ifdef VERBOSE_STATUS_INFO
#ifdef _OPENMP
  #pragma omp critical (stderr)
#endif
  {
    std::cerr << "  Matched " << max_f << " copies of character " << c << " (" << ((char)c) << ")";
    std::cerr << " from sequences of length " << range_type(ref_buffer.size(), seq_buffer.size()) << std::endl;
  }
#endif
  return max_f;
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
