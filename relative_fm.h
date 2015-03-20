#ifndef _RELATIVE_FM_RELATIVE_FM_H
#define _RELATIVE_FM_RELATIVE_FM_H

#include "simple_fm.h"

namespace relative
{

//------------------------------------------------------------------------------

struct align_parameters
{
  const static uint64_t BLOCK_SIZE      = 1024; // Split the BWTs into blocks of this size or less.
  const static uint64_t MAX_LENGTH      = 32;   // Maximum length of a pattern used to split the BWTs.
  const static uint64_t SEED_LENGTH     = 4;    // Generate all patterns of this length before parallelizing.
  const static int      MAX_D           = 8192; // Number of diagonals to consider in Myers' algorithm.
  const static uint64_t SA_SAMPLE_RATE  = 0;    // Sample rate for relative FM.
  const static uint64_t ISA_SAMPLE_RATE  = 0;   // Sample rate for relative FM.

  // Default size for on demand LCS buffer; enough for d = 100.
  const static uint64_t BUFFER_SIZE = 103 * 101;

  // Default sample rate with a BWT-invariant subsequence.
  const static uint64_t SECONDARY_SA_SAMPLE_RATE = 257;
  const static uint64_t SECONDARY_ISA_SAMPLE_RATE = 512;

  align_parameters() :
    block_size(BLOCK_SIZE), max_length(MAX_LENGTH), seed_length(SEED_LENGTH), max_d(MAX_D),
    sa_sample_rate(SA_SAMPLE_RATE), isa_sample_rate(ISA_SAMPLE_RATE),
    sorted_alphabet(true), preallocate(false), invariant(false)
  {
  }

  uint64_t block_size, max_length, seed_length;
  int      max_d;
  uint64_t sa_sample_rate, isa_sample_rate;
  bool     sorted_alphabet; // comp values are sorted by character values, making partitioning a bit faster.
  bool     preallocate;     // Use preallocated arrays in Myers' algorithm.
  bool     invariant;       // Find a BWT-invariant subsequence that supports relative SA samples.
};

//------------------------------------------------------------------------------

template<class ReferenceBWTType, class SequenceType>
void
getComplement(const ReferenceBWTType& bwt, SequenceType& output, const bit_vector& positions, uint64_t n);

template<class ReferenceType>
void
alignBWTs(const ReferenceType& ref, const ReferenceType& seq,
  bit_vector& ref_lcs, bit_vector& seq_lcs, uint64_t& lcs,
  const align_parameters& parameters, bool print);

template<class ReferenceType>
void
alignBWTs(const ReferenceType& ref, const ReferenceType& seq,
  bit_vector& bwt_ref_lcs, bit_vector& bwt_seq_lcs,
  bit_vector& text_ref_lcs, bit_vector& text_seq_lcs,
  uint64_t& lcs,
  int_vector<0>& sa_samples, int_vector<0>& isa_samples,
  const align_parameters& parameters, bool print);

//------------------------------------------------------------------------------

const std::string RELATIVE_FM_EXTENSION = ".rfm";

/*
  The two indexes must have identical or sorted alphabets.
*/
template<class ReferenceBWTType = bwt_type, class SequenceType = bwt_type>
class RelativeFM
{
public:
  typedef SimpleFM<ReferenceBWTType> reference_type;

//------------------------------------------------------------------------------

  RelativeFM(const reference_type& ref, const reference_type& seq,
    const align_parameters parameters = align_parameters(), bool print = false) :
    reference(ref)
  {
    this->m_size = 0;
    this->sa_sample_rate = parameters.sa_sample_rate;
    this->isa_sample_rate = parameters.isa_sample_rate;

    uint64_t lcs_length = 0;
    bit_vector bwt_ref_lcs, bwt_seq_lcs, text_ref_lcs, text_seq_lcs;
    if(parameters.invariant)
    {
      if(!(ref.supportsLocate(true)) || !(seq.supportsLocate(true)))
      {
        std::cerr << "RelativeFM::RelativeFM(): Both indexes must have SA samples with this construction option!" << std::endl;
        return;
      }
      alignBWTs(ref, seq, bwt_ref_lcs, bwt_seq_lcs, text_ref_lcs, text_seq_lcs, lcs_length,
        this->sa_samples, this->isa_samples, parameters, print);
    }
    else
    {
      alignBWTs(ref, seq, bwt_ref_lcs, bwt_seq_lcs, lcs_length, parameters, print);
      this->buildSamples(seq);
    }

    getComplement(ref.bwt, this->ref_minus_lcs, bwt_ref_lcs, lcs_length);
    getComplement(seq.bwt, this->seq_minus_lcs, bwt_seq_lcs, lcs_length);
    util::assign(this->bwt_lcs, LCS(bwt_ref_lcs, bwt_seq_lcs, lcs_length));
    util::clear(bwt_ref_lcs); util::clear(bwt_seq_lcs);
    if(parameters.invariant)
    {
      util::assign(this->text_lcs, LCS(text_ref_lcs, text_seq_lcs, lcs_length));
      util::clear(text_ref_lcs); util::clear(text_seq_lcs);
    }
    this->alpha = seq.alpha;
    this->m_size = seq.size();
  }

  RelativeFM(const reference_type& ref, const std::string& base_name) :
    reference(ref)
  {
    std::string filename = base_name + RELATIVE_FM_EXTENSION;
    std::ifstream input(filename.c_str(), std::ios_base::binary);
    if(!input)
    {
      std::cerr << "RelativeFM::RelativeFM(): Cannot open input file " << filename << std::endl;
      return;
    }
    this->loadFrom(input);
    input.close();
  }

  RelativeFM(const SimpleFM<>& ref, std::istream& input) :
    reference(ref)
  {
    this->loadFrom(input);
  }

  ~RelativeFM()
  {
  }

//------------------------------------------------------------------------------

  inline uint64_t size() const { return this->m_size; }
  inline uint64_t sequences() const { return this->alpha.C[1]; }

  uint64_t reportSize(bool print = false) const
  {
    uint64_t ref_bytes = size_in_bytes(this->ref_minus_lcs);
    uint64_t seq_bytes = size_in_bytes(this->seq_minus_lcs);
    uint64_t bwt_bytes = ref_bytes + seq_bytes;
    uint64_t bwt_lcs_bytes = size_in_bytes(this->bwt_lcs);
    uint64_t text_lcs_bytes = size_in_bytes(this->text_lcs);
    uint64_t sa_bytes = size_in_bytes(this->sa_samples);
    uint64_t isa_bytes = size_in_bytes(this->isa_samples);

  #ifdef REPORT_RUNS
    uint64_t ref_runs = 0, seq_runs = 0;
    uint64_t ref_gap0 = 0, ref_gap1 = 0, ref_run = 0, ref_delta = 0;
    uint64_t seq_gap0 = 0, seq_gap1 = 0, seq_run = 0, seq_delta = 0;
    if(print)
    {
      countRuns(this->bwt_lcs.ref, ref_runs, ref_gap0, ref_gap1, ref_run, ref_delta);
      countRuns(this->bwt_lcs.seq, seq_runs, seq_gap0, seq_gap1, seq_run, seq_delta);
    }
  #endif

    uint64_t bytes = bwt_bytes + bwt_lcs_bytes + text_lcs_bytes + size_in_bytes(this->alpha)
                   + sizeof(this->m_size)
                   + sizeof(this->sa_sample_rate) + sa_bytes
                   + sizeof(this->isa_sample_rate) + isa_bytes;

    if(print)
    {
  #ifdef VERBOSE_OUTPUT
      printSize("ref_minus_lcs", ref_bytes, this->size());
      printSize("seq_minus_lcs", seq_bytes, this->size());
      printSize("bwt_lcs", bwt_lcs_bytes, this->size());
      if(this->text_lcs.size() > 0) { printSize("text_lcs", text_lcs_bytes, this->size()); }
      if(this->sa_sample_rate > 0) { printSize("SA samples", sa_bytes, this->size()); }
      if(this->isa_sample_rate > 0) { printSize("ISA samples", isa_bytes, this->size()); }
  #endif
  #ifdef REPORT_RUNS
      std::cout << std::string(16, ' ') << "Ref: " << ref_runs << " runs (gap0 "
        << (inMegabytes(ref_gap0) / 8) << " MB, gap1 " << (inMegabytes(ref_gap1) / 8)
        << " MB, run " << (inMegabytes(ref_run) / 8) << " MB, delta "
        << inMegabytes(ref_delta) << " MB)" << std::endl;
      std::cout << std::string(16, ' ') << "Seq: " << seq_runs << " runs (gap0 "
        << (inMegabytes(seq_gap0) / 8) << " MB, gap1 " << (inMegabytes(seq_gap1) / 8)
        << " MB, run " << (inMegabytes(seq_run) / 8) << " MB, delta "
        << inMegabytes(seq_delta) << " MB)" << std::endl;
  #endif
      printSize("Relative FM", bytes, this->size());
      std::cout << std::endl;
    }

    return bytes;
  }

  void writeTo(const std::string& base_name) const
  {
    std::string filename = base_name + RELATIVE_FM_EXTENSION;
    std::ofstream output(filename.c_str(), std::ios_base::binary);
    if(!output)
    {
      std::cerr << "RelativeFM::writeTo(): Cannot open output file " << filename << std::endl;
      return;
    }
    this->writeTo(output);
    output.close();
  }

  void writeTo(std::ostream& output) const
  {
    this->ref_minus_lcs.serialize(output);
    this->seq_minus_lcs.serialize(output);
    this->bwt_lcs.serialize(output);
    this->text_lcs.serialize(output);
    this->alpha.serialize(output);
    write_member(this->sa_sample_rate, output);
    this->sa_samples.serialize(output);
    write_member(this->isa_sample_rate, output);
    this->isa_samples.serialize(output);
  }

//------------------------------------------------------------------------------

  template<class Iter> range_type find(Iter begin, Iter end) const
  {
    range_type res(0, this->size() - 1);
    while(begin != end)
    {
      --end;
      if(!hasChar(this->alpha, *end)) { return range_type(1, 0); }
      uint64_t begin = cumulative(this->alpha, *end);
      res.first = begin + this->rank(res.first, *end);
      res.second = begin + this->rank(res.second + 1, *end) - 1;
      if(length(res) == 0) { return range_type(1, 0); }
    }
    return res;
  }

  inline range_type LF(uint64_t i) const
  {
    uint64_t lcs_bits = this->bwt_lcs.seq_rank(i + 1);  // Up to position i in seq.
    uint64_t ref_pos = (lcs_bits > 0 ? this->bwt_lcs.ref_select(lcs_bits) : 0);
    uint8_t ref_c = 0, seq_c = 0;
    bool lcs_pos = this->bwt_lcs.seq[i];

    if(lcs_pos)
    {
      ref_c = this->reference.bwt[ref_pos];
      seq_c = this->alpha.char2comp[this->reference.alpha.comp2char[ref_c]];
    }
    else
    {
      seq_c = this->seq_minus_lcs[i - lcs_bits];
      ref_c = this->reference.alpha.char2comp[this->alpha.comp2char[seq_c]];
    }

    return this->LF(i, lcs_bits, ref_pos, ref_c, seq_c, lcs_pos);
  }

  bool supportsLocate(bool print = false) const
  {
    if(this->text_lcs.size() == 0 && this->sa_sample_rate == 0)
    {
      if(print)
      {
        std::cerr << "RelativeFM::supportsLocate(): The index does not contain text_lcs or SA samples." << std::endl;
      }
      return false;
    }
    return true;
  }

  // Call supportsLocate() first.
  uint64_t locate(uint64_t i) const
  {
    if(this->text_lcs.size() == 0) { return relative::locate(*this, i); }
    if(i >= this->size()) { return 0; }

    // Find an SA sample or an LCS position in BWT(seq).
    uint64_t steps = 0;
    while(this->bwt_lcs.seq[i] == 0)
    {
      if(this->sa_sample_rate > 0 && i % this->sa_sample_rate == 0)
      {
        return this->sa_samples[i / this->sa_sample_rate] + steps;
      }
      range_type temp = this->LF_complement(i);
      if(temp.second == 0) { return steps; }
      i = temp.first; steps++;
    }

    // Map the LCS position into the corresponding position in BWT(ref).
    i = this->bwt_lcs.ref_select(this->bwt_lcs.seq_rank(i) + 1);

    // Locate the position in ref and convert it into a position in seq.
    // Note that it is the previous position in the reference that is an LCS position.
    i = this->reference.locate(i);
    i = this->text_lcs.seq_select(this->text_lcs.ref_rank(i)) + 1;

    return i + steps;
  }

  bool supportsExtract(bool print = false) const
  {
    if(this->text_lcs.size() == 0)
    {
      if(this->isa_sample_rate == 0)
      {
        if(print)
        {
          std::cerr << "RelativeFM::supportsExtract(): The index does not contain text_lcs or ISA samples." << std::endl;
        }
        return false;
      }
    }
    else if(this->isa_sample_rate == 0 && !(this->reference.supportsExtract(false)))
    {
      if(print)
      {
        std::cerr << "RelativeFM::supportsExtract(): The index contains text_lcs, but the reference does not support extract." << std::endl;
      }
      return false;
    }
    return true;
  }

  // Call supportsExtract() first.
  inline std::string extract(range_type range) const
  {
    return relative::extract(*this, range);
  }

  // Call supportsExtract() first.
  inline std::string extract(uint64_t from, uint64_t to) const
  {
    return relative::extract(*this, range_type(from, to));
  }

  // Returns ISA[i]. Call supportsExtract() first.
  inline uint64_t inverse(uint64_t i) const
  {
    uint64_t bwt_pos = 0, text_pos = this->size() - 1;  // The end of seq.
    if(this->isa_sample_rate > 0) // Is there an ISA sample before the end?
    {
      text_pos = ((i + this->isa_sample_rate - 1) / this->isa_sample_rate) * this->isa_sample_rate;
      if(text_pos >= this->size()) { text_pos = this->size() - 1; }
      else { bwt_pos = this->isa_samples[text_pos / this->isa_sample_rate]; }
    }
    if(this->text_lcs.size() > 0) // Is there an LCS position before text_pos?
    {
      uint64_t lcs_pos = this->text_lcs.seq_rank(i);
      uint64_t next_pos = this->text_lcs.seq_select(lcs_pos + 1);
      if(next_pos + 1 < text_pos)
      {
        text_pos = next_pos + 1;  // +1 because the matching BWT positions are one step forward.
        uint64_t ref_pos = this->text_lcs.ref_select(lcs_pos + 1);
        bwt_pos = this->reference.inverse(ref_pos + 1);
        bwt_pos = this->bwt_lcs.seq_select(this->bwt_lcs.ref_rank(bwt_pos) + 1);
      }
    }

    // Move backwards from the starting position.
    while(text_pos > i)
    {
      bwt_pos = this->LF(bwt_pos).first; text_pos--;
    }

    return bwt_pos;
  }

//------------------------------------------------------------------------------

  const reference_type& reference;
  SequenceType          ref_minus_lcs, seq_minus_lcs;
  LCS                   bwt_lcs, text_lcs;
  Alphabet              alpha;
  uint64_t              m_size;

  uint64_t              sa_sample_rate, isa_sample_rate;
  int_vector<0>         sa_samples, isa_samples;

//------------------------------------------------------------------------------

private:
  void loadFrom(std::istream& input)
  {
    this->ref_minus_lcs.load(input);
    this->seq_minus_lcs.load(input);
    this->bwt_lcs.load(input);
    this->text_lcs.load(input);
    this->alpha.load(input);
    this->m_size = this->reference.size() + this->seq_minus_lcs.size() - this->ref_minus_lcs.size();
    read_member(this->sa_sample_rate, input);
    this->sa_samples.load(input);
    read_member(this->isa_sample_rate, input);
    this->isa_samples.load(input);
  }

  void buildSamples(const reference_type& seq)
  {
    bool sample_sa = (this->sa_sample_rate > 0), sample_isa = (this->isa_sample_rate > 0);

    if(this->sa_sample_rate == seq.sa_sample_rate)
    {
      this->sa_samples = seq.sa_samples;
      sample_sa = false;
    }
    if(this->isa_sample_rate == seq.isa_sample_rate)
    {
      this->isa_samples = seq.isa_samples;
      sample_isa = false;
    }
    if(sample_sa)
    {
      util::assign(this->sa_samples,
        int_vector<0>((this->size() + this->sa_sample_rate - 1) / this->sa_sample_rate, 0, bitlength(this->size() - 1)));
    }
    if(sample_isa)
    {
      util::assign(this->isa_samples,
        int_vector<0>((this->size() + this->isa_sample_rate - 1) / this->isa_sample_rate, 0, bitlength(this->size() - 1)));
    }

    if(!sample_sa && !sample_isa) { return; }
    uint64_t text_pos = this->size(), bwt_pos = 0;
    while(text_pos > 0)
    {
      text_pos--;
      if(sample_sa && bwt_pos % this->sa_sample_rate == 0)
      {
        this->sa_samples[bwt_pos / this->sa_sample_rate] = text_pos;
      }
      if(sample_isa && text_pos % this->isa_sample_rate == 0)
      {
        this->isa_samples[text_pos / this->isa_sample_rate] = bwt_pos;
      }
      bwt_pos = seq.LF(bwt_pos).first;
    }
  }

//------------------------------------------------------------------------------

  /*
    The straightforward implementation counts the number of c's up to position i.
    To get the proper result, we need to check whether position i is in LCS or its complement,
    and ignore the last position when doing rank in that sequence.
    Note that c is a real character.
  */
  inline uint64_t rank(uint64_t i, uint8_t c) const
  {
    uint64_t res = 0;
    uint8_t ref_c = this->reference.alpha.char2comp[c], seq_c = this->alpha.char2comp[c];

    uint64_t lcs_bits = this->bwt_lcs.seq_rank(i + 1); // Number of LCS bits up to i in seq.
    bool check_lcs = (this->bwt_lcs.seq[i] == 1); // Is position i in LCS.
    if(lcs_bits < i + 1) // There are i + 1 - lcs_bits non-LCS bits in seq.
    {
      res += this->seq_minus_lcs.rank(i + check_lcs - lcs_bits, seq_c);
    }
    if(lcs_bits > 0)
    {
      uint64_t ref_pos = this->bwt_lcs.ref_select(lcs_bits);  // Select is 1-based.
      res += this->reference.bwt.rank(ref_pos + 1 - check_lcs, ref_c);
      if(lcs_bits < ref_pos + 1) // At least one non-LCS bit in ref.
      {
        res -= this->ref_minus_lcs.rank(ref_pos + 1 - lcs_bits, ref_c);
      }
    }

    return res;
  }

  /*
    The actual implementation of LF. The first function is a special case for positions
    that are known to be in seq_minus_lcs, while the second one contains the core
    functionality.
  */
  inline range_type LF_complement(uint64_t i) const
  {
    uint64_t lcs_bits = this->bwt_lcs.seq_rank(i);  // Up to position i in seq.
    uint64_t ref_pos = (lcs_bits > 0 ? this->bwt_lcs.ref_select(lcs_bits) : 0);
    uint8_t seq_c = this->seq_minus_lcs[i - lcs_bits];
    uint8_t ref_c = this->reference.alpha.char2comp[this->alpha.comp2char[seq_c]];
    return this->LF(i, lcs_bits, ref_pos, ref_c, seq_c, false);
  }

  // This function contains the core functionality of all varieties of LF.
  inline range_type LF(uint64_t i, uint64_t lcs_bits, uint64_t ref_pos,
    uint8_t ref_c, uint8_t seq_c, bool lcs_pos) const
  {
    uint64_t res = this->alpha.C[seq_c];
    if(lcs_bits < i + lcs_pos)
    {
      res += this->seq_minus_lcs.rank(i + lcs_pos - lcs_bits, seq_c);
    }
    if(lcs_bits > 0)
    {
      res += this->reference.bwt.rank(ref_pos + 1 - lcs_pos, ref_c);
      if(lcs_bits < ref_pos + 1)  // At least one non-LCS position in ref.
      {
        res -= this->ref_minus_lcs.rank(ref_pos + 1 - lcs_bits, ref_c);
      }
    }
    return range_type(res, seq_c);
  }

//------------------------------------------------------------------------------

  RelativeFM();
  RelativeFM(const RelativeFM&);
  RelativeFM(RelativeFM&&);
  RelativeFM& operator=(const RelativeFM&);
  RelativeFM& operator==(RelativeFM&);
};

//------------------------------------------------------------------------------

template<class ReferenceBWTType, class SequenceType>
void
getComplement(const ReferenceBWTType& bwt, SequenceType& output, const bit_vector& positions, uint64_t n)
{
#ifdef VERBOSE_STATUS_INFO
  double timestamp = readTimer();
#endif
  int_vector<8> buffer(bwt.size() - n, 0);
  for(uint64_t i = 0, j = 0; i < bwt.size(); i++)
  {
    if(positions[i] == 0) { buffer[j] = bwt[i]; j++; }
  }
  directConstruct(output, buffer);
#ifdef VERBOSE_STATUS_INFO
  #pragma omp critical (stderr)
  {
    std::cerr << "Extracted complement of length " << buffer.size() << " in "
              << (readTimer() - timestamp) << " seconds" << std::endl;
  }
#endif
}

//------------------------------------------------------------------------------

/*
  Ranges in ref and seq matching the pattern.
*/
struct record_type
{
  range_type  left, right;
  std::string pattern;
  bool        onlyNs, endmarker;

  record_type(range_type lt, range_type rt, const std::string& p, bool n, bool e) :
    left(lt), right(rt), pattern(p), onlyNs(n), endmarker(e)
  {
  }
};

// Sort by longest average length in descending order.
struct record_comparator
{
  inline bool operator() (const record_type& a, const record_type& b)
  {
    return (length(a.left) + length(a.right) > length(b.left) + length(b.right));
  }
};

inline std::ostream& operator<<(std::ostream& stream, const record_type& record)
{
  return stream << "(" << record.left << ", " << record.right << ")";
}

struct short_record_type
{
public:
  short_record_type(range_type lt, range_type rt, bool n, bool e) :
    left(lt), right(rt)
  {
    this->right.first = (this->right.first << 1) | n;
    this->right.second = (this->right.second << 1) | e;
  }

  short_record_type(const record_type& source) :
    left(source.left), right(source.right)
  {
    this->right.first = (this->right.first << 1) | source.onlyNs;
    this->right.second = (this->right.second << 1) | source.endmarker;
  }

  inline range_type first() const { return this->left; }
  inline range_type second() const
  {
    return range_type(this->right.first >> 1, this->right.second >> 1);
  }
  inline bool onlyNs() const { return this->right.first & 1; }
  inline bool endmarker() const { return this->right.second & 1; }

private:
  range_type left, right;
};

struct short_record_comparator
{
  inline bool operator() (const short_record_type& a, const short_record_type& b)
  {
    return (a.first() < b.first());
  }
};

/*
  Eugene W. Myers: An O(ND) Difference Algorithm and Its Variations. Algorithmica, 1986.

  The implementation assumes that offsets fit into int. If OpenMP is used, bitvectors
  should have the same length as the corresponding ranges, plus some padding to make the
  word boundaries match with the original bitvectors. Without OpenMP, they are the
  global bitvectors, and the function updates the range specified by ref_range/seq_range.
  The return value is (LCS length, heuristic losses).

  max_d specifies the maximum diagonal in LCS computation. If further diagonals would be
  needed, only the most frequent character in the ranges will be matched. Maximal memory
  usage will be around 4 * MAX_D * MAX_D bytes (+ SimpleFM for the reference and the
  target sequence and the current intervals as std::vector<uint8_t>).

  To use a preallocated buffer, set parameters.preallocated and pass a pointer to a buffer
  of size at least (max_d + 3) * (max_d + 1).
*/
range_type greedyLCS(const std::vector<uint8_t>& ref_extract, const std::vector<uint8_t>& seq_extract,
  bit_vector& ref_lcs, bit_vector& seq_lcs,
  short_record_type& record,
  const uint8_t* alphabet, uint64_t sigma,
  const align_parameters& parameters,
  std::vector<int>* buf);

//------------------------------------------------------------------------------

template<class ReferenceType>
void
processSubtree(const record_type& root, uint8_t* alphabet, uint64_t sigma,
  const ReferenceType& ref, const ReferenceType& seq,
  std::vector<short_record_type>& results,
  const align_parameters& parameters)
{
  std::stack<record_type> record_stack; record_stack.push(root);
  while(!(record_stack.empty()))
  {
    record_type curr = record_stack.top(); record_stack.pop();
    if(length(curr.left) <= parameters.block_size || length(curr.right) <= parameters.block_size || curr.endmarker)
    {
      results.push_back(short_record_type(curr));
      continue;
    }
    if(curr.pattern.length() >= parameters.max_length)
    {
#ifdef VERBOSE_STATUS_INFO
      #pragma omp critical (stderr)
      {
        std::cerr << "Pattern length became " << curr.pattern.length() << " on ranges "
                  << std::make_pair(curr.left, curr.right) << std::endl;
        std::cerr << "  Range lengths: " << std::make_pair(length(curr.left), length(curr.right)) << std::endl;
        std::cerr << "  The pattern was: " << curr.pattern << std::endl;
      }
#endif
      results.push_back(short_record_type(curr));
      continue;
    }

    std::string pattern = curr.pattern + " ";
    for(uint64_t i = sigma; i > 0; i--)
    {
      pattern[pattern.length() - 1] = (unsigned char)(alphabet[i - 1]);
      range_type left = ref.find(pattern.begin(), pattern.end()); if(isEmpty(left)) { continue; }
      range_type right = seq.find(pattern.begin(), pattern.end()); if(isEmpty(right)) { continue; }
      bool onlyNs = (curr.onlyNs && alphabet[i - 1] == 'N'), endmarker = (i == 1);
      record_stack.push(record_type(left, right, pattern, onlyNs, endmarker));
    }
  }
}

template<class ReferenceType>
void
generatePatterns(std::vector<record_type>& record_array, uint8_t* alphabet, uint64_t sigma,
  const ReferenceType& ref, const ReferenceType& seq,
  uint64_t seed_length)
{
  std::stack<record_type> record_stack;
  record_stack.push(record_type(range_type(0, ref.size() - 1), range_type(0, seq.size() - 1), "", true, false));
  record_array.clear();
  while(!(record_stack.empty()))
  {
    record_type curr = record_stack.top(); record_stack.pop();
    if(curr.pattern.length() >= seed_length) { record_array.push_back(curr); continue; }
    std::string pattern = curr.pattern + " ";
    for(uint64_t i = sigma; i > 0; i--)
    {
      pattern[pattern.length() - 1] = (unsigned char)(alphabet[i - 1]);
      range_type left = ref.find(pattern.begin(), pattern.end()); if(isEmpty(left)) { continue; }
      range_type right = seq.find(pattern.begin(), pattern.end()); if(isEmpty(right)) { continue; }
      bool onlyNs = (curr.onlyNs && alphabet[i - 1] == 'N'), endmarker = (i == 1);
      record_stack.push(record_type(left, right, pattern, onlyNs, endmarker));
    }
  }

  // Load balancing works better when long ranges are processed first.
  record_comparator comp;
  sequentialSort(record_array.begin(), record_array.end(), comp);
}

/*
  Copy the last length(range) bits from source to target[range]. We assume that source
  has been paddes so that word boundaries are in the same positions in both bitvectors.
*/
void copyBits(const bit_vector& source, bit_vector& target, range_type range);

template<class ReferenceType>
void
alignBWTs(const ReferenceType& ref, const ReferenceType& seq,
  bit_vector& ref_lcs, bit_vector& seq_lcs, uint64_t& lcs,
  const align_parameters& parameters, bool print)
{
  if(print)
  {
    std::cout << "Reference size: " << ref.size() << std::endl;
    std::cout << "Target size: " << seq.size() << std::endl;
  }
  util::clear(ref_lcs); util::clear(seq_lcs);

  uint64_t sigma = 0;
  uint8_t alphabet[256];
  if(parameters.sorted_alphabet)  // Use the intersection of the alphabets.
  {
    for(uint64_t c = 0; c < 256; c++)
    {
      if(hasChar(ref.alpha, c) && hasChar(seq.alpha, c))
      {
        alphabet[sigma] = c; sigma++;
      }
    }
  }
  else  // Use the alphabet from seq.
  {
    sigma = seq.alpha.sigma;
    for(uint64_t c = 0; c < sigma; c++) { alphabet[c] = seq.alpha.comp2char[c]; }
  }

  // Partition the BWTs.
  double timestamp = readTimer();
  uint64_t intersection = 0;
  std::vector<short_record_type> results;
  {
    std::vector<record_type> record_array;
    generatePatterns(record_array, alphabet, sigma, ref, seq, parameters.seed_length);
    #pragma omp parallel for schedule(dynamic, 1)
    for(uint64_t i = 0; i < record_array.size(); i++)
    {
      std::vector<short_record_type> result_buffer;
      processSubtree(record_array[i], alphabet, sigma, ref, seq, result_buffer, parameters);
      #pragma omp critical (alignbwts)
      {
        results.insert(results.end(), result_buffer.begin(), result_buffer.end());
      }
    }
    short_record_comparator comp;
    parallelMergeSort(results.begin(), results.end(), comp);
  }
  if(print)
  {
    uint64_t ref_missed = 0, seq_missed = 0, ref_expect = 0, seq_expect = 0, ref_overlap = 0, seq_overlap = 0;
    for(uint64_t i = 0; i < results.size(); i++)
    {
      range_type first = results[i].first(), second = results[i].second();
      intersection += std::min(length(first), length(second));
      if(first.first < ref_expect) { ref_overlap++; }
      ref_missed += first.first - ref_expect; ref_expect = first.second + 1;
      if(second.first < seq_expect) { seq_overlap++; }
      seq_missed += second.first - seq_expect; seq_expect = second.second + 1;
    }
    ref_missed += ref.size() - ref_expect; seq_missed += seq.size() - seq_expect;

    std::cout << "Found " << results.size() << " ranges with intersection length " << intersection
              << " in " << (readTimer() - timestamp) << " seconds" << std::endl;
    std::cout << "Partitioning misses: reference " << ref_missed
              << ", target " << seq_missed << std::endl;
    std::cout << "Partitioning losses: reference " << (ref.size() - intersection - ref_missed)
              << ", target " << (seq.size() - intersection - seq_missed) << std::endl;
    if(ref_overlap > 0 || seq_overlap > 0)
    {
      std::cout << "Overlapping ranges: reference " << ref_overlap << ", target " << seq_overlap << std::endl;
      std::cout << "(There is a bug somewhere)" << std::endl;
    }
  }

  // Find the approximate LCS using the partitioning.
  timestamp = readTimer(); lcs = 0;
  uint64_t losses = 0, processed = 0, percentage = 1;
  util::assign(ref_lcs, bit_vector(ref.size(), 0));
  util::assign(seq_lcs, bit_vector(seq.size(), 0));
  uint64_t threads = omp_get_max_threads();
  uint64_t chunk = std::max((uint64_t)1, results.size() / (threads * threads));
  std::vector<std::vector<int>*> buffers(threads, 0);
  if(parameters.preallocate)
  {
    for(uint64_t i = 0; i < threads; i++)
    {
      buffers[i] = new std::vector<int>;
      buffers[i]->reserve((parameters.max_d + 3) * (parameters.max_d + 1));
    }
  }
  #pragma omp parallel for schedule(dynamic, chunk)
  for(uint64_t i = 0; i < results.size(); i++)
  {
    // Ensure that ref_buffer and seq_buffer are aligned in the same way as ref_lcs and seq_lcs.
    range_type ref_range = results[i].first(), seq_range = results[i].second();
    bit_vector ref_buffer(length(ref_range) + ref_range.first % 64, 0);
    bit_vector seq_buffer(length(seq_range) + seq_range.first % 64, 0);
    std::vector<uint8_t> ref_extract; ref.extractBWT(ref_range, ref_extract);
    std::vector<uint8_t> seq_extract; seq.extractBWT(seq_range, seq_extract);

    range_type temp = greedyLCS(ref_extract, seq_extract, ref_buffer, seq_buffer, results[i],
      alphabet, sigma, parameters, buffers[omp_get_thread_num()]);
    #pragma omp critical (alignbwts)
    {
      copyBits(ref_buffer, ref_lcs, ref_range);
      copyBits(seq_buffer, seq_lcs, seq_range);
      lcs += temp.first; losses += temp.second;
    }
    #pragma omp critical (stderr)
    {
      processed++;
      if(processed >= (percentage / 100.0) * results.size())
      {
        double curr = readTimer();
        std::cerr << "Processed " << processed << " / " << results.size() << " ranges in "
                  << (curr - timestamp) << " seconds" << std::endl;
        percentage++;
      }
    }
  }
  if(parameters.preallocate)
  {
    for(uint64_t i = 0; i < threads; i++) { delete buffers[i]; buffers[i] = 0; }
  }
  if(print)
  {
    std::cout << "Found a common subsequence of length " << lcs << " in "
              << (readTimer() - timestamp) << " seconds" << std::endl;
    std::cout << "LCS losses: exact " << (intersection - lcs - losses) << ", heuristics " << losses << std::endl;
    std::cout << std::endl;
  }
}

//------------------------------------------------------------------------------

struct range_pair
{
  range_type ref_range, seq_range;  // Semiopen ranges [begin, end).

  range_pair(range_type ref, range_type seq) : ref_range(ref), seq_range(seq) {}
};

/*
  We have implicitly two arrays over the reference, where range ref_range has values seq_range.
  This function finds the longest increasing subsequence over the reference, where each value can be
  from either of the two arrays. Return value is subsequence length, the positions of the
  subsequence are marked in ref_lcs, and the values in the subsequence are marked in seq_lcs.

  The function clears left_matches and right_matches.
*/
uint64_t increasingSubsequence(std::vector<range_pair>& left_matches, std::vector<range_pair>& right_matches,
  bit_vector& ref_lcs, bit_vector& seq_lcs, uint64_t ref_len, uint64_t seq_len);

// Sort the matches and add a guard to the end.
void sortMatches(std::vector<range_pair>& matches, uint64_t ref_len, uint64_t seq_len);

// Determine the number of positions covered by a left/right match.
uint64_t countCoverage(std::vector<range_pair>& left_matches, std::vector<range_pair>& right_matches, uint64_t ref_len);

template<class ReferenceType>
struct Match
{
  const ReferenceType& seq;
  uint64_t bwt_pos;
  uint64_t text_pos;  // Ending position of the run in the text.
  uint64_t distance;  // Distance from SA[bwt_pos] to text_pos.
  uint64_t length;    // Length of the run.
  bool following;

  Match(const ReferenceType& _seq) :
    seq(_seq)
  {
    this->reset();
  }

  void reset()
  {
    this->bwt_pos = this->seq.size(); this->text_pos = this->seq.size();
    this->distance = 0; this->length = 0;
    this->following = false;
  }

  void follow(uint64_t new_pos)
  {
    this->bwt_pos = new_pos; this->text_pos = this->seq.size();
    this->distance = 0; this->length = 1;
    this->following = true;
    this->setSample();
  }

  inline void setSample()
  {
    if(this->bwt_pos % this->seq.sa_sample_rate == 0)
    {
      this->text_pos = this->seq.sa_samples[this->bwt_pos / this->seq.sa_sample_rate] + this->distance;
    }
  }

  // c is the character that should match the one used for LF.
  void advance(uint64_t c)
  {
    if(this->bwt_pos >= this->seq.size()) { return; }

    range_type temp = this->seq.LF(this->bwt_pos);
    if(!hasChar(this->seq.alpha, c) || this->seq.alpha.char2comp[c] != temp.second)
    {
      this->following = false;
    }
    if(temp.second == 0)
    {
      this->text_pos = this->distance;
      this->following = false;
    }
    else
    {
      this->bwt_pos = temp.first; this->distance++;
      this->setSample();
    }
  }

  void addRun(std::vector<range_pair>& vec, uint64_t ref_starting_pos)
  {
    if(this->length > 1)
    {
      range_type ref_run(ref_starting_pos, ref_starting_pos + this->length - 1);
      this->findSample();
      range_type seq_run(this->text_pos + 1 - this->length, this->text_pos);
      vec.push_back(range_pair(ref_run, seq_run));
    }
    this->reset();
  }

  void findSample()
  {
    while(this->text_pos >= this->seq.size()) { this->advance(0); }
  }
};

template<class ReferenceType>
void
findMatches(const ReferenceType& ref, const ReferenceType& seq,
  std::vector<range_pair>& left_matches, std::vector<range_pair>& right_matches,
  const bit_vector& merge_vector)
{
  uint64_t bwt_pos = 0, text_pos = ref.size() - 1;  // SA(ref)[bwt_pos] = text_pos
  Match<ReferenceType> left(seq), right(seq);
  right.follow(0);  // The empty suffix of seq is the right match of the empty suffix of ref.
  bit_vector::select_0_type ref_select(&merge_vector);
  bit_vector::rank_1_type seq_rank(&merge_vector);

  while(text_pos > 0)
  {
    // Advance to the previous position.
    range_type prev = ref.LF(bwt_pos); bwt_pos = prev.first; text_pos--;
    prev.second = ref.alpha.comp2char[prev.second];
    uint64_t mutual_pos = ref_select(bwt_pos + 1);

    // We stop following, if the character is wrong, the match has reached the beginning of
    // the sequence, there is no longer a match, or the old match is no longer adjacent.
    if(left.following)
    {
      left.advance(prev.second);
      uint64_t new_match = seq.size();
      if(merge_vector[mutual_pos - 1] == 0) { left.following = false; }
      else if((new_match = seq_rank(mutual_pos - 1)) != left.bwt_pos) { left.following = false; }
      if(left.following) { left.length++; }
      else
      {
        left.addRun(left_matches, text_pos + 1);
        if(new_match < seq.size()) { left.follow(new_match); }
      }
    }
    else if(mutual_pos > 0 && merge_vector[mutual_pos - 1] == 1)
    {
      left.follow(seq_rank(mutual_pos - 1));
    }

    if(right.following)
    {
      right.advance(prev.second);
      uint64_t new_match = seq.size();
      if(mutual_pos + 1 < seq.size() && merge_vector[mutual_pos + 1] == 0) { right.following = false; }
      else if((new_match = seq_rank(mutual_pos + 1)) != right.bwt_pos) { right.following = false; }
      if(right.following) { right.length++; }
      else
      {
        right.addRun(right_matches, text_pos + 1);
        if(new_match < seq.size()) { right.follow(new_match); }
      }
    }
    else if(mutual_pos + 1 < seq.size() && merge_vector[mutual_pos + 1] == 1)
    {
      right.follow(seq_rank(mutual_pos + 1));
    }
  }

  left.addRun(left_matches, 0); right.addRun(right_matches, 0);
}

template<class ReferenceType>
void
permuteVector(const bit_vector& source, bit_vector& target, const ReferenceType& index,
  uint64_t sa_sample_rate, int_vector<0>& sa_samples,
  uint64_t isa_sample_rate, int_vector<0>& isa_samples)
{
  bool sample_sa = (sa_sample_rate > 0);
  bool sample_isa = (isa_sample_rate > 0);
  if(sample_sa)
  {
    util::assign(sa_samples,
      int_vector<0>((index.size() + sa_sample_rate - 1) / sa_sample_rate, 0, bitlength(index.size() - 1)));
  }
  if(sample_isa)
  {
    util::assign(isa_samples,
      int_vector<0>((index.size() + isa_sample_rate - 1) / isa_sample_rate, 0, bitlength(index.size() - 1)));
  }
  util::assign(target, bit_vector(index.size(), 0));

  uint64_t text_pos = index.size(), bwt_pos = 0;
  while(text_pos > 1)
  {
    text_pos--; // Invariant: SA[bwt_pos] == text_pos
    target[bwt_pos] = source[text_pos - 1];
    if(sample_sa && bwt_pos % sa_sample_rate == 0) { sa_samples[bwt_pos / sa_sample_rate] = text_pos; }
    if(sample_isa && text_pos % isa_sample_rate == 0) { isa_samples[text_pos / isa_sample_rate] = bwt_pos; }
    bwt_pos = index.LF(bwt_pos).first;
  }
  target[bwt_pos] = source[index.size() - 1];
  if(sample_isa) { isa_samples[0] = bwt_pos; }
}

template<class ReferenceType>
void
alignBWTs(const ReferenceType& ref, const ReferenceType& seq,
  bit_vector& bwt_ref_lcs, bit_vector& bwt_seq_lcs,
  bit_vector& text_ref_lcs, bit_vector& text_seq_lcs,
  uint64_t& lcs,
  int_vector<0>& sa_samples, int_vector<0>& isa_samples,
  const align_parameters& parameters, bool print)
{
  if(print)
  {
    std::cout << "Reference size: " << ref.size() << std::endl;
    std::cout << "Target size: " << seq.size() << std::endl;
  }
  util::clear(bwt_ref_lcs); util::clear(bwt_seq_lcs);
  util::clear(text_ref_lcs); util::clear(text_seq_lcs);

  // Build a mapping from the comp values of seq to the comp values of ref.
  int_vector<64> comp_mapping(seq.alpha.sigma);
  bit_vector in_ref(seq.alpha.sigma);
  if(parameters.sorted_alphabet)  // Use first ref_c >= seq_c.
  {
    uint64_t ref_comp = 0;
    for(uint64_t seq_comp = 0; seq_comp < seq.alpha.sigma; seq_comp++)
    {
      uint64_t seq_c = seq.alpha.comp2char[seq_comp];
      while(ref_comp < ref.alpha.sigma && ref.alpha.comp2char[ref_comp] < seq_c) { ref_comp++; }
      comp_mapping[seq_comp] = ref_comp;
      in_ref[seq_comp] = (ref_comp < ref.alpha.sigma && ref.alpha.comp2char[ref_comp] == seq_c);
    }
  }
  else  // Assume an identical alphabet.
  {
    for(uint64_t i = 0; i < seq.alpha.sigma; i++) { comp_mapping[i] = i; in_ref[i] = true; }
  }

  // Build the merging bitvector.
  double timestamp = readTimer();
  bit_vector merge_vector(ref.size() + seq.size());
  {
    uint64_t ref_pos = ref.sequences(), seq_pos = 0; // Number of suffixes smaller than the current one.
    while(true)
    {
      merge_vector[ref_pos + seq_pos] = 1;
      range_type temp = seq.LF(seq_pos);
      if(temp.second == 0) { break; }
      seq_pos = temp.first;
      if(in_ref[temp.second])
      {
        temp.second = comp_mapping[temp.second];
        ref_pos = ref.alpha.C[temp.second] + ref.bwt.rank(ref_pos, temp.second);
      }
      else
      {
        ref_pos = ref.alpha.C[comp_mapping[temp.second]];
      }
    }
  }
  if(print)
  {
    std::cout << "Built the merging bitvector in " << (readTimer() - timestamp) << " seconds" << std::endl;
  }

  // Build left match and right match arrays.
  std::vector<range_pair> left_matches, right_matches;
  timestamp = readTimer();
  {
    findMatches(ref, seq, left_matches, right_matches, merge_vector);
    util::clear(merge_vector);
    sortMatches(left_matches, ref.size(), seq.size());
    sortMatches(right_matches, ref.size(), seq.size());
    if(print)
    {
      uint64_t total_matches = countCoverage(left_matches, right_matches, ref.size());
      std::cout << "Matched " << total_matches << " positions in "
                << (readTimer() - timestamp) << " seconds" << std::endl;
    }
  }

  // Find an increasing subsequence in the arrays and mark it in text_lcs bitvectors.
  timestamp = readTimer();
  lcs = increasingSubsequence(left_matches, right_matches, text_ref_lcs, text_seq_lcs, ref.size(), seq.size());
  if(print)
  {
    std::cout << "Found a common subsequence of length " << lcs << " in "
              << (readTimer() - timestamp) << " seconds" << std::endl;
  }

  // Traverse both CSAs to build the bwt_lcs bitvectors.
  // Sample SA at the same time if necessary.
  timestamp = readTimer();
  {
    const ReferenceType* index[2] = { &ref, &seq };
    bit_vector* text_lcs[2] = { &text_ref_lcs, &text_seq_lcs };
    bit_vector* bwt_lcs[2] = { &bwt_ref_lcs, &bwt_seq_lcs };
    uint64_t sa_rate[2] = { 0, parameters.sa_sample_rate };
    uint64_t isa_rate[2] = { 0, parameters.isa_sample_rate };
    #pragma omp parallel for
    for(uint64_t i = 0; i < 2; i++)
    {
      permuteVector(*(text_lcs[i]), *(bwt_lcs[i]), *(index[i]),
        sa_rate[i], sa_samples, isa_rate[i], isa_samples);
    }
  }
  if(print)
  {
    std::cout << "Built the bwt_lcs bitvectors ";
    if(parameters.sa_sample_rate > 0 || parameters.isa_sample_rate > 0) { std::cout << "and samples "; }
    std::cout << "in " << (readTimer() - timestamp) << " seconds" << std::endl;
  }
}

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_RELATIVE_FM_H
