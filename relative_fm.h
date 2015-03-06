#ifndef _RELATIVE_FM_RELATIVE_FM_H
#define _RELATIVE_FM_RELATIVE_FM_H

#include "simple_fm.h"

namespace relative
{

//------------------------------------------------------------------------------

struct align_parameters
{
  const static uint64_t BLOCK_SIZE  = 1024; // Split the BWTs into blocks of this size or less.
  const static uint64_t MAX_LENGTH  = 32;   // Maximum length of a pattern used to split the BWTs.
  const static uint64_t SEED_LENGTH = 4;    // Generate all patterns of this length before parallelizing.
  const static int      MAX_D       = 8192; // Number of diagonals to consider in Myers' algorithm.

  // Default size for on demand LCS buffer; enough for d = 100.
  const static uint64_t BUFFER_SIZE = 103 * 101;

  align_parameters() :
    block_size(BLOCK_SIZE), max_length(MAX_LENGTH), seed_length(SEED_LENGTH), max_d(MAX_D),
    sorted_alphabet(true), preallocate(false)
  {
  }

  uint64_t block_size, max_length, seed_length;
  int      max_d;
  bool     sorted_alphabet; // comp values are sorted by character values.
  bool     preallocate;     // Use preallocated arrays in Myers' algorithm.
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

//------------------------------------------------------------------------------

class LCS
{
public:
  typedef uint64_t size_type;

#ifdef USE_HYBRID_BITVECTORS
  typedef hyb_vector<> vector_type;
#else
  typedef rrr_vector<63> vector_type;
#endif

  LCS();
  LCS(const bit_vector& a, const bit_vector& b, size_type _lcs_size);
  LCS(const LCS& l);
  LCS(LCS&& l);
  ~LCS();

  void swap(LCS& l);
  LCS& operator=(const LCS& l);
  LCS& operator=(LCS&& l);

  uint64_t serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  inline uint64_t size() const { return this->lcs_size; }
  inline uint64_t ref_size() const { return this->ref.size(); }
  inline uint64_t seq_size() const { return this->seq.size(); }

  /*
    To find how many LCS bits are in ref before the 0-based position i, use ref_rank(i).
    To find the 1-based LCS bit i in seq, use seq_select(i).
  */

  vector_type                ref;
  vector_type::rank_1_type   ref_rank;
#ifdef USE_HYBRID_BITVECTORS
  inline uint64_t ref_select(uint64_t i) const { return this->select(this->ref, this->ref_rank, i); }
#else
  vector_type::select_1_type  ref_select;
#endif

  vector_type                seq;
  vector_type::rank_1_type   seq_rank;
#ifdef USE_HYBRID_BITVECTORS
  inline uint64_t seq_select(uint64_t i) const { return this->select(this->seq, this->seq_rank, i); }
#else
  vector_type::select_1_type  seq_select;
#endif

  size_type lcs_size;

private:
  void copy(const LCS& l);
  void set_vectors();

#ifdef USE_HYBRID_BITVECTORS
  uint64_t select(const vector_type& vec, const vector_type::rank_1_type& rank, uint64_t i) const;
#endif
};  // class LCS

//------------------------------------------------------------------------------

const std::string RELATIVE_FM_EXTENSION = ".rfm";

template<class ReferenceBWTType = bwt_type, class SequenceType = bwt_type>
class RelativeFM
{
public:
  typedef SimpleFM<ReferenceBWTType> reference_type;

//------------------------------------------------------------------------------

  /*
    If the alphabet is not sorted, only characters in seq.alpha will be considered for partitioning
    the BWTs.
  */
  RelativeFM(const reference_type& ref, const reference_type& seq,
    const align_parameters parameters = align_parameters(), bool print = false) :
    reference(ref)
  {
    this->m_size = seq.size();

    uint64_t lcs_length = 0;
    bit_vector ref_lcs, seq_lcs;
    alignBWTs(ref, seq, ref_lcs, seq_lcs, lcs_length, parameters, print);

    getComplement(ref.bwt, this->ref_minus_lcs, ref_lcs, lcs_length);
    getComplement(seq.bwt, this->seq_minus_lcs, seq_lcs, lcs_length);
    util::assign(this->bwt_lcs, LCS(ref_lcs, seq_lcs, lcs_length));

    this->alpha = seq.alpha;
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

    uint64_t bytes = bwt_bytes + bwt_lcs_bytes + text_lcs_bytes + size_in_bytes(this->alpha) + sizeof(this->m_size);

    if(print)
    {
  #ifdef VERBOSE_OUTPUT
      printSize("ref_minus_lcs", ref_bytes, this->size());
      printSize("seq_minus_lcs", seq_bytes, this->size());
      printSize("bwt_lcs", bwt_lcs_bytes, this->size());
      if(this->text_lcs.size() > 0) { printSize("text_lcs", text_lcs_bytes, this->size()); }
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

  bool supportsLocate(bool print = false) const
  {
    if(this->text_lcs.size() == 0)
    {
      if(print)
      {
        std::cerr << "RelativeFM::supportsLocate(): The index does not contain text_lcs." << std::endl;
      }
      return false;
    }
    return true;
  }

  // Call supportsLocate() first.
  uint64_t locate(uint64_t i) const { return i; }

//------------------------------------------------------------------------------

  const reference_type&       reference;
  SequenceType                ref_minus_lcs, seq_minus_lcs;
  LCS                         bwt_lcs, text_lcs;
  Alphabet                    alpha;
  uint64_t                    m_size;

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
  }

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
  #ifdef _OPENMP
  #pragma omp critical (stderr)
  #endif
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

// The interpretation of ref_lcs and seq_lcs depends on whether OpenMP is used.
range_type mostFrequentChar(std::vector<uint8_t>& ref_buffer, std::vector<uint8_t>& seq_buffer,
  bit_vector& ref_lcs, bit_vector& seq_lcs,
  range_type ref_range, range_type seq_range);

// No LCS, just greedily match the runs.
// The interpretation of ref_lcs and seq_lcs depends on whether OpenMP is used.
range_type matchRuns(std::vector<uint8_t>& ref_buffer, std::vector<uint8_t>& seq_buffer,
  bit_vector& ref_lcs, bit_vector& seq_lcs,
  range_type ref_range, range_type seq_range,
  const uint8_t* alphabet, uint64_t sigma);

// Diagonal i has 2i+3 elements: -(i+1) to (i+1).
inline int offsetFor(int pos, int d) { return (d + 3) * d + pos + 1; }

/*
  Eugene W. Myers: An O(ND) Difference Algorithm and Its Variations. Algorithmica, 1986.

  The implementation assumes that offsets fit into int. If OpenMP is used, bitvectors
  should have the same length as the corresponding ranges. Without OpenMP, they are the
  global bitvectors, and the function updates the range specified by ref_range/seq_range.
  The return value is (LCS length, heuristic losses).

  max_d specifies the maximum diagonal in LCS computation. If further diagonals would be
  needed, only the most frequent character in the ranges will be matched. Maximal memory
  usage will be around 4 * MAX_D * MAX_D bytes (+ SimpleFM for the reference and the
  target sequence and the current intervals as std::vector<uint8_t>).

  To use a preallocated buffer, set parameters.preallocated and pass a pointer to a buffer
  of size at least (max_d + 3) * (max_d + 1).
*/
template<class ReferenceType>
range_type
greedyLCS(const ReferenceType& ref, const ReferenceType& seq,
  bit_vector& ref_lcs, bit_vector& seq_lcs,
  short_record_type& record,
  const uint8_t* alphabet, uint64_t sigma,
  const align_parameters& parameters,
  std::vector<int>* buf)
{
  range_type ref_range = record.first(), seq_range = record.second();
  bool onlyNs = record.onlyNs(), endmarker = record.endmarker();
  if(isEmpty(ref_range) || isEmpty(seq_range)) { return range_type(0, 0); }

  std::vector<uint8_t> ref_buffer; ref.extractBWT(ref_range, ref_buffer); int ref_len = length(ref_range);
  std::vector<uint8_t> seq_buffer; seq.extractBWT(seq_range, seq_buffer); int seq_len = length(seq_range);

  if(endmarker)
  {
    return matchRuns(ref_buffer, seq_buffer, ref_lcs, seq_lcs, ref_range, seq_range, alphabet, sigma);
  }
  else if(onlyNs || abs(ref_len - seq_len) > parameters.max_d)
  {
    return mostFrequentChar(ref_buffer, seq_buffer, ref_lcs, seq_lcs, ref_range, seq_range);
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
      while(x < ref_len && y < seq_len && ref_buffer[x] == seq_buffer[y]) { x++; y++; }
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
#ifdef _OPENMP
      #pragma omp critical (stderr)
#endif
      std::cerr << "MAX_D exceeded on ranges " << std::make_pair(ref_range, seq_range) << std::endl;
#endif
      if(!(parameters.preallocate)) { delete buf; buf = 0; }
      return mostFrequentChar(ref_buffer, seq_buffer, ref_lcs, seq_lcs, ref_range, seq_range);
    }
  }

  // Extract the LCS.
  uint64_t lcs = 0;
  d--;  // The last diagonal.
#ifdef _OPENMP
  uint64_t ref_padding = ref_range.first % 64, seq_padding = seq_range.first % 64;
#else
  uint64_t ref_padding = ref_range.first, seq_padding = seq_range.first;
#endif
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
#ifdef _OPENMP
      #pragma omp critical (stderr)
#endif
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

  // Build the union of the alphabets.
  uint64_t sigma = 0;
  uint8_t alphabet[256];
  if(parameters.sorted_alphabet)
  {
    for(uint64_t c = 0; c < 256; c++)
    {
      if(hasChar(ref.alpha, c) || hasChar(seq.alpha, c))
      {
        alphabet[sigma] = c; sigma++;
      }
    }
  }
  else
  {
    sigma = seq.alpha.sigma;
    for(uint64_t c = 0; c < sigma; c++) { alphabet[c] = seq.alpha.comp2char[c]; }
  }

  // Partition the BWTs.
  double timestamp = readTimer();
  uint64_t intersection = 0;
  std::vector<short_record_type> results;
  {
#ifdef _OPENMP
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
#else
    record_type root(range_type(0, ref.size() - 1), range_type(0, seq.size() - 1), "", true, false);
    processSubtree(root, alphabet, sigma, ref, seq, results, parameters);
#endif
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
#ifdef _OPENMP
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
    range_type temp = greedyLCS(ref, seq, ref_buffer, seq_buffer, results[i],
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
#else
  std::vector<int>* buffer = 0;
  if(parameters.preallocate)
  {
    buffer = new std::vector<int>;
    buffer.reserve((parameters.max_d + 3) * (parameters.max_d + 1));
  }
  for(uint64_t i = 0; i < results.size(); i++)
  {
    range_type temp = greedyLCS(ref, seq, ref_buffer, seq_buffer, results[i], alphabet, sigma, parameters, buffer);
    lcs += temp.first; losses += temp.second;
    processed++;
    if(processed >= (percentage / 100.0) * results.size())
    {
      double curr = readTimer();
      std::cerr << "Processed " << processed << " / " << results.size() << " ranges in "
                << (curr - timestamp) << " seconds" << std::endl;
      percentage++;
    }
  }
  if(parameters.preallocate) { delete buffer; buffer = 0; }
#endif
  if(print)
  {
    std::cout << "Found a common subsequence of length " << lcs << " in "
              << (readTimer() - timestamp) << " seconds" << std::endl;
    std::cout << "LCS losses: exact " << (intersection - lcs - losses) << ", heuristics " << losses << std::endl;
    std::cout << std::endl;
  }
}

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_RELATIVE_FM_H
