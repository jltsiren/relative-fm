#ifndef _RELATIVE_FM_RELATIVE_FM_H
#define _RELATIVE_FM_RELATIVE_FM_H

#include "simple_fm.h"

namespace relative
{

//------------------------------------------------------------------------------

template<class ReferenceBWTType, class SequenceType>
void
getComplement(const ReferenceBWTType& bwt, SequenceType& output, const bit_vector& positions, uint64_t n);

// The interpretation of ref_lcs and seq_lcs depends on whether OpenMP is used.
uint64_t
mostFrequentChar(std::vector<uint8_t>& ref_buffer, std::vector<uint8_t>& seq_buffer,
  bit_vector& ref_lcs, bit_vector& seq_lcs,
  range_type ref_range, range_type seq_range);

template<class ReferenceType>
uint64_t
greedyLCS(const ReferenceType& ref, const ReferenceType& seq,
  bit_vector& ref_lcs, bit_vector& seq_lcs,
  range_type ref_range, range_type seq_range, bool onlyNs);

template<class ReferenceType>
void
alignBWTs(const ReferenceType& ref, const ReferenceType& seq,
  bit_vector& ref_lcs, bit_vector& seq_lcs,
  uint64_t block_size, uint64_t max_depth, uint64_t& lcs, bool sorted_alphabet, bool print);

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
  const static uint64_t BLOCK_SIZE  = 1024; // Split the BWTs into blocks of this size or less.
  const static uint64_t MAX_DEPTH   = 32;   // Maximum length of a pattern used to split the BWTs.

  typedef SimpleFM<ReferenceBWTType> reference_type;

  /*
    Maximum diagonal in LCS computation. If further diagonals would be needed, only the most frequent
    character in the ranges will be matched. Maximal memory usage will be around
    4 * MAX_D * MAX_D bytes (+ SimpleFM for the reference and the target sequence and the
    current intervals as std::vector<uint8_t>).
  */
  const static int MAX_D = 10000;

//------------------------------------------------------------------------------

  /*
    If the alphabet is not sorted, only characters in seq.alpha will be considered for partitioning
    the BWTs.
  */
  RelativeFM(const reference_type& ref, const reference_type& seq, bool sorted_alphabet = true, bool print = false):
    reference(ref)
  {
    this->m_size = seq.bwt.size();

    uint64_t lcs_length = 0;
    bit_vector ref_lcs, seq_lcs;
    alignBWTs(ref, seq, ref_lcs, seq_lcs, BLOCK_SIZE, MAX_DEPTH, lcs_length, sorted_alphabet, print);

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
    uint64_t lcs_bytes = size_in_bytes(this->bwt_lcs);

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

    uint64_t bytes = bwt_bytes + lcs_bytes + size_in_bytes(this->alpha);

    if(print)
    {
  #ifdef VERBOSE_OUTPUT
      printSize("ref_minus_lcs", ref_bytes, this->size());
      printSize("seq_minus_lcs", seq_bytes, this->size());
      printSize("bwt_lcs", lcs_bytes, this->size());
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
  #else
      printSize("BWT", bwt_bytes, this->size());
      printSize("LCS", bitvector_bytes, this->size());
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

//------------------------------------------------------------------------------

  const reference_type&       reference;
  SequenceType                ref_minus_lcs, seq_minus_lcs;
  LCS                         bwt_lcs;
  Alphabet                    alpha;
  uint64_t                    m_size;

//------------------------------------------------------------------------------

private:
  void loadFrom(std::istream& input)
  {
    this->ref_minus_lcs.load(input);
    this->seq_minus_lcs.load(input);
    this->bwt_lcs.load(input);
    this->alpha.load(input);
    this->m_size = this->reference.bwt.size() + this->seq_minus_lcs.size() - this->ref_minus_lcs.size();
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
  int_vector<8> buffer(bwt.size() - n, 0);
  for(uint64_t i = 0, j = 0; i < bwt.size(); i++)
  {
    if(positions[i] == 0) { buffer[j] = bwt[i]; j++; }
  }
  directConstruct(output, buffer);
}

//------------------------------------------------------------------------------

/*
  Eugene W. Myers: An O(ND) Difference Algorithm and Its Variations. Algorithmica, 1986.

  The implementation assumes that offsets fit into int. If OpenMP is used, bitvectors
  should have the same length as the corresponding ranges. Without OpenMP, they are the
  global bitvectors, and the function updates the range specified by ref_range/seq_range.

  FIXME Space optimizations have not been implemented yet.
*/
template<class ReferenceType>
uint64_t
greedyLCS(const ReferenceType& ref, const ReferenceType& seq,
  bit_vector& ref_lcs, bit_vector& seq_lcs,
  range_type ref_range, range_type seq_range,
  bool onlyNs)
{
  if(isEmpty(ref_range) || isEmpty(seq_range)) { return 0; }

  std::vector<uint8_t> ref_buffer; ref.extractBWT(ref_range, ref_buffer); int ref_len = length(ref_range);
  std::vector<uint8_t> seq_buffer; seq.extractBWT(seq_range, seq_buffer); int seq_len = length(seq_range);

  if(onlyNs || abs(ref_len - seq_len) > RelativeFM<>::MAX_D)
  {
    return mostFrequentChar(ref_buffer, seq_buffer, ref_lcs, seq_lcs, ref_range, seq_range);
  }

  // v[k] stores how many characters of ref have been processed on diagonal k.
  std::vector<std::vector<int> > store;
  std::vector<int> v(3, 0); // v[mapToUint(1)] = 0;
  bool found = false;
  for(int d = 0; !found && d <= RelativeFM<>::MAX_D; d++)
  {
    for(int k = -d; k <= d; k += 2)
    {
      int x = 0;
      if(k == -d || (k != d && v[mapToUint(k - 1)] < v[mapToUint(k + 1)]))
      {
        x = v[mapToUint(k + 1)];
      }
      else
      {
        x = v[mapToUint(k - 1)] + 1;
      }
      int y = x - k;
      while(x < ref_len && y < seq_len && ref_buffer[x] == seq_buffer[y]) { x++; y++; }
      v[mapToUint(k)] = x;
      if(x >= ref_len && y >= seq_len) { found = true; }
    }
    store.push_back(v);
    v.resize(v.size() + 2, 0);

    // Too much memory required, match just the most frequent characters.
    if(d >= RelativeFM<>::MAX_D && !found)
    {
#ifdef VERBOSE_STATUS_INFO
#ifdef _OPENMP
      #pragma omp critical (stderr)
#endif
      std::cerr << "MAX_D exceeded on ranges " << std::make_pair(ref_range, seq_range) << std::endl;
#endif
      return mostFrequentChar(ref_buffer, seq_buffer, ref_lcs, seq_lcs, ref_range, seq_range);
    }
  }
  {
    std::vector<int> temp; v.swap(temp);  // Delete the contents.
  }

  // Extract the LCS.
  uint64_t lcs = 0;
  for(int d = store.size() - 1, k = ref_len - seq_len; d >= 0; d--)
  {
    int x_lim = 0, next_k = 0;
    if(d == 0) { x_lim = 0; }
    else if(k == -d || (k != d && store[d - 1][mapToUint(k - 1)] < store[d - 1][mapToUint(k + 1)]))
    {
      x_lim = store[d - 1][mapToUint(k + 1)]; next_k = k + 1;
    }
    else { x_lim = store[d - 1][mapToUint(k - 1)] + 1; next_k = k - 1; }
    for(int x = store[d][mapToUint(k)]; x > x_lim; x--)
    {
#ifdef _OPENMP
      ref_lcs[x - 1] = 1;
      seq_lcs[x - 1 - k] = 1;
#else
      ref_lcs[ref_range.first + x - 1] = 1;
      seq_lcs[seq_range.first + x - 1 - k] = 1;
#endif
      lcs++;
    }
    k = next_k;
  }

  return lcs;
}

//------------------------------------------------------------------------------

const static uint64_t ALIGN_BWTS_SEED_LENGTH = 3; // Generate all patterns of this length.

/*
  Ranges in ref and seq matching the pattern.
*/
struct record_type
{
  range_type  left, right;
  std::string pattern;
  bool        onlyNs;

  record_type(range_type lt, range_type rt, const std::string& p, bool n) :
    left(lt), right(rt), pattern(p), onlyNs(n)
  {
  }
};

struct short_record_type
{
  range_type left, right; // right.second also stores onlyNs information.

  short_record_type(range_type lt, range_type rt, bool n) :
    left(lt), right(rt)
  {
    this->right.second = (this->right.second << 1) | n;
  }

  inline range_type first() const { return this->left; }
  inline range_type second() const
  {
    return range_type(this->right.first, this->right.second >> 1);
  }
  inline bool onlyNs() const { return this->right.second & 1; }
};

struct short_record_comparator
{
  inline bool operator() (const short_record_type& a, const short_record_type& b)
  {
    return (a.left < b.left);
  }
};

inline std::ostream& operator<<(std::ostream& stream, const record_type& record)
{
  return stream << "(" << record.left << ", " << record.right << ")";
}

template<class ReferenceType>
void
processSubtree(const record_type& root, uint8_t* alphabet, uint64_t sigma,
  const ReferenceType& ref, const ReferenceType& seq,
  uint64_t block_size, uint64_t max_depth,
  std::vector<short_record_type>& results)
{
  std::stack<record_type> record_stack; record_stack.push(root);
  while(!(record_stack.empty()))
  {
    record_type curr = record_stack.top(); record_stack.pop();
    std::string pattern = curr.pattern + " ";

    for(uint64_t i = sigma; i > 0; i--)
    {
      pattern[pattern.length() - 1] = (unsigned char)(alphabet[i - 1]);
      range_type left = ref.find(pattern.begin(), pattern.end()); if(isEmpty(left)) { continue; }
      range_type right = seq.find(pattern.begin(), pattern.end()); if(isEmpty(right)) { continue; }

      bool N_pattern = (curr.onlyNs && alphabet[i - 1] == 'N');
      if(length(left) > block_size && length(right) > block_size && pattern.length() < max_depth)
      {
        record_stack.push(record_type(left, right, pattern, N_pattern));
      }
      else
      {
#ifdef VERBOSE_STATUS_INFO
        if(pattern.length() >= max_depth)
        {
#ifdef _OPENMP
          #pragma omp critical (stderr)
#endif
          {
            std::cerr << "Pattern length became " << pattern.length() << " on ranges "
                      << std::make_pair(left, right) << std::endl;
            std::cerr << "  Range lengths: " << std::make_pair(length(left), length(right)) << std::endl;
            std::cerr << "  The pattern was: " << pattern << std::endl;
          }
        }
#endif
        results.push_back(short_record_type(left, right, N_pattern));
      }
    }
  }
}

template<class ReferenceType>
void
generatePatterns(std::vector<record_type>& record_array, uint8_t* alphabet, uint64_t sigma,
  const ReferenceType& ref, const ReferenceType& seq,
  uint64_t len)
{
  std::stack<record_type> record_stack;
  record_stack.push(record_type(range_type(0, ref.bwt.size() - 1), range_type(0, seq.bwt.size() - 1), "", true));
  while(!(record_stack.empty()))
  {
    record_type curr = record_stack.top(); record_stack.pop();
    std::string pattern = curr.pattern + " ";

    for(uint64_t i = sigma; i > 0; i--)
    {
      pattern[pattern.length() - 1] = (unsigned char)(alphabet[i - 1]);
      range_type left = ref.find(pattern.begin(), pattern.end()); if(isEmpty(left)) { continue; }
      range_type right = seq.find(pattern.begin(), pattern.end()); if(isEmpty(right)) { continue; }

      bool N_pattern = (curr.onlyNs && alphabet[i - 1] == 'N');
      if(pattern.length() < len)
      {
        record_stack.push(record_type(left, right, pattern, N_pattern));
      }
      else
      {
        record_array.push_back(record_type(left, right, pattern, N_pattern));
      }
    }
  }
}

template<class ReferenceType>
void
alignBWTs(const ReferenceType& ref, const ReferenceType& seq,
  bit_vector& ref_lcs, bit_vector& seq_lcs,
  uint64_t block_size, uint64_t max_depth, uint64_t& lcs, bool sorted_alphabet, bool print)
{
  if(print)
  {
    std::cout << "Reference size: " << ref.bwt.size() << std::endl;
    std::cout << "Target size: " << seq.bwt.size() << std::endl;
  }
  util::clear(ref_lcs); util::clear(seq_lcs);

  // Build the union of the alphabets.
  uint64_t sigma = 0;
  uint8_t alphabet[256];
  if(sorted_alphabet)
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
  std::vector<short_record_type> results;
  {
#ifdef _OPENMP
    std::vector<record_type> record_array;
    generatePatterns(record_array, alphabet, sigma, ref, seq, ALIGN_BWTS_SEED_LENGTH);
    #pragma omp parallel for schedule(dynamic, 1)
    for(uint64_t i = 0; i < record_array.size(); i++)
    {
      std::vector<short_record_type> result_buffer;
      processSubtree(record_array[i], alphabet, sigma, ref, seq, block_size, max_depth, result_buffer);
      #pragma omp critical (alignbwts)
      {
        results.insert(results.end(), result_buffer.begin(), result_buffer.end());
      }
    }
    short_record_comparator comp;
    parallelSort(results.begin(), results.end(), comp);
#else
    record_type root(range_type(0, ref.bwt.size() - 1), range_type(0, seq.bwt.size() - 1), "", true);
    processSubtree(root, alphabet, sigma, ref, seq, block_size, max_depth, results);
#endif
  }
  if(print)
  {
    std::cout << "Found " << results.size() << " ranges in " << (readTimer() - timestamp) << " seconds" << std::endl;
  }

  // Find the approximate LCS using the partitioning.
  timestamp = readTimer(); lcs = 0;
  util::assign(ref_lcs, bit_vector(ref.bwt.size(), 0));
  util::assign(seq_lcs, bit_vector(seq.bwt.size(), 0));
#ifdef _OPENMP
  #ifdef VERBOSE_STATUS_INFO
  uint64_t processed = 0, percentage = 1;
  double prev = readTimer();
  #endif
  #pragma omp parallel for schedule(dynamic)
  for(uint64_t i = 0; i < results.size(); i++)
  {
    range_type ref_range = results[i].first(), seq_range = results[i].second();
    bit_vector ref_buffer(length(ref_range), 0), seq_buffer(length(seq_range), 0);
    uint64_t cur_lcs = greedyLCS(ref, seq, ref_buffer, seq_buffer, ref_range, seq_range, results[i].onlyNs());
    #pragma omp critical (alignbwts)
    {
      for(uint64_t j = ref_range.first; j <= ref_range.second; j += 64)
      {
        uint64_t bits = std::min((uint64_t)64, ref_range.second + 1 - j);
        uint64_t temp = ref_buffer.get_int(j - ref_range.first, bits);
        ref_lcs.set_int(j, temp, bits);
      }
      for(uint64_t j = seq_range.first; j <= seq_range.second; j += 64)
      {
        uint64_t bits = std::min((uint64_t)64, seq_range.second + 1 - j);
        uint64_t temp = seq_buffer.get_int(j - seq_range.first, bits);
        seq_lcs.set_int(j, temp, bits);
      }
      lcs += cur_lcs;
    }
  #ifdef VERBOSE_STATUS_INFO
    #pragma omp critical (stderr)
    {
      processed++;
      if(processed >= (percentage / 100.0) * results.size())
      {
        double curr = readTimer();
        std::cerr << "Processed " << processed << " / " << results.size() << " ranges in "
                  << (curr - prev) << " seconds" << std::endl;
        percentage++; prev = curr;
      }
    }
  #endif
  }
#else
  for(uint64_t i = 0; i < results.size(); i++)
  {
    lcs += greedyLCS(ref, seq, ref_buffer, seq_buffer,
      results[i].first(), results[i].second(), results[i].onlyNs());
  }
#endif
  if(print)
  {
    std::cout << "Found a common subsequence of length " << lcs << " in "
              << (readTimer() - timestamp) << " seconds" << std::endl;
    std::cout << std::endl;
  }
}

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_RELATIVE_FM_H
