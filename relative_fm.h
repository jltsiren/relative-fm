#ifndef _RELATIVE_FM_RELATIVE_FM_H
#define _RELATIVE_FM_RELATIVE_FM_H

#include "simple_fm.h"

namespace relative
{

//------------------------------------------------------------------------------

template<class ReferenceBWTType, class SequenceType>
void
getComplement(const ReferenceBWTType& bwt, SequenceType& output, const bit_vector& positions, uint64_t n);

std::vector<std::pair<int, int> >
mostFrequentChar(std::vector<uint8_t>& ref_buffer, std::vector<uint8_t>& seq_buffer);

template<class ReferenceType>
std::vector<std::pair<int, int> >
greedyLCS(const ReferenceType& ref, const ReferenceType& seq, range_type ref_range, range_type seq_range, bool onlyNs);

template<class ReferenceType>
std::pair<bit_vector, bit_vector>
alignBWTs(const ReferenceType& ref, const ReferenceType& seq, uint64_t block_size, uint max_depth, uint64_t& lcs, bool sorted_alphabet, bool print);

//------------------------------------------------------------------------------

const std::string RELATIVE_FM_EXTENSION = ".rfm";

template<class ReferenceBWTType = bwt_type, class SequenceType = bwt_type>
class RelativeFM
{
public:
  const static uint64_t BLOCK_SIZE = 1024;  // Split the BWTs into blocks of this size or less.
  const static uint     MAX_DEPTH  = 32;    // Maximum length of a pattern used to split the BWTs.

  typedef SimpleFM<ReferenceBWTType> reference_type;

  /*
    Maximum diagonal in LCS computation. If further diagonals would be needed, only the most frequent
    character in the ranges will be matched. Maximal memory usage will be around
    4 * MAX_D * MAX_D bytes (+ SimpleFM for the reference and the target sequence).
  */
  const static int MAX_D = 50000;

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
    auto lcs_vecs = alignBWTs(ref, seq, BLOCK_SIZE, MAX_DEPTH, lcs_length, sorted_alphabet, print);
    getComplement(ref.bwt, this->ref_minus_lcs, lcs_vecs.first, lcs_length);
    getComplement(seq.bwt, this->seq_minus_lcs, lcs_vecs.second, lcs_length);

    this->ref_lcs = lcs_vecs.first;
    this->seq_lcs = lcs_vecs.second;
  #ifndef USE_HYBRID_BITVECTORS
    util::init_support(this->ref_select, &(this->ref_lcs));
  #endif
    util::init_support(this->seq_rank, &(this->seq_lcs));

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

  uint64_t size() const { return this->m_size; }
  uint64_t sequences() const { return this->alpha.C[1]; }

  uint64_t reportSize(bool print = false) const
  {
    uint64_t ref_bytes = size_in_bytes(this->ref_minus_lcs);
    uint64_t seq_bytes = size_in_bytes(this->seq_minus_lcs);
    uint64_t bwt_bytes = ref_bytes + seq_bytes;

    uint64_t reflcs_bytes = size_in_bytes(this->ref_lcs);
  #ifndef USE_HYBRID_BITVECTORS
    reflcs_bytes +=size_in_bytes(this->ref_select);
  #endif
    uint64_t seqlcs_bytes = size_in_bytes(this->seq_lcs) + size_in_bytes(this->seq_rank);
    uint64_t bitvector_bytes = reflcs_bytes + seqlcs_bytes;

  #ifdef REPORT_RUNS
    uint64_t ref_runs = 0, seq_runs = 0;
    uint64_t ref_gap0 = 0, ref_gap1 = 0, ref_run = 0, ref_delta = 0;
    uint64_t seq_gap0 = 0, seq_gap1 = 0, seq_run = 0, seq_delta = 0;
    if(print)
    {
      countRuns(this->ref_lcs, ref_runs, ref_gap0, ref_gap1, ref_run, ref_delta);
      countRuns(this->seq_lcs, seq_runs, seq_gap0, seq_gap1, seq_run, seq_delta);
    }
  #endif

    uint64_t bytes = bwt_bytes + bitvector_bytes + size_in_bytes(this->alpha);

    if(print)
    {
  #ifdef VERBOSE_OUTPUT
      printSize("ref_minus_lcs", ref_bytes, this->size());
      printSize("seq_minus_lcs", seq_bytes, this->size());
      printSize("ref_lcs", reflcs_bytes, this->size());
  #ifdef REPORT_RUNS
      std::cout << std::string(16, ' ') << ref_runs << " runs (gap0 "
        << (inMegabytes(ref_gap0) / 8) << " MB, gap1 " << (inMegabytes(ref_gap1) / 8)
        << " MB, run " << (inMegabytes(ref_run) / 8) << " MB, delta "
        << inMegabytes(ref_delta) << " MB)" << std::endl;
  #endif
      printSize("seq_lcs", seqlcs_bytes, this->size());
  #ifdef REPORT_RUNS
      std::cout << std::string(16, ' ') << seq_runs << " runs (gap0 "
        << (inMegabytes(seq_gap0) / 8) << " MB, gap1 " << (inMegabytes(seq_gap1) / 8)
        << " MB, run " << (inMegabytes(seq_run) / 8) << " MB, delta "
        << inMegabytes(seq_delta) << " MB)" << std::endl;
  #endif
  #else
      printSize("BWT", bwt_bytes, this->size());
      printSize("Bitvectors", bitvector_bytes, this->size());
  #ifdef REPORT_RUNS
      std::cout << std::string(16, ' ') << ref_runs << " runs (gap0 "
        << (inMegabytes(ref_gap0) / 8) << " MB, gap1 " << (inMegabytes(ref_gap1) / 8)
        << " MB, run " << (inMegabytes(ref_run) / 8) << " MB, delta "
        << inMegabytes(ref_delta) << " MB)" << std::endl;
      std::cout << std::string(16, ' ') << seq_runs << " runs (gap0 "
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
    this->ref_lcs.serialize(output);
    this->seq_lcs.serialize(output);
    this->alpha.serialize(output);
  }

//------------------------------------------------------------------------------

  template<class Iter> range_type find(Iter begin, Iter end) const
  {
    range_type res(0, this->size() - 1);
    while(begin != end)
    {
      --end;
      uint64_t begin = cumulative(this->alpha, *end);
      res.first = begin + this->rank(res.first, *end);
      res.second = begin + this->rank(res.second + 1, *end) - 1;
      if(length(res) == 0) { return range_type(1, 0); }
    }
    return res;
  }

//------------------------------------------------------------------------------

#ifdef USE_HYBRID_BITVECTORS
  typedef hyb_vector<> vector_type;
#else
  typedef rrr_vector<63> vector_type;
#endif

  const reference_type&       reference;
  SequenceType                ref_minus_lcs, seq_minus_lcs;
  vector_type                 ref_lcs, seq_lcs;
  Alphabet                    alpha;
  uint64_t                    m_size;

  vector_type::rank_1_type    seq_rank;

#ifdef USE_HYBRID_BITVECTORS
  uint64_t ref_select(uint64_t i) const
  {
    if(i == 0) { return 0; }

    // Find the last position, where rank < i.
    uint64_t low = 0, high = this->ref_lcs.size();
    vector_type::rank_1_type ref_rank(&(this->ref_lcs));
    while(low < high)
    {
      uint64_t mid = low + (high - low + 1) / 2;
      if(ref_rank(mid) >= i) { high = mid - 1; }
      else { low = mid; }
    }

    return low;
  }
#else
  vector_type::select_1_type  ref_select;
#endif

//------------------------------------------------------------------------------

private:
  void loadFrom(std::istream& input)
  {
    this->ref_minus_lcs.load(input);
    this->seq_minus_lcs.load(input);
    this->ref_lcs.load(input);
    this->seq_lcs.load(input);
  #ifndef USE_HYBRID_BITVECTORS
    util::init_support(this->ref_select, &(this->ref_lcs));
  #endif
    util::init_support(this->seq_rank, &(this->seq_lcs));
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

    uint64_t lcs_bits = this->seq_rank(i + 1); // Number of LCS bits up to i in seq.
    bool check_lcs = (this->seq_lcs[i] == 1); // Is position i in LCS.
    if(lcs_bits < i + 1) // There are i + 1 - lcs_bits non-LCS bits in seq.
    {
      res += this->seq_minus_lcs.rank(i + check_lcs - lcs_bits, seq_c);
    }
    if(lcs_bits > 0)
    {
      uint64_t ref_pos = this->ref_select(lcs_bits);  // Select is 1-based.
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

  The implementation assumes that offsets fit into int.

  FIXME Space optimizations have not been implemented yet.
*/
template<class ReferenceType>
std::vector<std::pair<int, int> >
greedyLCS(const ReferenceType& ref, const ReferenceType& seq, range_type ref_range, range_type seq_range, bool onlyNs)
{
  std::vector<std::pair<int, int> > res;
  if(isEmpty(ref_range) || isEmpty(seq_range)) { return res; }

  std::vector<uint8_t> ref_buffer; ref.extractBWT(ref_range, ref_buffer); int ref_len = length(ref_range);
  std::vector<uint8_t> seq_buffer; seq.extractBWT(seq_range, seq_buffer); int seq_len = length(seq_range);

  if(onlyNs || abs(ref_len - seq_len) > RelativeFM<>::MAX_D)
  {
#ifdef VERBOSE_STATUS_INFO
      std::cout << "Reverting to heuristic on ranges " << std::make_pair(ref_range, seq_range) << std::endl;
#endif
      return mostFrequentChar(ref_buffer, seq_buffer);
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
      std::cout << "MAX_D exceeded on ranges " << std::make_pair(ref_range, seq_range) << std::endl;
#endif
      return mostFrequentChar(ref_buffer, seq_buffer);
    }
  }
  {
    std::vector<int> temp; v.swap(temp);  // Delete the contents.
  }

  // Extract the LCS.
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
      res.push_back(std::make_pair(x - 1, x - 1 - k));
    }
    k = next_k;
    std::vector<int> temp; store[d].swap(temp); // Delete the contents.
  }

  std::reverse(res.begin(), res.end());
  return res;
}

//------------------------------------------------------------------------------

/*
  Ranges in ref and seq matching the pattern.
*/
struct record_type
{
  range_type  left, right;
  std::string pattern;

  record_type(range_type l, range_type r, std::string p) : left(l), right(r), pattern(p) {}

  inline bool onlyNs() const
  {
    for(auto c : this->pattern) { if(c != 'N') { return false; } }
    return true;
  }
};

void verifyRanges(std::vector<record_type>& ranges, uint64_t ref_len, uint64_t seq_len);

struct record_comparator
{
  inline bool operator()(const record_type& a, const record_type& b) const
  {
    return ((a.left < b.left) || (a.left == b.left && a.right < b.right));
  }
};

inline std::ostream& operator<<(std::ostream& stream, const record_type& record)
{
  return stream << "(" << record.left << ", " << record.right << ", " << record.pattern << ")";
}

template<class ReferenceType>
std::pair<bit_vector, bit_vector>
alignBWTs(const ReferenceType& ref, const ReferenceType& seq, uint64_t block_size, uint max_depth, uint64_t& lcs, bool sorted_alphabet, bool print)
{
  if(print)
  {
    std::cout << "Reference size: " << ref.bwt.size() << std::endl;
    std::cout << "Target size: " << seq.bwt.size() << std::endl;
  }

  // Build the union of the alphabets.
  uint sigma = 0;
  uint8_t alphabet[256];
  if(sorted_alphabet)
  {
    for(uint c = 0; c < 256; c++)
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
    for(uint c = 0; c < sigma; c++) { alphabet[c] = seq.alpha.comp2char[c]; }
  }

  // Partition the BWTs.
  // FIXME: Could be faster, but does it matter?
  std::stack<record_type> record_stack;
  record_stack.push(record_type(range_type(0, ref.bwt.size() - 1), range_type(0, seq.bwt.size() - 1), ""));
  std::vector<record_type> ranges;
  while(!(record_stack.empty()))
  {
    record_type curr = record_stack.top(); record_stack.pop();
    range_type expect(curr.left.first, curr.right.first);

    for(uint i = 0; i < sigma; i++)
    {
      std::stringstream ss; ss << curr.pattern << (unsigned char)(alphabet[i]);
      std::string pattern = ss.str();

      range_type left = ref.find(pattern.begin(), pattern.end());
      if(isEmpty(left)) { left = range_type(expect.first, expect.first - 1); }
      else { expect.first = left.second + 1; }

      range_type right = seq.find(pattern.begin(), pattern.end());
      if(isEmpty(right)) { right = range_type(expect.second, expect.second - 1); }
      else { expect.second = right.second + 1; }

      if(length(left) > block_size && length(right) > block_size && pattern.length() < max_depth)
      {
        record_stack.push(record_type(left, right, pattern));
      }
      else if(length(left) > 0 || length(right) > 0)
      {
        ranges.push_back(record_type(left, right, pattern));
      }
    }
  }

  // Sort the ranges and check that they cover the entire BWTs.
  // verifyRanges() fails when the alphabets are different, but it does not matter.
  record_comparator comparator;
  std::sort(ranges.begin(), ranges.end(), comparator);
#ifdef VERIFY_CONSTRUCTION
  verifyRanges(ranges, ref.bwt.size(), seq.bwt.size());
#endif
  if(print) { std::cout << "Number of ranges: " << ranges.size() << std::endl; }

  // Find the approximate LCS using the partitioning.
  lcs = 0;
  bit_vector ref_lcs(ref.bwt.size(), 0), seq_lcs(seq.bwt.size(), 0);
  for(auto curr : ranges)
  {
    auto block = greedyLCS(ref, seq, curr.left, curr.right, curr.onlyNs());
    lcs += block.size();
    for(auto xy : block)
    {
      ref_lcs[curr.left.first + xy.first] = 1;
      seq_lcs[curr.right.first + xy.second] = 1;
    }
  }
  if(print)
  {
    std::cout << "Length of approximate LCS: " << lcs << std::endl;
    std::cout << std::endl;
  }

  return std::make_pair(ref_lcs, seq_lcs);
}

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_RELATIVE_FM_H
