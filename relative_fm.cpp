#include <cstdlib>

#include "relative_fm.h"

//------------------------------------------------------------------------------

const std::string RelativeFM::EXTENSION = ".rfm";

bwt_type
getComplement(const bwt_type& bwt, const bit_vector& positions, uint64_t n)
{
  int_vector<8> buffer(bwt.size() - n, 0);
  for(uint64_t i = 0, j = 0; i < bwt.size(); i++)
  {
    if(positions[i] == 0) { buffer[j] = bwt[i]; j++; }
  }

  bwt_type temp;
  construct_im(temp, buffer);
  return temp;
}

RelativeFM::RelativeFM(const SimpleFM<>& ref, const SimpleFM<>& seq, bool print) :
  reference(ref)
{
  this->size = seq.bwt.size();

  uint64_t lcs_length = 0;
  auto lcs_vecs = alignBWTs(ref, seq, BLOCK_SIZE, MAX_DEPTH, lcs_length, print);
  this->ref_minus_lcs = getComplement(ref.bwt, lcs_vecs.first, lcs_length);
  this->seq_minus_lcs = getComplement(seq.bwt, lcs_vecs.second, lcs_length);

#ifdef USE_SPARSE_BITVECTORS
  lcs_vecs.first.flip();
  lcs_vecs.second.flip();
#endif
  this->ref_lcs = lcs_vecs.first;
  this->seq_lcs = lcs_vecs.second;
  this->buildRankSelect();

  this->alpha = seq.alpha;
}

RelativeFM::RelativeFM(const SimpleFM<>& ref, const std::string& base_name) :
  reference(ref)
{
  std::string filename = base_name + EXTENSION;
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "RelativeFM::RelativeFM(): Cannot open input file " << filename << std::endl;
    return;
  }
  this->loadFrom(input);
  input.close();
}

RelativeFM::RelativeFM(const SimpleFM<>& ref, std::istream& input) :
  reference(ref)
{
  this->loadFrom(input);
}

RelativeFM::~RelativeFM()
{
}

void
countRuns(const RelativeFM::vector_type& vec, uint64_t& runs, uint64_t& gap0, uint64_t& gap1, uint64_t& run, uint64_t& delta)
{
  bool in_run = false;
  uint64_t prev0 = 0, prev1 = 0, prevr = 0, onebits = 0;
  for(uint64_t i = 0; i < vec.size(); i++)
  {
    if(vec[i] == 1)
    {
      gap1 += bitlength(i + 1 - prev1);
      if(!in_run)
      {
        run += bitlength(i + 1 - prevr);
        runs++; prevr = i + 1;
      }
      in_run = true; prev1 = i + 1;
      onebits++;
    }
    else
    {
      gap0 += bitlength(i + 1 - prev0);
      if(in_run)
      {
        run += bitlength(i + 1 - prevr);
        prevr = i + 1;
      }
      in_run = false; prev0 = i + 1;
    }
  }

  int_vector<64> buffer(onebits);
  for(uint64_t i = 0, j = 0; i < vec.size(); i++)
  {
    if(vec[i] == 1) { buffer[j] = i; j++; }
  }
  enc_vector<> delta_vec(buffer);
  delta = size_in_bytes(delta_vec);
}

uint64_t
RelativeFM::reportSize(bool print) const
{
  uint64_t ref_bytes = size_in_bytes(this->ref_minus_lcs);
  uint64_t seq_bytes = size_in_bytes(this->seq_minus_lcs);
  uint64_t bwt_bytes = ref_bytes + seq_bytes;

  uint64_t reflcs_bytes = size_in_bytes(this->ref_lcs) + size_in_bytes(this->ref_select);
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
    printSize("ref_minus_lcs", ref_bytes, this->size);
    printSize("seq_minus_lcs", seq_bytes, this->size);
    printSize("ref_lcs", reflcs_bytes, this->size);
#ifdef REPORT_RUNS
    std::cout << std::string(16, ' ') << ref_runs << " runs (gap0 "
      << (inMegabytes(ref_gap0) / 8) << " MB, gap1 " << (inMegabytes(ref_gap1) / 8)
      << " MB, run " << (inMegabytes(ref_run) / 8) << " MB, delta "
      << inMegabytes(ref_delta) << " MB)" << std::endl;
#endif
    printSize("seq_lcs", seqlcs_bytes, this->size);
#ifdef REPORT_RUNS
    std::cout << std::string(16, ' ') << seq_runs << " runs (gap0 "
      << (inMegabytes(seq_gap0) / 8) << " MB, gap1 " << (inMegabytes(seq_gap1) / 8)
      << " MB, run " << (inMegabytes(seq_run) / 8) << " MB, delta "
      << inMegabytes(seq_delta) << " MB)" << std::endl;
#endif
#else
    printSize("BWT", bwt_bytes, this->size);
    printSize("Bitvectors", bitvector_bytes, this->size);
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
    printSize("Relative FM", bytes, this->size);
    std::cout << std::endl;
  }

  return bytes;
}

void
RelativeFM::writeTo(const std::string& base_name) const
{
  std::string filename = base_name + EXTENSION;
  std::ofstream output(filename.c_str(), std::ios_base::binary);
  if(!output)
  {
    std::cerr << "RelativeFM::writeTo(): Cannot open output file " << filename << std::endl;
    return;
  }
  this->writeTo(output);
  output.close();
}

void
RelativeFM::writeTo(std::ostream& output) const
{
  this->ref_minus_lcs.serialize(output);
  this->seq_minus_lcs.serialize(output);
  this->ref_lcs.serialize(output);
  this->seq_lcs.serialize(output);
  this->alpha.serialize(output);
}

void
RelativeFM::loadFrom(std::istream& input)
{
  this->ref_minus_lcs.load(input);
  this->seq_minus_lcs.load(input);
  this->ref_lcs.load(input);
  this->seq_lcs.load(input);
  this->buildRankSelect();
  this->alpha.load(input);
  this->size = this->reference.bwt.size() + this->seq_minus_lcs.size() - this->ref_minus_lcs.size();
}

void
RelativeFM::buildRankSelect()
{
  util::init_support(this->ref_select, &(this->ref_lcs));
  util::init_support(this->seq_rank, &(this->seq_lcs));
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

  bool onlyNs() const
  {
    for(auto c : this->pattern) { if(c != 'N') { return false; } }
    return true;
  }
};

struct
{
  bool operator()(const record_type& a, const record_type& b) const
  {
    return ((a.left < b.left) || (a.left == b.left && a.right < b.right));
  }
} record_comparator;

std::ostream& operator<<(std::ostream& stream, const record_type& record)
{
  return stream << "(" << record.left << ", " << record.right << ", " << record.pattern << ")";
}

void
verifyRanges(std::vector<record_type>& ranges, uint64_t ref_len, uint64_t seq_len)
{
  range_type expect(0, 0);
  for(auto curr : ranges)
  {
    if(curr.left.first != expect.first ||  curr.right.first != expect.second)
    {
      std::cerr << "verifyRanges(): Expected " << expect << ", got " << curr << "!" << std::endl;
    }
    expect = range_type(curr.left.second + 1, curr.right.second + 1);
  }
  if(expect.first != ref_len || expect.second != seq_len)
  {
    std::cerr << "verifyRanges(): Final expect was " << expect << "!" << std::endl;
  }
}

//------------------------------------------------------------------------------

uint
mapToUint(int val)
{
  if(val >= 0) { return 2 * val; }
  else         { return 2 * (-val) - 1; }
}

std::vector<std::pair<int, int> >
mostFrequentChar(std::vector<uint8_t>& ref_buffer, std::vector<uint8_t>& seq_buffer)
{
  std::vector<std::pair<int, int> > freqs(256, std::make_pair(0, 0));
  for(auto c : ref_buffer) { freqs[c].first++; }
  for(auto c : seq_buffer) { freqs[c].second++; }
  for(auto& f : freqs) { f.first = std::min(f.first, f.second); }
  auto max_pos = std::max_element(freqs.begin(), freqs.end());
  uint max_c = max_pos - freqs.begin(), max_f = max_pos->first;

  std::vector<std::pair<int, int> > res; res.reserve(max_f);
  auto ref_p = ref_buffer.begin(), seq_p = seq_buffer.begin();
  while(res.size() < max_f)
  {
    while(*ref_p != max_c) { ++ref_p; }
    while(*seq_p != max_c) { ++seq_p; }
    res.push_back(std::make_pair(ref_p - ref_buffer.begin(), seq_p - seq_buffer.begin()));
    ++ref_p; ++seq_p;
  }

#ifdef VERBOSE_OUTPUT
  std::cout << "  Matched " << max_f << " copies of character " << max_c << " (" << ((char)max_c) << ")";
  std::cout << " from sequences of length " << range_type(ref_buffer.size(), seq_buffer.size()) << std::endl;
#endif

  return res;
}

/*
  Eugene W. Myers: An O(ND) Difference Algorithm and Its Variations. Algorithmica, 1986.

  The implementation assumes that offsets fit into int.

  FIXME Space optimizations have not been implemented yet.
*/
std::vector<std::pair<int, int> >
greedyLCS(const SimpleFM<>& ref, const SimpleFM<>& seq, range_type ref_range, range_type seq_range, bool onlyNs)
{
  std::vector<std::pair<int, int> > res;
  if(isEmpty(ref_range) || isEmpty(seq_range)) { return res; }

  std::vector<uint8_t> ref_buffer; ref.extractBWT(ref_range, ref_buffer); int ref_len = length(ref_range);
  std::vector<uint8_t> seq_buffer; seq.extractBWT(seq_range, seq_buffer); int seq_len = length(seq_range);

  if(onlyNs || abs(ref_len - seq_len) > RelativeFM::MAX_D)
  {
#ifdef VERBOSE_OUTPUT
      std::cout << "Reverting to heuristic on ranges " << std::make_pair(ref_range, seq_range) << std::endl;
#endif
      return mostFrequentChar(ref_buffer, seq_buffer);
  }

  // v[k] stores how many characters of ref have been processed on diagonal k.
  std::vector<std::vector<int> > store;
  std::vector<int> v(3, 0); // v[mapToUint(1)] = 0;
  bool found = false;
  for(int d = 0; !found && d <= RelativeFM::MAX_D; d++)
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
    if(d >= RelativeFM::MAX_D && !found)
    {
#ifdef VERBOSE_OUTPUT
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

std::pair<bit_vector, bit_vector>
alignBWTs(const SimpleFM<>& ref, const SimpleFM<>& seq, uint64_t block_size, uint max_depth, uint64_t& lcs, bool print)
{
  if(print)
  {
    std::cout << "Reference size: " << ref.bwt.size() << std::endl;
    std::cout << "Target size: " << seq.bwt.size() << std::endl;
  }

  uint sigma = 0;
  uint8_t alphabet[256];
  for(uint c = 0; c < 256; c++)
  {
    if(hasChar(ref.alpha, c) || hasChar(seq.alpha, c))
    {
      alphabet[sigma] = c; sigma++;
    }
  }

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

  std::sort(ranges.begin(), ranges.end(), record_comparator);
  verifyRanges(ranges, ref.bwt.size(), seq.bwt.size());
  if(print) { std::cout << "Number of ranges: " << ranges.size() << std::endl; }

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
