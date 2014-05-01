#include "relative_fm.h"

//------------------------------------------------------------------------------

const std::string SimpleFM::EXTENSION = ".bwt";

SimpleFM::SimpleFM(const std::string& base_name)
{
  std::string filename = base_name + EXTENSION;
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "SimpleFM::SimpleFM(): Cannot open input file " << filename << "!" << std::endl;
    return;
  }
  this->loadFrom(input);
  input.close();
}

SimpleFM::SimpleFM(std::ifstream& input)
{
  this->loadFrom(input);
}

SimpleFM::~SimpleFM()
{
}

uint64_t
SimpleFM::reportSize(bool print) const
{
  uint64_t bwt_bytes = size_in_bytes(this->bwt);
  uint64_t bytes = sizeof(*this) + bwt_bytes + size_in_bytes(this->alpha);

  if(print)
  {
    printSize("BWT", bwt_bytes, this->bwt.size());
    printSize("Simple FM", bytes, this->bwt.size());
    std::cout << std::endl;
  }

  return bytes;
}

void
SimpleFM::writeTo(const std::string& base_name) const
{
  std::string filename = base_name + EXTENSION;
  std::ofstream output(filename.c_str(), std::ios_base::binary);
  if(!output)
  {
    std::cerr << "SimpleFM::writeTo(): Cannot open output file " << filename << "!" << std::endl;
    return;
  }
  this->writeTo(output);
  output.close();
}

void
SimpleFM::writeTo(std::ofstream& output) const
{
  this->bwt.serialize(output);
  this->alpha.serialize(output);
}

void
SimpleFM::loadFrom(std::ifstream& input)
{
  this->bwt.load(input);
  this->alpha.load(input);
}

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

RelativeFM::RelativeFM(const SimpleFM& ref, const SimpleFM& seq, bool print) :
  reference(ref)
{
  this->size = seq.bwt.size();

  uint64_t lcs_length = 0;
  auto lcs_vecs = alignBWTs(ref, seq, BLOCK_SIZE, MAX_DEPTH, lcs_length, print);
  this->ref_minus_lcs = getComplement(ref.bwt, lcs_vecs.first, lcs_length);
  this->seq_minus_lcs = getComplement(seq.bwt, lcs_vecs.second, lcs_length);

  this->ref_lcs = lcs_vecs.first;
  this->seq_lcs = lcs_vecs.second;
  this->buildRankSelect();

  this->alpha = seq.alpha;
}

RelativeFM::RelativeFM(const SimpleFM& ref, const std::string& base_name) :
  reference(ref)
{
  std::string filename = base_name + EXTENSION;
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "RelativeFM::RelativeFM(): Cannot open input file " << filename << "!" << std::endl;
    return;
  }
  this->loadFrom(input);
  input.close();
}

RelativeFM::RelativeFM(const SimpleFM& ref, std::ifstream& input) :
  reference(ref)
{
  this->loadFrom(input);
}

RelativeFM::~RelativeFM()
{
  delete ref_rank; ref_rank = 0;
  delete ref_select; ref_select = 0;
  delete seq_rank; seq_rank = 0;
  delete seq_select; seq_select = 0;
}

uint64_t
RelativeFM::reportSize(bool print) const
{
  uint64_t bwt_bytes = size_in_bytes(this->ref_minus_lcs) + size_in_bytes(this->seq_minus_lcs);

  uint64_t bitvector_bytes = size_in_bytes(this->ref_lcs) + size_in_bytes(seq_lcs);
  if(this->ref_rank != 0) { bitvector_bytes += size_in_bytes(*(this->ref_rank)); }
  if(this->ref_select != 0) { bitvector_bytes += size_in_bytes(*(this->ref_select)); }
  if(this->seq_rank != 0) { bitvector_bytes += size_in_bytes(*(this->seq_rank)); }
  if(this->seq_select != 0) { bitvector_bytes += size_in_bytes(*(this->seq_select)); }

#ifdef REPORT_RUNS
  uint64_t ref_runs = 0, seq_runs = 0;
  if(print)
  {
    bool in_run = false;
    for(uint64_t i = 0; i < this->ref_lcs.size(); i++)
    {
      if(this->ref_lcs[i] == 1) { if(!in_run) { ref_runs++; } in_run = true; }
      else { in_run = false;}
    }

    in_run = false;
    for(uint64_t i = 0; i < this->seq_lcs.size(); i++)
    {
      if(this->seq_lcs[i] == 1) { if(!in_run) { seq_runs++; } in_run = true; }
      else { in_run = false;}
    }
  }
#endif

  uint64_t bytes = sizeof(*this) + bwt_bytes + bitvector_bytes + size_in_bytes(this->alpha);

  if(print)
  {
    printSize("BWT", bwt_bytes, this->size);
    printSize("Bitvectors", bitvector_bytes, this->size);
#ifdef REPORT_RUNS
    std::cout << std::string(16, ' ') << ref_runs << " runs in ref_lcs, " << seq_runs << " runs in seq_lcs" << std::endl;
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
    std::cerr << "RelativeFM::writeTo(): Cannot open output file " << filename << "!" << std::endl;
    return;
  }
  this->writeTo(output);
  output.close();
}

void
RelativeFM::writeTo(std::ofstream& output) const
{
  this->ref_minus_lcs.serialize(output);
  this->seq_minus_lcs.serialize(output);
  this->ref_lcs.serialize(output);
  this->seq_lcs.serialize(output);
  this->alpha.serialize(output);
}

void
RelativeFM::loadFrom(std::ifstream& input)
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
  this->ref_rank = new vector_type::rank_1_type(&(this->ref_lcs));
  this->ref_select = new vector_type::select_1_type(&(this->ref_lcs));
  this->seq_rank = new vector_type::rank_1_type(&(this->seq_lcs));
  this->seq_select = new vector_type::select_1_type(&(this->seq_lcs));
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

std::vector<uint8_t>
extract(const bwt_type& bwt, range_type range)
{
  std::vector<uint8_t> res;

  for(uint64_t i = 0; i < length(range); i++)
  {
    res.push_back(bwt[range.first + i]);
  }

  return res;
}

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
greedyLCS(const bwt_type& ref, const bwt_type& seq, range_type ref_range, range_type seq_range)
{
  std::vector<std::pair<int, int> > res;
  if(isEmpty(ref_range) || isEmpty(seq_range)) { return res; }

  std::vector<uint8_t> ref_buffer = extract(ref, ref_range); int ref_len = length(ref_range);
  std::vector<uint8_t> seq_buffer = extract(seq, seq_range); int seq_len = length(seq_range);

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
alignBWTs(const SimpleFM& ref, const SimpleFM& seq, uint64_t block_size, uint max_depth, uint64_t& lcs, bool print)
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

      range_type left = ref.find(pattern.rbegin(), pattern.rend());
      if(isEmpty(left)) { left = range_type(expect.first, expect.first - 1); }
      else { expect.first = left.second + 1; }

      range_type right = seq.find(pattern.rbegin(), pattern.rend());
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
    auto block = greedyLCS(ref.bwt, seq.bwt, curr.left, curr.right);
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
