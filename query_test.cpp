#include "relative_fm.h"
#include "rlz_fm.h"
#include "rlz_vector.h"
#include "sequence.h"


template<class Index>
void
testIndex(std::string name, const Index& index, std::vector<std::string>& patterns, uint64_t chars);

int
main(int argc, char** argv)
{
  if(argc < 4)
  {
    std::cerr << "Usage: query_test ref seq patterns" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  std::cout << "Query test" << std::endl;
  std::cout << std::endl;

  std::cout << "Reference: " << argv[1] << std::endl;
  std::cout << std::endl;
  SimpleFM<> ref(argv[1]);
  ref.reportSize(true);
  std::cout << std::endl;

  std::cout << "Sequence: " << argv[2] << std::endl;
  std::cout << std::endl;
  SimpleFM<> plain(argv[2]);
  plain.reportSize(true);
  SimpleFM<wt_huff<rrr_vector<63> > > rrr(argv[2]);
  rrr.reportSize(true);
  SimpleFM<Sequence> seq(argv[2]);
  seq.reportSize(true);
  SimpleFM<RLSequence> rle(argv[2]);
  rle.reportSize(true);
  RelativeFM rel(ref, argv[2]);
  rel.reportSize(true);

#ifdef TEST_RLZ_INDEXES
  RLZFM rlz(ref, argv[2]);
  rlz.reportSize(true);

  SimpleFM<wt_huff<rlz_vector> > rlzv(argv[2]);
  bit_vector::rank_1_type b_r(&(ref.bwt.bv));
  bit_vector::select_1_type b_s1(&(ref.bwt.bv));
  bit_vector::select_0_type b_s0(&(ref.bwt.bv));
  {
    rlz_vector& temp = const_cast<rlz_vector&>(rlzv.bwt.bv);
    temp.compress(ref.bwt.bv, b_r, b_s1, b_s0);
  }
  rlzv.reportSize(true);
#endif
  std::cout << std::endl;

  std::cout << "Patterns: " << argv[3] << std::endl;
  std::cout << std::endl;
  std::vector<std::string> patterns;
  uint64_t chars = readRows(argv[3], patterns, true);
  std::cout << patterns.size() << " patterns of total length " << chars << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  testIndex("FM<plain>", plain, patterns, chars);
  testIndex("FM<rrr>", rrr, patterns, chars);
  testIndex("FM<seq>", seq, patterns, chars);
  testIndex("FM<rle>", rle, patterns, chars);
  testIndex("Relative FM", rel, patterns, chars);

#ifdef TEST_RLZ_INDEXES
  testIndex("RLZ FM", rlz, patterns, chars);
  testIndex("FM<rlz>", rlzv, patterns, chars);
#endif

  double memory = inMegabytes(memoryUsage());
  std::cout << std::endl;
  std::cout << "Memory usage: " << memory << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
}


template<class Index>
void
testIndex(std::string name, const Index& index, std::vector<std::string>& patterns, uint64_t chars)
{
  double start = readTimer();
  uint64_t found = 0, matches = 0;
  for(auto pattern : patterns)
  {
    range_type res = index.find(pattern.begin(), pattern.end());
    if(length(res) > 0) { found++; matches += length(res); }
  }
  double seconds = readTimer() - start;
  printTime(name, found, matches, chars, seconds);
}
