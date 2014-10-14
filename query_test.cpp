#include "relative_fm.h"
#include "rlz_fm.h"
#include "rlz_vector.h"


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
  SimpleFM<> seq(argv[2]);
  seq.reportSize(true);
  RelativeFM rel(ref, argv[2]);
  rel.reportSize(true);
  RLZFM rlz(ref, argv[2]);
  rlz.reportSize(true);

  SimpleFM<wt_huff<rlz_vector>> rlzv(argv[2]);
  bit_vector::rank_1_type b_r(&(ref.bwt.bv));
  bit_vector::select_1_type b_s1(&(ref.bwt.bv));
  bit_vector::select_0_type b_s0(&(ref.bwt.bv));
  {
    rlz_vector& temp = const_cast<rlz_vector&>(rlzv.bwt.bv);
    temp.compress(ref.bwt.bv, b_r, b_s1, b_s0);
  }
  rlzv.reportSize(true);
  std::cout << std::endl;

  std::cout << "Patterns: " << argv[3] << std::endl;
  std::cout << std::endl;
  std::vector<std::string> patterns;
  uint64_t chars = readRows(argv[3], patterns, true);
  std::cout << patterns.size() << " patterns of total length " << chars << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  {
    double start = readTimer();
    uint64_t found = 0, matches = 0;
    for(auto pattern : patterns)
    {
      range_type res = seq.find(pattern.rbegin(), pattern.rend());
      if(length(res) > 0) { found++; matches += length(res); }
    }
    double seconds = readTimer() - start;
    printTime("Simple FM", found, matches, chars, seconds);
  }

  {
    double start = readTimer();
    uint64_t found = 0, matches = 0;
    for(auto pattern : patterns)
    {
      range_type res = rel.find(pattern.rbegin(), pattern.rend());
      if(length(res) > 0) { found++; matches += length(res); }
    }
    double seconds = readTimer() - start;
    printTime("Relative FM", found, matches, chars, seconds);
  }

  {
    double start = readTimer();
    uint64_t found = 0, matches = 0;
    for(auto pattern : patterns)
    {
      range_type res = rlz.find(pattern.rbegin(), pattern.rend());
      if(length(res) > 0) { found++; matches += length(res); }
    }
    double seconds = readTimer() - start;
    printTime("RLZ FM", found, matches, chars, seconds);
  }

  {
    double start = readTimer();
    uint64_t found = 0, matches = 0;
    for(auto pattern : patterns)
    {
      range_type res = rlzv.find(pattern.rbegin(), pattern.rend());
      if(length(res) > 0) { found++; matches += length(res); }
    }
    double seconds = readTimer() - start;
    printTime("RLZ vector", found, matches, chars, seconds);
  }

  double memory = inMegabytes(memoryUsage());
  std::cout << std::endl;
  std::cout << "Memory usage: " << memory << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
}
