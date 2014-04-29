#include "relative_fm.h"


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
  SimpleFM ref(argv[1]);
  ref.reportSize(true);
  std::cout << std::endl;

  std::cout << "Sequence: " << argv[2] << std::endl;
  std::cout << std::endl;
  SimpleFM seq(argv[2]);
  seq.reportSize(true);
  RelativeFM rel(ref, argv[2]);
  rel.reportSize(true);
  std::cout << std::endl;

  std::cout << "Patters: " << argv[3] << std::endl;
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
    double seconds = readTimer();
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
    double seconds = readTimer();
    printTime("Relative FM", found, matches, chars, seconds);
  }

  double memory = inMegabytes(memoryUsage());
  std::cout << std::endl;
  std::cout << "Memory usage: " << memory << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
}
