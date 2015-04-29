/*
  Copyright (c) 2015 Genome Research Ltd.
  Copyright (c) 2014 Jouni Siren

  Author: Jouni Siren <jouni.siren@iki.fi>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#include "rlz_fm.h"

using namespace relative;


int
main(int argc, char** argv)
{
  if(argc < 3)
  {
    std::cerr << "Usage: build_rlzfm ref seq1 [seq2 ...]" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  std::cout << "Relative Lempel-Ziv FM-index builder" << std::endl;
  std::cout << std::endl;
  std::cout << "Reference: " << argv[1] << std::endl;
  std::cout << std::endl;

  SimpleFM<> ref(argv[1]);
  ref.reportSize(true);
  std::cout << std::endl;
  csa_wt<> csa;
  {
    int_vector<8> buffer(ref.size()); ref.extractBWT(buffer);
    reverseIndex(buffer, csa);
  }

  for(int arg = 2; arg < argc; arg++)
  {
    std::cout << "Target: " << argv[arg] << std::endl;
    SimpleFM<> seq(argv[arg]);
    double start = readTimer();
    RLZFM rel(ref, seq, &csa);
    double seconds = readTimer() - start;
    std::cout << "Index built in " << seconds << " seconds" << std::endl;
    std::cout << std::endl;

    rel.writeTo(argv[arg]);
    seq.reportSize(true);
    rel.reportSize(true);
    std::cout << std::endl;
  }

  double memory = inMegabytes(memoryUsage());
  std::cout << "Memory usage: " << memory << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
}
