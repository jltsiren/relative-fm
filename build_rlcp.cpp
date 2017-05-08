/*
  Copyright (c) 2017 Genome Research Ltd.

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

#include "new_relative_lcp.h"

using namespace relative;

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 3)
  {
    std::cerr << "Usage: build_rlcp reference target1 [target2 ...]" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  std::cout << "Building the relative LCP array" << std::endl;
  std::cout << std::endl;

  std::string ref_name = std::string(argv[1]) + LCP_EXTENSION;
  printHeader("Reference"); std::cout << ref_name << std::endl;
  SLArray reference;
  sdsl::load_from_file(reference, ref_name);
  printHeader("Length"); std::cout << reference.size() << std::endl;
  printSize("LCP array", sdsl::size_in_bytes(reference), reference.size());
  std::cout << std::endl;

  for(int i = 2; i < argc; i++)
  {
    std::string target_name = std::string(argv[i]) + LCP_EXTENSION;
    printHeader("Target"); std::cout << target_name << std::endl;
    SLArray target;
    sdsl::load_from_file(target, target_name);
    printHeader("Length"); std::cout << target.size() << std::endl;

    double start = readTimer();
    NewRelativeLCP rlcp(reference, target);
    double seconds = readTimer() - start;
    printHeader("Construction"); std::cout << seconds << " seconds" << std::endl;
    double phrase_length = rlcp.size() / (double)(rlcp.phrases());
    printHeader("Phrases"); std::cout << rlcp.phrases() << " (average length " << phrase_length << ")" << std::endl;
    rlcp.reportSize(true);
    rlcp.writeTo(argv[i]);

    std::cout << std::endl;
  }

  std::cout << "Memory usage: " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------
