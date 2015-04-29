/*
  Copyright (c) 2015 Genome Research Ltd.

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

#include <string>

#include "utils.h"

using namespace relative;

//------------------------------------------------------------------------------

const double INDEL_RATE = 0.1;
const double INDEL_EXTEND = 0.8;

inline double
probability(std::mt19937_64& rng)
{
  return rng() / (rng.max() + 1.0);
}

inline uint64_t
indelLength(std::mt19937_64& rng)
{
  uint64_t len = 1;
  double prob = probability(rng);
  while(prob >= (1.0 - INDEL_EXTEND))
  {
    prob = (prob - (1.0 - INDEL_EXTEND)) / INDEL_EXTEND;
    len++;
  }
  return len;
}

const std::string ALPHABET = "ACGT";

inline uint64_t
randomChar(std::mt19937_64& rng)
{
  return ALPHABET[rng() % ALPHABET.size()];
}

inline uint64_t
substitute(std::mt19937_64& rng, uint64_t old_char)
{
  uint64_t res = 0;
  do
  {
    res = randomChar(rng);
  }
  while(res == old_char);
  return res;
}

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 4)
  {
    std::cerr << "Usage: mutate source target rate [seed]" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  std::cout << "Sequence mutator" << std::endl;
  std::cout << std::endl;

  int_vector_buffer<8> source, target;
  uint64_t seed = 0xDEADBEEF;
  double rate = 0.0;
  {
    std::string source_name = argv[1];
    std::cout << "Source: " << source_name << std::endl;
    source = int_vector_buffer<8>(source_name, std::ios::in, MEGABYTE, 8, true);
    std::cout << "Size: " << source.size() << std::endl;

    std::string target_name = argv[2];
    std::cout << "Target: " << target_name << std::endl;
    target = int_vector_buffer<8>(target_name, std::ios::out, MEGABYTE, 8, true);

    rate = std::stod(argv[3]);
    std::cout << "Mutation rate: " << rate << std::endl;

    if(argc > 4) { seed = fnv1a_hash((uint64_t)std::stoul(argv[4]), seed); }
    std::cout << "Seed: " << seed << std::endl;

    std::cout << std::endl;
  }


  std::mt19937_64 rng(seed);
  uint64_t substitutions = 0, insertions = 0, insertion_total = 0, deletions = 0, deletion_total = 0;
  for(uint64_t i = 0; i < source.size(); i++)
  {
    bool at_n = (source[i] == 'N'), after_n = (i > 0 && source[i - 1] == 'N');
    double mutation = probability(rng);
    if(mutation < rate)
    {
      mutation /= rate;
      if(mutation < INDEL_RATE / 2) // Insertion
      {
        insertions++;
        uint64_t len = indelLength(rng);
        insertion_total += len;
        if(at_n && after_n) // Insert N's inside a runs of N's.
        {
          for(uint64_t j = 0; j < len; j++) { target.push_back('N'); }
        }
        else
        {
          for(uint64_t j = 0; j < len; j++) { target.push_back(randomChar(rng)); }
        }
        target.push_back(source[i]);
      }
      else if(mutation < INDEL_RATE) // Deletion
      {
        deletions++;
        uint64_t len = indelLength(rng);
        len = std::min(len, source.size() - i);
        deletion_total += len;
        i += len - 1;
      }
      else  // Substitution
      {
        if(at_n)  // N's cannot be substituted with other characters.
        {
          target.push_back(source[i]);
        }
        else
        {
          substitutions++;
          target.push_back(substitute(rng, source[i]));
        }
      }
    }
    else { target.push_back(source[i]); }
  }

  std::cout << "Target size: " << target.size() << std::endl;
  std::cout << "Substitutions: " << substitutions << std::endl;
  std::cout << "Insertions: " << insertions << ", total size " << insertion_total << std::endl;
  std::cout << "Deletions: " << deletions << ", total size " << deletion_total << std::endl;
  std::cout << std::endl;

  target.close();
  return 0;
}

//------------------------------------------------------------------------------
