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

#ifndef _RELATIVE_FM_RLZ_FM_H
#define _RELATIVE_FM_RLZ_FM_H

#include "rlz.h"
#include "simple_fm.h"

namespace relative
{

//------------------------------------------------------------------------------

/*
  A classic FM-index based on relative Lempel-Ziv parsing of the BWT.
  Neither BWT should not contain (real) character value 1.
*/
class RLZFM
{
public:
  const static std::string EXTENSION;       // .rlzfm

  const static uint64_t    BLOCK_SIZE = 4;  // How many phrases in one rank() block.

  /*
    The last parameter is an optional index for the reverse reference.
  */
  RLZFM(const SimpleFM<>& ref, const SimpleFM<>& seq, const csa_wt<>* csa = 0);
  RLZFM(const SimpleFM<>& ref, const std::string& base_name);
  RLZFM(const SimpleFM<>& ref, std::istream& input);
  ~RLZFM();

  uint64_t size() const { return this->blocks.sum(); }
  uint64_t sequences() const { return this->alpha.C[1]; }

  uint64_t reportSize(bool print = false) const;
  void writeTo(const std::string& base_name) const;
  void writeTo(std::ostream& output) const;

  // Do not use with character value 0.
  template<class Iter> range_type find(Iter begin, Iter end) const;

  inline bool supportsLocate(bool print = false) const
  {
    if(print)
    {
      std::cerr << "RLZFM::supportsLocate(): The index does not support locate()." << std::endl;
    }
    return false;
  }

  // Call supportsLocate() first.
  inline uint64_t locate(uint64_t i) const { return i; }

  inline bool supportsExtract(bool print = false) const
  {
    if(print)
    {
      std::cerr << "RLZFM::supportsExtract(): The index does not support extract()." << std::endl;
    }
    return false;
  }

  // Call supportsExtract first.
  inline std::string extract(range_type range) const
  {
    std::cerr << "RLZFM::extract(): Trying to extract range " << range << std::endl;
    return std::string();
  }

  inline std::string extract(uint64_t from, uint64_t to) const { return this->extract(range_type(from, to)); }

  const SimpleFM<>&       reference;

  Alphabet                alpha;
  relative_encoder        phrases;
  CumulativeArray         blocks;
  CumulativeArray*        block_rank;
  int_vector<8>           mismatches; // FIXME this could be packed

private:
  // Counts the number of occurrences of (real character) c in reference[ref_pos, ref_pos + phrase_length - 1].
  inline uint64_t countOf(uint64_t ref_pos, uint64_t phrase_length, uint8_t c) const
  {
    if(phrase_length == 0) { return 0; }
    c = this->reference.alpha.char2comp[c];
    return this->reference.bwt.rank(ref_pos + phrase_length, c) - this->reference.bwt.rank(ref_pos, c);
  }

  uint64_t rank(uint64_t i, uint8_t c) const; // c is real character.

  void loadFrom(std::istream& input);

  RLZFM();
  RLZFM(const RLZFM&);
  RLZFM(RLZFM&&);
  RLZFM& operator=(const RLZFM&);
  RLZFM& operator=(RLZFM&&);
};

template<class Iter>
range_type
RLZFM::find(Iter begin, Iter end) const
{
  range_type res(0, this->size() - 1);
  while(begin != end)
  {
    --end;
    uint64_t pos = cumulative(this->alpha, *end);
    res.first = pos + this->rank(res.first, *end);
    res.second = pos + this->rank(res.second + 1, *end) - 1;
    if(length(res) == 0) { return range_type(1, 0); }
  }
  return res;
}

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_RLZ_FM_H
