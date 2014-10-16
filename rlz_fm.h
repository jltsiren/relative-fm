#ifndef _RELATIVE_FM_RLZ_FM_H
#define _RELATIVE_FM_RLZ_FM_H


#include "relative_fm.h"
#include "rlz.h"

//------------------------------------------------------------------------------

/*
  A classic FM-index based on relative Lempel-Ziv parsing of the BWT.
  Neither BWT should not contain character value 1.
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

  uint64_t reportSize(bool print = false) const;
  void writeTo(const std::string& base_name) const;
  void writeTo(std::ostream& output) const;

  uint64_t size() const { return this->blocks.v.size(); }

  // Do not use with character value 0.
  template<class Iter> range_type find(Iter begin, Iter end) const;

  const SimpleFM<>&       reference;

  alphabet_type           alpha;
  relative_encoder        phrases;
  rlz_helper              blocks;
  rlz_helper*             block_rank;
  int_vector<8>           mismatches;

private:
  // Counts the number of occurrences of c in reference[ref_pos, ref_pos + phrase_length - 1].
  inline uint64_t countOf(uint64_t ref_pos, uint64_t phrase_length, uint8_t c) const
  {
    if(phrase_length == 0) { return 0; }
    return this->reference.bwt.rank(ref_pos + phrase_length, c) - this->reference.bwt.rank(ref_pos, c);
  }

  uint64_t rank(uint64_t i, uint8_t c) const;

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

#endif // _RELATIVE_FM_RLZ_FM_H
