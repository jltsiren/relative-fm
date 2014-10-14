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
  int_vector<0>           phrases;    // Phrase starts encoded as (ref_pos - text_pos).
  bit_vector              phrase_rle; // phrase_rle[i] is set, if phrase i is stored.
  bit_vector::rank_1_type phrase_rank;
  rlz_helper              blocks;
  rlz_helper*             block_rank;
  int_vector<8>           mismatches;

private:
  // Converts text positions into reference positions, assuming that they are
  // within copied substrings.
  inline uint64_t refPos(uint64_t phrase, uint64_t text_pos) const
  {
    uint64_t temp = this->phrases[this->phrase_rank(phrase + 1) - 1];
    if(temp & 1) { return text_pos - (temp >> 1); }
    return (temp >> 1) + text_pos;
  }

  // Encodes (val - ref) as unsigned integer.
  inline static uint64_t relativeEncoding(uint64_t val, uint64_t ref)
  {
    if(val >= ref) { return (val - ref) << 1; }
    return ((ref - val) << 1) | 1;
  }

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
    if(!hasChar(this->alpha, *begin)) { return range_type(1, 0); }
    uint64_t pos = cumulative(this->alpha, *begin);
    res.first = pos + this->rank(res.first, *begin);
    res.second = pos + this->rank(res.second + 1, *begin) - 1;
    if(isEmpty(res)) { return range_type(1, 0); }
    ++begin;
  }
  return res;
}

//------------------------------------------------------------------------------

#endif // _RELATIVE_FM_RLZ_FM_H
