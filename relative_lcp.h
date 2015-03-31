#ifndef _RELATIVE_FM_RELATIVE_LCP_H
#define _RELATIVE_FM_RELATIVE_LCP_H

#include "rlz.h"
#include "support.h"

namespace relative
{

//------------------------------------------------------------------------------

class RelativeLCP
{
public:
  typedef SLArray       lcp_type;
  typedef int_vector<0> dlcp_type;
  typedef int_vector<0> index_type;

  const static uint64_t BRANCHING_FACTOR = 4; // For the minimal LCP value tree.
  const static std::string EXTENSION; // .rlcp

//------------------------------------------------------------------------------

  // Reference is an LCP array.
  RelativeLCP(const lcp_type& ref, const lcp_type& seq,
              const index_type& ref_sa, bool print = false);
  RelativeLCP(const lcp_type& ref, const std::string& base_name);
  RelativeLCP(const lcp_type& ref, std::istream& input);
  ~RelativeLCP();

  inline uint64_t size() const { return this->blocks.sum(); }

  uint64_t reportSize(bool print = false) const;
  void writeTo(const std::string& base_name) const;
  void writeTo(std::ostream& output) const;

//------------------------------------------------------------------------------

  inline uint64_t operator[] (uint64_t i) const
  {
    if(i >= this->size()) { return 0; }

    bool is_last = false; // Last character in the phrase.
    uint64_t phrase = this->blocks.inverse(i, is_last); // Phrase is 0-based.
    if(is_last) { return this->samples[phrase + 1]; }

    // Starting position of the phrase in seq and ref.
    uint64_t seq_pos = this->blocks.sum(phrase);
    uint64_t ref_pos = this->phrases.decode(phrase, seq_pos);

    uint64_t res = this->samples[phrase];
    if(ref_pos > 0) { res -= this->reference[ref_pos - 1]; }
    res += this->reference[ref_pos + i - seq_pos];
    return res;
  }

  // FIXME implement extracting multiple cells.

  // FIXME what these should return?
  // FIXME implement
  uint64_t previousSmaller(uint64_t i) const;
  uint64_t nextSmaller(uint64_t i) const;
  uint64_t minimumValue(range_type range) const;
  uint64_t minimumValue(uint64_t from, uint64_t to) const;

//------------------------------------------------------------------------------

  /*
    We encode (val - prev) as a positive integer. SA construction does not work with
    character value 0.
  */
  inline static uint64_t encode(uint64_t val, uint64_t prev)
  {
    if(val >= prev) { return 2 * (val - prev) + 1; }
    else            { return 2 * (prev - val); }
  }

  inline static uint64_t decode(uint64_t code, uint64_t prev)
  {
    if(code & 1) { return (code - 1) / 2 + prev; }
    else         { return prev - code / 2; }
  }

  // Use endmarker == true, if you are going to build an index for the DLCP array.
  static void differentialArray(const lcp_type& lcp, dlcp_type& dlcp, bool endmarker);

//------------------------------------------------------------------------------

  const lcp_type&  reference;
  relative_encoder phrases;
  CumulativeArray  blocks;
  SLArray          samples;
  SLArray          tree;

//------------------------------------------------------------------------------

private:
  void loadFrom(std::istream& input);
  void absoluteSamples(const lcp_type& lcp, const std::vector<uint64_t>& lengths);

//------------------------------------------------------------------------------

  RelativeLCP();
  RelativeLCP(const RelativeLCP&);
  RelativeLCP(RelativeLCP&&);
  RelativeLCP& operator=(const RelativeLCP&);
  RelativeLCP& operator==(RelativeLCP&);
};

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_RELATIVE_LCP_H
