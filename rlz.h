/*
  Copyright (c) 2015, 2016 Genome Research Ltd.
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

#ifndef _RELATIVE_FM_RLZ_H
#define _RELATIVE_FM_RLZ_H

#include <vector>

#include <sdsl/suffix_arrays.hpp>

#include "support.h"

namespace relative
{

//------------------------------------------------------------------------------

template<class IntVector, bool differential>
struct CharAt
{
  inline static size_type at(const IntVector& seq, size_type i)
  {
    if(i >= seq.size()) { return 0; }
    return seq[i];
  }
};

template<class IntVector>
struct CharAt<IntVector, true>
{
  inline static size_type at(const IntVector& seq, size_type i)
  {
    if(i >= seq.size()) { return 0; }
    size_type prev = (i > 0 ? seq[i - 1] : 0), curr = seq[i];
    return DiffEncoderNZ::encode(curr, prev);
  }
};

//------------------------------------------------------------------------------

// FIXME Handle this in a more general way.
const size_type MAX_RLZ_PHRASE_LENGTH = 1024;

/*
  This version is intended for integer sequences. It parses the reference using a prebuilt
  suffix array. The text must not contain character value 0, while the reference has it only
  as an implicit endmarker.
  If mismatches == 0, the algorithm does not output them.
  If differential == true, the parsing is based on differential values instead of absolute
  values.
*/
template<class IntVector, bool differential>
void
relativeLZ(const IntVector& text, const IntVector& reference, const sdsl::int_vector<0>& sa,
  std::vector<size_type>& starts, std::vector<size_type>& lengths, std::vector<size_type>* mismatches)
{
  if(text.size() == 0) { return; }

  starts.clear(); lengths.clear();
  if(mismatches != 0) { mismatches->clear(); }
  size_type text_pos = 0;
  while(text_pos < text.size())
  {
    size_type sp = 0, ep = sa.size() - 1, matched = 0;
    size_type limit = std::min(text.size() - text_pos, MAX_RLZ_PHRASE_LENGTH);
    while(matched < limit)
    {
      size_type low = sp, high = ep, next = CharAt<IntVector, differential>::at(text, text_pos + matched);
      while(low < high) // Find the first suffix that matches the next character.
      {
        size_type mid = low + (high - low) / 2;
        size_type val = CharAt<IntVector, differential>::at(reference, sa[mid] + matched);
        if(val < next) { low = mid + 1; }
        else if(val > next) { high = mid - 1; }
        else { high = mid; }
      }
      if(CharAt<IntVector, differential>::at(reference, sa[low] + matched) != next) { break; }
      sp = low;

      high = ep;
      while(low < high) // Find the last suffix that matches the next character.
      {
        size_type mid = low + (high + 1 - low) / 2;
        size_type val = CharAt<IntVector, differential>::at(reference, sa[mid] + matched);
        if(val > next) { high = mid - 1; }
        else { low = mid; }
      }
      ep = high; matched++;
    }

    starts.push_back(sa[sp]); // FIXME: Find the match closest to text_pos?
    if(text_pos + matched < text.size()) { matched++; } // Add the mismatch.
    lengths.push_back(matched);
    if(mismatches != 0) { mismatches->push_back(text[text_pos + matched - 1]); }
    text_pos += matched;
  }
}

//------------------------------------------------------------------------------

} // namespace relative

#endif // _RELATIVE_FM_RLZ_H
