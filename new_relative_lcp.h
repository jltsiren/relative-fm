/*
  Copyright (c) 2015, 2016, 2017 Genome Research Ltd.

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

#ifndef _NEW_RELATIVE_FM_RELATIVE_LCP_H
#define _NEW_RELATIVE_FM_RELATIVE_LCP_H

#include <rlz/lcp/api.hpp>

#include "support.h"

namespace relative
{

//------------------------------------------------------------------------------

class NewRelativeLCP
{
public:
  typedef SLArray              lcp_type;
  typedef lcp_type::size_type  size_type;
  typedef lcp_type::value_type value_type;

  typedef rlz::lcp::LcpIndex<lcp_type::iterator, lcp_type::value_type> rlcp_type;

  const static size_type BRANCHING_FACTOR = 64; // For the range minima tree.
  const static std::string EXTENSION; // .rlcp

//------------------------------------------------------------------------------

  // Reference is an LCP array.
  NewRelativeLCP(const lcp_type& ref, const lcp_type& seq);
  NewRelativeLCP(const lcp_type& ref, const std::string& base_name);
  NewRelativeLCP(const lcp_type& ref, std::istream& input);
  ~NewRelativeLCP();

  size_type reportSize(bool print = false) const;
  void writeTo(const std::string& base_name) const;
  void writeTo(std::ostream& output) const;

//------------------------------------------------------------------------------

  inline size_type size() const { return this->array.size(); }

  // For range minima tree.
  inline size_type phrases() const { return this->array.phrases(); }
  inline size_type values() const { return this->tree.size(); }
  inline size_type levels() const { return this->offsets.size() - 1; }
  inline size_type branching() const { return this->branching_factor; }

  inline value_type operator[] (size_type i) const { return this->array(i); }

//------------------------------------------------------------------------------

  /*
    The return value is (res, LCP[res]) or notFound(). RMQ always returns the leftmost
    minimum value.

    FIXME implement
  */

  range_type psv(size_type pos) const;
  range_type psev(size_type pos) const;

  range_type nsv(size_type pos) const;
  range_type nsev(size_type pos) const;

  range_type rmq(size_type sp, size_type ep) const;
  range_type rmq(range_type range) const;

  // Returned when a psv/nsv/rmq query cannot find a suitable value.
  inline range_type notFound() const { return range_type(this->size() + this->values(), this->size()); }

//------------------------------------------------------------------------------

  const lcp_type&      reference;
  rlcp_type            array;

  size_type            branching_factor;
  lcp_type             tree;
  sdsl::int_vector<64> offsets;

//------------------------------------------------------------------------------

private:
  void loadFrom(std::istream& input);

//------------------------------------------------------------------------------

  NewRelativeLCP(const NewRelativeLCP&) = delete;
  NewRelativeLCP& operator=(const NewRelativeLCP&) = delete;
};

//------------------------------------------------------------------------------

} // namespace relative

#endif // _NEW_RELATIVE_FM_RELATIVE_LCP_H
