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

#include "relative_fm.h"

using namespace relative;

//------------------------------------------------------------------------------

template<class Index>
void testIndex(std::string name, Index& index, std::vector<std::string>& patterns, size_type chars, size_type tags);

template<class ReferenceBWTType, class SequenceType>
void testIndex(std::string ref_name, std::string seq_name,
  std::string name, std::vector<std::string>& patterns, size_type chars, size_type tags);

const size_type SEQ_UNKNOWN         = 0;
const size_type SEQ_WT_PLAIN        = 1;
const size_type SEQ_WT_RRR          = 2;
const size_type SEQ_REPRESENTATIONS = 3; // Largest representation identifier + 1.

const size_type TAG_BUILD_INDEXES   = 0x01;
const size_type TAG_NATIVE_FORMAT   = 0x02;
const size_type TAG_LOCATE          = 0x20;
const size_type TAG_VERIFY          = 0x40;

size_type seqRepresentation(char_type type);
std::string indexName(size_type representation);
std::string indexName(size_type ref_representation, size_type seq_representation);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 4)
  {
    std::cerr << "Usage: query_test [indexes|tags] ref seq patterns" << std::endl;
    std::cerr << std::endl;

    std::cerr << "Index type is identified by a string, where X is the representation for" << std::endl;
    std::cerr << "sequences in the index, and Y is the representation in the reference." << std::endl;
    std::cerr << "Tags toggle an option for all subsequent indexes." << std::endl;
    std::cerr << std::endl;

    std::cerr << "Index types:" << std::endl;
    std::cerr << "  sX   SimpleFM" << std::endl;
    std::cerr << "  rYX  RelativeFM" << std::endl;
    std::cerr << std::endl;

    std::cerr << "Tags:" << std::endl;
    std::cerr << "  b    Build relative indexes (default: off)" << std::endl;
    std::cerr << "  p    Load SimpleFMs in plain format (default: on)" << std::endl;
    std::cerr << "  n    Load SimpleFMs in native format (default: off)" << std::endl;
    std::cerr << "  f    Execute find() queries (default: on)" << std::endl;
    std::cerr << "  L    Execute locate() queries (default: off)" << std::endl;
    std::cerr << "  V    Verify the results of locate() using extract() (default: off)" << std::endl;
    // FIXME Add an option to write the built index?
    std::cerr << std::endl;

    std::cerr << "Sequence representations:" << std::endl;
    std::cerr << "  p    Plain bitvectors in a wavelet tree" << std::endl;
    std::cerr << "  r    RRR bitvectors in a wavelet tree" << std::endl;
    std::cerr << std::endl;

    return 1;
  }
  int ref_arg = argc - 3, seq_arg = argc - 2, pat_arg = argc - 1;

  std::cout << "Query test" << std::endl;
  std::cout << std::endl;

  std::cout << "Reference: " << argv[ref_arg] << std::endl;
  std::cout << "Sequence: " << argv[seq_arg] << std::endl;
  std::cout << "Patterns: " << argv[pat_arg] << std::endl;
  std::cout << std::endl;

  std::vector<std::string> patterns;
  size_type chars = readRows(argv[pat_arg], patterns, true);
  std::cout << "Read " << patterns.size() << " patterns of total length " << chars << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  size_type tags = 0;
  for(int i = 1; i < ref_arg; i++)
  {
    size_type ref_enc = 0, seq_enc = 0;
    switch(argv[i][0])
    {

    // SimpleFM
    case 's':
      seq_enc = seqRepresentation(argv[i][1]);
      switch(seq_enc)
      {
      case SEQ_WT_PLAIN:
        {
          SimpleFM<> seq(argv[seq_arg], tags & TAG_NATIVE_FORMAT);
          testIndex(indexName(seq_enc), seq, patterns, chars, tags);
        }
        break;
      case SEQ_WT_RRR:
        {
          SimpleFM<sdsl::wt_huff<sdsl::rrr_vector<63> > > seq(argv[seq_arg], tags & TAG_NATIVE_FORMAT);
          testIndex(indexName(seq_enc), seq, patterns, chars, tags);
        }
        break;
      default:
        std::cerr << "query_test: Invalid index type: " << argv[i] << std::endl;
        break;
      }
      break;  // SimpleFM

    // RelativeFM
    case 'r':
      ref_enc = seqRepresentation(argv[i][1]); seq_enc = seqRepresentation(argv[i][2]);
      switch(ref_enc * SEQ_REPRESENTATIONS + seq_enc)
      {
      case SEQ_WT_PLAIN * SEQ_REPRESENTATIONS + SEQ_WT_PLAIN:
        testIndex<bwt_type, bwt_type>(argv[ref_arg], argv[seq_arg], indexName(ref_enc, seq_enc), patterns, chars, tags);
        break;
      case SEQ_WT_PLAIN * SEQ_REPRESENTATIONS + SEQ_WT_RRR:
        testIndex<bwt_type, sdsl::wt_huff<sdsl::rrr_vector<63> > >(argv[ref_arg], argv[seq_arg], indexName(ref_enc, seq_enc), patterns, chars, tags);
        break;

      case SEQ_WT_RRR * SEQ_REPRESENTATIONS + SEQ_WT_PLAIN:
        testIndex<sdsl::wt_huff<sdsl::rrr_vector<63> >, bwt_type>(argv[ref_arg], argv[seq_arg], indexName(ref_enc, seq_enc), patterns, chars, tags);
        break;
      case SEQ_WT_RRR * SEQ_REPRESENTATIONS + SEQ_WT_RRR:
        testIndex<sdsl::wt_huff<sdsl::rrr_vector<63> >, sdsl::wt_huff<sdsl::rrr_vector<63> > >(argv[ref_arg], argv[seq_arg], indexName(ref_enc, seq_enc), patterns, chars, tags);
        break;

      default:
        std::cerr << "query_test: Invalid index type: " << argv[i] << std::endl;
        break;
      }
      break;  // RelativeFM

    // Tags and the default case.
    case 'b':
      tags ^= TAG_BUILD_INDEXES;
      if(tags & TAG_BUILD_INDEXES) { std::cout << "Index construction activated." << std::endl; }
      else { std::cout << "Index construction deactivated." << std::endl; }
      std::cout << std::endl;
      break;
    case 'p':
      tags &= ~TAG_NATIVE_FORMAT;
      std::cout << "Loading indexes in plain format." << std::endl;
      std::cout << std::endl;
      break;
    case 'n':
      tags |= TAG_NATIVE_FORMAT;
      std::cout << "Loading indexes in native format." << std::endl;
      std::cout << std::endl;
      break;
    case 'f':
      tags &= ~TAG_LOCATE;
      std::cout << "Executing find() queries." << std::endl;
      std::cout << std::endl;
      break;
    case 'L':
      tags |= TAG_LOCATE;
      std::cout << "Executing locate() queries." << std::endl;
      std::cout << std::endl;
      break;
    case 'V':
      tags ^= TAG_VERIFY;
      if(tags & TAG_VERIFY) { std::cout << "Verifying the results with extract() queries." << std::endl; }
      else { std::cout << "Switching off verification." << std::endl; }
      std::cout << std::endl;
      break;
    default:
      std::cerr << "query_test: Invalid index type: " << argv[i] << std::endl;
      break;
    }
  }

  double memory = inMegabytes(memoryUsage());
  std::cout << std::endl;
  std::cout << "Memory usage: " << memory << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

template<class Index>
void
testIndex(std::string name, Index& index, std::vector<std::string>& patterns, size_type chars, size_type tags)
{
  bool locate = tags & TAG_LOCATE, verify = tags & TAG_VERIFY;
  if(!locate) { verify = false; }
  if(locate && !(index.supportsLocate(true))) { locate = false; }
  if(verify && !(index.supportsExtract(true))) { verify = false; }

  double start = readTimer();
  size_type found = 0, matches = 0, failed = 0;
  size_type hash = FNV_OFFSET_BASIS;
  for(auto pattern : patterns)
  {
    range_type res = index.find(pattern.begin(), pattern.end());
    if(Range::length(res) > 0)
    {
      found++; matches += Range::length(res);
      if(locate)
      {
        if(verify)
        {
          bool fail = false;
          size_type temp = 0;
          if(res.first > 0) // The suffix before the range.
          {
            temp = index.locate(res.first - 1);
            if(index.extract(temp, temp + pattern.length() - 1) >= pattern) { fail = true; }
          }
          { // The first suffix.
            temp = index.locate(res.first); hash = fnv1a_hash(index.locate(res.first), hash);
            if(index.extract(temp, temp + pattern.length() - 1) != pattern) { fail = true; }
          }
          for(size_type i = res.first + 2; i <= res.second; i++) // The middle suffixes.
          {
            hash = fnv1a_hash(index.locate(i - 1), hash);
          }
          if(res.second > res.first)  // The last suffix.
          {
            temp = index.locate(res.second); hash = fnv1a_hash(index.locate(res.second), hash);
            if(index.extract(temp, temp + pattern.length() - 1) != pattern) { fail = true; }
          }
          if(res.second + 1 < index.size()) // The suffix after the range.
          {
            temp = index.locate(res.second + 1);
            if(index.extract(temp, temp + pattern.length() - 1) <= pattern) { fail = true; }
          }
          if(fail)
          {
#ifdef VERBOSE_STATUS_INFO
            if(failed == 0)
            {
              std::cerr << "testIndex(): Verification failed for the following patterns:" << std::endl;
            }
            std::cerr << pattern << std::endl;
#endif
            failed++;
          }
        }
        else
        {
          for(size_type i = res.first; i <= res.second; i++) { hash = fnv1a_hash(index.locate(i), hash); }
        }
      }
    }
  }
  double seconds = readTimer() - start;

  printSize(name, index.reportSize(), index.size());
  printTime(name, found, matches, chars, seconds, locate);
  if(locate) { std::cout << "Hash of located positions: " << hash << std::endl; }
  if(verify && failed > 0) { std::cout << "Verification failed for " << failed << " patterns." << std::endl; }
  std::cout << std::endl;
}

template<class ReferenceBWTType, class SequenceType>
void
testIndex(std::string ref_name, std::string seq_name,
  std::string name, std::vector<std::string>& patterns, size_type chars, size_type tags)
{
  SimpleFM<ReferenceBWTType> ref(ref_name, tags & TAG_NATIVE_FORMAT);

  if(tags & TAG_BUILD_INDEXES)
  {
    align_parameters parameters;
    SimpleFM<ReferenceBWTType> seq(seq_name, tags & TAG_NATIVE_FORMAT);
    RelativeFM<ReferenceBWTType, SequenceType> rfm(ref, seq, parameters);
    testIndex(name, rfm, patterns, chars, tags);
  }
  else
  {
    RelativeFM<ReferenceBWTType, SequenceType> seq(ref, seq_name);
    testIndex(name, seq, patterns, chars, tags);
  }
}

size_type
seqRepresentation(char_type type)
{
  switch(type)
  {
    case 'p': return SEQ_WT_PLAIN;
    case 'r': return SEQ_WT_RRR;
  }
  return SEQ_UNKNOWN;
}

std::string
encodingName(size_type representation)
{
  switch(representation)
  {
    case SEQ_WT_PLAIN:   return "plain";
    case SEQ_WT_RRR:     return "rrr";
  }
  return "unknown";
}

std::string
indexName(size_type representation)
{
  switch(representation)
  {
    case SEQ_WT_PLAIN:   return "SimpleFM<plain>";
    case SEQ_WT_RRR:     return "SimpleFM<rrr>";
  }
  return "Unknown";
}

std::string
indexName(size_type ref_representation, size_type seq_representation)
{
  return "RFM<" + encodingName(ref_representation) + "," + encodingName(seq_representation) + ">";
}

//------------------------------------------------------------------------------
