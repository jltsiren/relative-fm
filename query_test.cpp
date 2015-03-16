#include "relative_fm.h"
#include "rlz_fm.h"
#include "rlz_vector.h"
#include "sequence.h"

using namespace relative;

//------------------------------------------------------------------------------

template<class Index>
void testIndex(std::string name, Index& index, std::vector<std::string>& patterns, uint64_t chars, uint64_t tags);

template<class ReferenceBWTType, class SequenceType>
void testIndex(std::string ref_name, std::string seq_name,
  std::string name, std::vector<std::string>& patterns, uint64_t chars, uint64_t tags);

const uint64_t SEQ_UNKNOWN    = 0;
const uint64_t SEQ_WT_PLAIN   = 1;
const uint64_t SEQ_WT_RRR     = 2;
const uint64_t SEQ_WT_RLZ     = 3;
const uint64_t SEQ_RLSEQUENCE = 4;
const uint64_t SEQ_REPRESENTATIONS = 5; // Largest representation identifier + 1.

const uint64_t TAG_ROPEBWT2_ALPHABET = 0x01;
const uint64_t TAG_BUILD_INDEXES     = 0x02;
const uint64_t TAG_NATIVE_FORMAT     = 0x04;
const uint64_t TAG_ROPEBWT2_FORMAT   = 0x08;
const uint64_t TAG_LOCATE            = 0x10;

uint64_t seqRepresentation(uint8_t type);
std::string indexName(uint64_t representation);
std::string indexName(uint64_t ref_representation, uint64_t seq_representation);

LoadMode getMode(uint64_t tags);
std::string modeName(LoadMode mode);

// Outputs warnings for suspicious tag combinations to cout.
void sanityCheck(uint64_t tags);

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
    std::cerr << "  lY   RLZFM" << std::endl;
    // FIXME Different encodings are currently not supported for RLZFM.
    std::cerr << std::endl;
    std::cerr << "Tags:" << std::endl;
    std::cerr << "  a    Use ropebwt2 alphabet (default: off)" << std::endl;
    std::cerr << "  b    Build relative indexes (default: off)" << std::endl;
    std::cerr << "  L    Execute locate() queries (default: off)" << std::endl;
    std::cerr << "  n    Load SimpleFMs in native format (default: off)" << std::endl;
    std::cerr << "  2    Load SimpleFMs in ropebwt2 format (default: off)" << std::endl;
    // FIXME Add an option to write the built index?
    std::cerr << std::endl;
    std::cerr << "Sequence representations:" << std::endl;
    std::cerr << "  p    Plain bitvectors in a wavelet tree" << std::endl;
    std::cerr << "  r    RRR bitvectors in a wavelet tree" << std::endl;
    std::cerr << "  l    RLZ bitvectors in a wavelet tree (only for SimpleFM)" << std::endl;
    std::cerr << "  R    RLSequence" << std::endl;
    std::cerr << std::endl;
    // FIXME RLZ bitvectors may need native format for both ref and seq.
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
  uint64_t chars = readRows(argv[pat_arg], patterns, true);
  std::cout << "Read " << patterns.size() << " patterns of total length " << chars << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  uint64_t tags = 0;
  LoadMode mode = mode_plain;
  for(int i = 1; i < ref_arg; i++)
  {
    uint64_t ref_enc = 0, seq_enc = 0;
    switch(argv[i][0])
    {

    // SimpleFM
    case 's':
      sanityCheck(tags);
      seq_enc = seqRepresentation(argv[i][1]);
      switch(seq_enc)
      {
      case SEQ_WT_PLAIN:
        {
          SimpleFM<> seq(argv[seq_arg], mode);
          testIndex(indexName(seq_enc), seq, patterns, chars, tags);
        }
        break;
      case SEQ_WT_RRR:
        {
          SimpleFM<wt_huff<rrr_vector<63> > > seq(argv[seq_arg], mode);
          testIndex(indexName(seq_enc), seq, patterns, chars, tags);
        }
        break;
      case SEQ_WT_RLZ:
        {
          SimpleFM<> ref(argv[ref_arg], mode);
          SimpleFM<wt_huff<rlz_vector> > seq(argv[seq_arg], mode);
          bit_vector::rank_1_type b_r(&(ref.bwt.bv));
          bit_vector::select_1_type b_s1(&(ref.bwt.bv));
          bit_vector::select_0_type b_s0(&(ref.bwt.bv));
          {
            rlz_vector& temp = const_cast<rlz_vector&>(seq.bwt.bv);
            temp.compress(ref.bwt.bv, b_r, b_s1, b_s0);
          }
          testIndex(indexName(seq_enc), seq, patterns, chars, tags);
        }
        break;
      case SEQ_RLSEQUENCE:
        {
          SimpleFM<RLSequence> seq(argv[seq_arg], mode);
          testIndex(indexName(seq_enc), seq, patterns, chars, tags);
          if(mode == mode_ropebwt2)
          {
            uint64_t hash_value = seq.bwt.hash(seq.alpha);
            std::cerr << "BWT imported from ropebwt2 (hash value " << hash_value << ")" << std::endl;
            for(uint64_t c = 0; c < seq.alpha.sigma; c++)
            {
              std::cerr << "  Character " << c << " (" << (char)(seq.alpha.comp2char[c])
                        << "), count " << (seq.alpha.C[c + 1] - seq.alpha.C[c]) << std::endl;
            }
          }
        }
        break;
      default:
        std::cerr << "query_test: Invalid index type: " << argv[i] << std::endl;
        break;
      }
      break;  // SimpleFM

    // RelativeFM
    case 'r':
      sanityCheck(tags);
      ref_enc = seqRepresentation(argv[i][1]); seq_enc = seqRepresentation(argv[i][2]);
      switch(ref_enc * SEQ_REPRESENTATIONS + seq_enc)
      {
      case SEQ_WT_PLAIN * SEQ_REPRESENTATIONS + SEQ_WT_PLAIN:
        testIndex<bwt_type, bwt_type>(argv[ref_arg], argv[seq_arg], indexName(ref_enc, seq_enc), patterns, chars, tags);
        break;
      case SEQ_WT_PLAIN * SEQ_REPRESENTATIONS + SEQ_WT_RRR:
        testIndex<bwt_type, wt_huff<rrr_vector<63> > >(argv[ref_arg], argv[seq_arg], indexName(ref_enc, seq_enc), patterns, chars, tags);
        break;
      case SEQ_WT_PLAIN * SEQ_REPRESENTATIONS + SEQ_RLSEQUENCE:
        testIndex<bwt_type, RLSequence>(argv[ref_arg], argv[seq_arg], indexName(ref_enc, seq_enc), patterns, chars, tags);
        break;

      case SEQ_WT_RRR * SEQ_REPRESENTATIONS + SEQ_WT_PLAIN:
        testIndex<wt_huff<rrr_vector<63> >, bwt_type>(argv[ref_arg], argv[seq_arg], indexName(ref_enc, seq_enc), patterns, chars, tags);
        break;
      case SEQ_WT_RRR * SEQ_REPRESENTATIONS + SEQ_WT_RRR:
        testIndex<wt_huff<rrr_vector<63> >, wt_huff<rrr_vector<63> > >(argv[ref_arg], argv[seq_arg], indexName(ref_enc, seq_enc), patterns, chars, tags);
        break;
      case SEQ_WT_RRR * SEQ_REPRESENTATIONS + SEQ_RLSEQUENCE:
        testIndex<wt_huff<rrr_vector<63> >, RLSequence>(argv[ref_arg], argv[seq_arg], indexName(ref_enc, seq_enc), patterns, chars, tags);
        break;

      case SEQ_RLSEQUENCE * SEQ_REPRESENTATIONS + SEQ_WT_PLAIN:
        testIndex<RLSequence, bwt_type>(argv[ref_arg], argv[seq_arg], indexName(ref_enc, seq_enc), patterns, chars, tags);
        break;
      case SEQ_RLSEQUENCE * SEQ_REPRESENTATIONS + SEQ_WT_RRR:
        testIndex<RLSequence, wt_huff<rrr_vector<63> > >(argv[ref_arg], argv[seq_arg], indexName(ref_enc, seq_enc), patterns, chars, tags);
        break;
      case SEQ_RLSEQUENCE * SEQ_REPRESENTATIONS + SEQ_RLSEQUENCE:
        testIndex<RLSequence, RLSequence>(argv[ref_arg], argv[seq_arg], indexName(ref_enc, seq_enc), patterns, chars, tags);
        break;

      default:
        std::cerr << "query_test: Invalid index type: " << argv[i] << std::endl;
        break;
      }
      break;  // RelativeFM

    // RLZFM
    case 'l':
      sanityCheck(tags);
      {
        SimpleFM<> ref(argv[ref_arg], mode);
        if(tags & TAG_BUILD_INDEXES)
        {
          SimpleFM<> seq(argv[seq_arg], mode);
          if(tags & TAG_ROPEBWT2_ALPHABET)
          {
            ref.alpha.assign(ROPEBWT2_ALPHABET); seq.alpha.assign(ROPEBWT2_ALPHABET);
          }
          RLZFM rlz(ref, seq);
          testIndex("RLZFM", rlz, patterns, chars, tags);
        }
        else
        {
          RLZFM seq(ref, argv[seq_arg]);
          testIndex("RLZFM", seq, patterns, chars, tags);
        }
      }
      break;  // RLZFM

    // Tags and the default case.
    case 'a':
      tags ^= TAG_ROPEBWT2_ALPHABET;
      if(tags & TAG_ROPEBWT2_ALPHABET) { std::cout << "Using ropebwt2 alphabet." << std::endl; }
      else { std::cout << "No longer using ropebwt2 alphabet" << std::endl; }
      std::cout << std::endl;
      break;
    case 'b':
      tags ^= TAG_BUILD_INDEXES;
      if(tags & TAG_BUILD_INDEXES) { std::cout << "Index construction activated." << std::endl; }
      else { std::cout << "Index construction deactivated." << std::endl; }
      std::cout << std::endl;
      break;
    case 'L':
      tags ^= TAG_LOCATE;
      if(tags & TAG_LOCATE) { std::cout << "Executing locate() queries." << std::endl; }
      else { std::cout << "Executing find() queries." << std::endl; }
      std::cout << std::endl;
      break;
    case 'n':
      tags ^= TAG_NATIVE_FORMAT; tags &= ~TAG_ROPEBWT2_FORMAT;
      mode = getMode(tags);
      std::cout << "Loading SimpleFMs in " << modeName(mode) << " format." << std::endl;
      std::cout << std::endl;
      break;
    case '2':
      tags ^= TAG_ROPEBWT2_FORMAT; tags &= ~TAG_NATIVE_FORMAT;
      mode = getMode(tags);
      std::cout << "Loading SimpleFMs in " << modeName(mode) << " format." << std::endl;
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
testIndex(std::string name, Index& index, std::vector<std::string>& patterns, uint64_t chars, uint64_t tags)
{
  if(tags & TAG_ROPEBWT2_ALPHABET)
  {
    index.alpha.assign(ROPEBWT2_ALPHABET);
  }

  bool locate = (tags & TAG_LOCATE);
  if(locate && !(index.supportsLocate(true))) { locate = false; }

  double start = readTimer();
  uint64_t found = 0, matches = 0;
  uint64_t hash = FNV_OFFSET_BASIS;
  for(auto pattern : patterns)
  {
    range_type res = index.find(pattern.begin(), pattern.end());
    if(length(res) > 0)
    {
      found++; matches += length(res);
      if(locate)
      {
        for(uint64_t i = res.first; i <= res.second; i++)
        {
          uint64_t temp = index.locate(i);
//          std::cerr << "SA[" << i << "] = " << temp << std::endl;
          hash = fnv1a_hash(temp, hash);
        }
      }
    }
  }
  double seconds = readTimer() - start;

  printSize(name, index.reportSize(), index.size());
  printTime(name, found, matches, chars, seconds, locate);
  if(locate) { std::cout << "Hash of located positions: " << hash << std::endl; }
  std::cout << std::endl;
}

template<class ReferenceBWTType, class SequenceType>
void
testIndex(std::string ref_name, std::string seq_name,
  std::string name, std::vector<std::string>& patterns, uint64_t chars, uint64_t tags)
{
  LoadMode mode = getMode(tags);
  SimpleFM<ReferenceBWTType> ref(ref_name, mode);
  if(tags & TAG_ROPEBWT2_ALPHABET) { ref.alpha.assign(ROPEBWT2_ALPHABET); }

  if(tags & TAG_BUILD_INDEXES)
  {
    align_parameters parameters; parameters.sorted_alphabet = !(tags & TAG_ROPEBWT2_ALPHABET);
    SimpleFM<ReferenceBWTType> seq(seq_name, mode);
    if(tags & TAG_ROPEBWT2_ALPHABET) { seq.alpha.assign(ROPEBWT2_ALPHABET); }
    RelativeFM<ReferenceBWTType, SequenceType> rfm(ref, seq, parameters);
    testIndex(name, rfm, patterns, chars, tags);
  }
  else
  {
    RelativeFM<ReferenceBWTType, SequenceType> seq(ref, seq_name);
    testIndex(name, seq, patterns, chars, tags);
  }
}

uint64_t
seqRepresentation(uint8_t type)
{
  switch(type)
  {
    case 'p': return SEQ_WT_PLAIN;
    case 'r': return SEQ_WT_RRR;
    case 'l': return SEQ_WT_RLZ;
    case 'R': return SEQ_RLSEQUENCE;
  }
  return SEQ_UNKNOWN;
}

std::string
encodingName(uint64_t representation)
{
  switch(representation)
  {
    case SEQ_WT_PLAIN:   return "plain";
    case SEQ_WT_RRR:     return "rrr";
    case SEQ_WT_RLZ:     return "rlz";
    case SEQ_RLSEQUENCE: return "rle";
  }
  return "unknown";
}

std::string
indexName(uint64_t representation)
{
  switch(representation)
  {
    case SEQ_WT_PLAIN:   return "SimpleFM<plain>";
    case SEQ_WT_RRR:     return "SimpleFM<rrr>";
    case SEQ_WT_RLZ:     return "SimpleFM<rlz>";
    case SEQ_RLSEQUENCE: return "SimpleFM<rle>";
  }
  return "Unknown";
}

std::string
indexName(uint64_t ref_representation, uint64_t seq_representation)
{
  return "RFM<" + encodingName(ref_representation) + "," + encodingName(seq_representation) + ">";
}

LoadMode
getMode(uint64_t tags)
{
  if(tags & TAG_ROPEBWT2_FORMAT) { return mode_ropebwt2; }
  else if(tags & TAG_NATIVE_FORMAT) { return mode_native; }
  else { return mode_plain; }
}

std::string
modeName(LoadMode mode)
{
  if(mode == mode_plain) { return "plain"; }
  else if(mode == mode_native) { return "native"; }
  else if(mode == mode_ropebwt2) { return "ropebwt2"; }
  else { return "unknown"; }
}

void
sanityCheck(uint64_t tags)
{
  if((tags & TAG_ROPEBWT2_FORMAT) && !(tags & TAG_ROPEBWT2_ALPHABET))
  {
    std::cout << "Warning: Using ropebwt2 format without ropebwt2 alphabet" << std::endl;
  }
}

//------------------------------------------------------------------------------
