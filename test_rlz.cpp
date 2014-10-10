#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>

#include "rlz.h"
#include "rlz_vector.h"
#include "utils.h"


// Which tests to run.
//#define RLZ_TESTS
//#define BV_TESTS
//#define WT_TESTS
#define CST_TESTS


// For RLZ tests.
const uint64_t RLZ_SIZE = 128 * (uint64_t)1048576;

// For bitvector tests.
const uint64_t BV_SIZE = 128 * (uint64_t)1048576;

// For wavelet tree tests.
const uint64_t WT_SIZE = 512 * (uint64_t)1048576;
const uint8_t  WT_ALPHABET = 4;

// For all tests.
const uint64_t CORRECTNESS_QUERIES = 100000;
const uint64_t TIMING_QUERIES = 10000000;
const uint64_t TOTAL_TIME_TO_NANOSECS = 100; // 10^9 / TIMING_QUERIES

// Bitvector generation parameters.
const uint64_t EXTRA_SIZE = 1024; // Some extra space to avoid problems in text generation.
const double   EXTENSION_PROB = 0.3;

//------------------------------------------------------------------------------

void testRLZ(int argc, char** argv);
void testBV(int argc, char** argv);
void testWT(int argc, char** argv);
void testCST(int argc, char** argv);

int
main(int argc, char** argv)
{
  std::cout << "Relative Lempel-Ziv tests" << std::endl;
  std::cout << std::endl;

#ifdef RLZ_TESTS
  testRLZ(argc, argv);
#endif

#ifdef BV_TESTS
  testBV(argc, argv);
#endif

#ifdef WT_TESTS
  testWT(argc, argv);
#endif

#ifdef CST_TESTS
  testCST(argc, argv);
#endif

  return 0;
}

//------------------------------------------------------------------------------

bit_vector generateReference(std::mt19937_64& rng, uint64_t size);
bit_vector generateVariant(std::mt19937_64& rng, const bit_vector& reference, double mutation_rate);

int_vector<8> generateReference(std::mt19937_64& rng, uint64_t size, uint8_t alphabet_size);
int_vector<8> generateVariant(std::mt19937_64& rng, const int_vector<8>& reference, double mutation_rate);

inline double probability(std::mt19937_64& rng)
{
  return rng() / (rng.max() + 1.0);
}

template<class A, class B, class C, class D>
uint64_t totalSize(const A& a, const B& b, const C& c, const D& d)
{
  return size_in_bytes(a) + size_in_bytes(b) + size_in_bytes(c) + size_in_bytes(d);
}

void
testRLZ(int argc, char** argv)
{
  std::cout << "Parsing tests" << std::endl;
  std::cout << std::endl;

  for(int arg = 1; arg < argc; arg++)
  {
    std::mt19937_64 rng(0xDEADBEEF);
    double prob = atof(argv[arg]);
    std::cout << "Mutation rate: " << prob << std::endl;

    bit_vector reference = generateReference(rng, RLZ_SIZE);
    bit_vector text = generateVariant(rng, reference, prob);
    std::cout << "Reference length " << reference.size() << ", text length " << text.size() << "." << std::endl;

    std::vector<uint64_t> starts, lengths;
    bit_vector mismatches;
    double start_time = readTimer();
    relativeLZ(text, reference, starts, lengths, mismatches);
    double seconds = readTimer() - start_time;
    std::cout << "Parsing took " << seconds << " seconds, " << inMegabytes(memoryUsage()) << " MB." << std::endl;

    uint64_t errors = 0;
    for(uint64_t phrase = 0, text_pos = 0; phrase < starts.size(); phrase++)
    {
      for(uint64_t ref_pos = starts[phrase]; ref_pos < starts[phrase] + lengths[phrase] - 1; text_pos++, ref_pos++)
      {
        if(text[text_pos] != reference[ref_pos]) { errors++; }
      }
      if(text[text_pos] != mismatches[phrase]) { errors++; } text_pos++;
    }
    std::cout << "Decompressed the text with " << errors << " error(s)." << std::endl;

    std::cout << std::endl;
  }
}

void
testBV(int argc, char** argv)
{
  std::cout << "Bitvector tests" << std::endl;
  std::cout << std::endl;

  for(int arg = 1; arg < argc; arg++)
  {
    std::mt19937_64 rng(0xDEADBEEF);
    double prob = atof(argv[arg]);
    std::cout << "Mutation rate: " << prob << std::endl;

    bit_vector reference = generateReference(rng, BV_SIZE);
    bit_vector text = generateVariant(rng, reference, prob);
    uint64_t onebits = util::cnt_one_bits(text);
    std::cout << "Reference length " << reference.size() << ", text length " << text.size() << "." << std::endl;

    bit_vector::rank_1_type   ref_rank(&reference);
    bit_vector::select_1_type ref_select_1(&reference);
    bit_vector::select_0_type ref_select_0(&reference);
    std::cout << "Reference bitvector: " <<
      inBPC(totalSize(reference, ref_rank, ref_select_1, ref_select_0), reference.size()) << " bpc" << std::endl;

    bit_vector::rank_1_type   text_rank(&text);
    bit_vector::select_1_type text_select_1(&text);
    bit_vector::select_0_type text_select_0(&text);
    std::cout << "Plain bitvector: " <<
      inBPC(totalSize(text, text_rank, text_select_1, text_select_0), text.size()) << " bpc" << std::endl;

    RLZVector relative(text, reference, ref_rank, ref_select_1, ref_select_0);
    std::cout << "Relative bitvector: " << inBPC(relative.reportSize(), text.size()) << " bpc" << std::endl;

    uint64_t errors = 0;
    for(uint64_t i = 0; i < CORRECTNESS_QUERIES; i++)
    {
      uint64_t pos = rand() % text.size();
      if(text_rank(pos) != relative.rank(pos)) { errors++; }
    }
    std::cout << "Completed " << CORRECTNESS_QUERIES << " rank queries with " << errors << " error(s)." << std::endl;

    errors = 0;
    for(uint64_t i = 0; i < CORRECTNESS_QUERIES; i++)
    {
      uint64_t pos = rand() % onebits + 1;
      if(text_select_1(pos) != relative.select_1(pos)) { errors++; }
    }
    std::cout << "Completed " << CORRECTNESS_QUERIES << " select_1 queries with " << errors << " error(s)." << std::endl;

    errors = 0;
    for(uint64_t i = 0; i < CORRECTNESS_QUERIES; i++)
    {
      uint64_t pos = rand() % (text.size() - onebits) + 1;
      if(text_select_0(pos) != relative.select_0(pos)) { errors++; }
    }
    std::cout << "Completed " << CORRECTNESS_QUERIES << " select_0 queries with " << errors << " error(s)." << std::endl;

    errors = 0;
    for(uint64_t i = 0; i < CORRECTNESS_QUERIES; i++)
    {
      uint64_t pos = rand() % text.size();
      if(text[pos] != relative[pos]) { errors++; }
    }
    std::cout << "Completed " << CORRECTNESS_QUERIES << " access queries with " << errors << " error(s)." << std::endl;

    {
      std::vector<uint64_t> queries(TIMING_QUERIES);
      for(uint64_t i = 0; i < TIMING_QUERIES; i++) { queries[i] = rand() % text.size(); }
      double start = readTimer();
      uint64_t checksum = 0;
      for(uint64_t i = 0; i < TIMING_QUERIES; i++) { checksum += relative.rank(queries[i]); }
      double seconds = readTimer() - start;
      std::cout << "rank(): " << (seconds * TOTAL_TIME_TO_NANOSECS) << " ns/query (checksum " << checksum << ")" << std::endl;
    }

    {
      std::vector<uint64_t> queries(TIMING_QUERIES);
      for(uint64_t i = 0; i < TIMING_QUERIES; i++) { queries[i] = rand() % onebits + 1; }
      double start = readTimer();
      uint64_t checksum = 0;
      for(uint64_t i = 0; i < TIMING_QUERIES; i++) { checksum += relative.select_1(queries[i]); }
      double seconds = readTimer() - start;
      std::cout << "select_1(): " << (seconds * TOTAL_TIME_TO_NANOSECS) << " ns/query (checksum " << checksum << ")" << std::endl;
    }

    {
      std::vector<uint64_t> queries(TIMING_QUERIES);
      for(uint64_t i = 0; i < TIMING_QUERIES; i++) { queries[i] = rand() % (text.size() - onebits) + 1; }
      double start = readTimer();
      uint64_t checksum = 0;
      for(uint64_t i = 0; i < TIMING_QUERIES; i++) { checksum += relative.select_0(queries[i]); }
      double seconds = readTimer() - start;
      std::cout << "select_0(): " << (seconds * TOTAL_TIME_TO_NANOSECS) << " ns/query (checksum " << checksum << ")" << std::endl;
    }

    {
      std::vector<uint64_t> queries(TIMING_QUERIES);
      for(uint64_t i = 0; i < TIMING_QUERIES; i++) { queries[i] = rand() % text.size(); }
      double start = readTimer();
      uint64_t checksum = 0;
      for(uint64_t i = 0; i < TIMING_QUERIES; i++) { checksum += relative[queries[i]]; }
      double seconds = readTimer() - start;
      std::cout << "access(): " << (seconds * TOTAL_TIME_TO_NANOSECS) << " ns/query (checksum " << checksum << ")" << std::endl;
    }

    std::cout << std::endl;
  }
}

template<class csa_type>
void
printSize(const csa_type& csa, std::string name)
{
  uint64_t bytes = size_in_bytes(csa);
  uint64_t sample_bytes = size_in_bytes(csa.sa_sample) + size_in_bytes(csa.isa_sample);

  std::cout << name << " CSA: " << inMegabytes(bytes) << " MB (" << inMegabytes(bytes - sample_bytes) << 
    " MB without samples)" << std::endl;
}

void
testWT(int argc, char** argv)
{
  std::cout << "Wavelet tree tests" << std::endl;
  std::cout << std::endl;

  for(int arg = 1; arg < argc; arg++)
  {
    std::mt19937_64 rng(0xDEADBEEF);
    double prob = atof(argv[arg]);
    std::cout << "Mutation rate: " << prob << std::endl;

    int_vector<8> reference = generateReference(rng, WT_SIZE, WT_ALPHABET);
    int_vector<8> text = generateVariant(rng, reference, prob);
    std::cout << "Reference length " << reference.size() << ", text length " << text.size() << "." << std::endl;

    csa_wt<> reference_csa;
    construct_im(reference_csa, reference);
    printSize(reference_csa, "Reference");

    csa_wt<> plain_csa;
    construct_im(plain_csa, text);
    printSize(plain_csa, "Plain");

    csa_wt<wt_huff<rlz_vector> > rlz_csa;
    construct_im(rlz_csa, text);
    // Ugly hack is about to begin.
    rlz_vector& rlz = const_cast<rlz_vector&>(rlz_csa.wavelet_tree.bv);
    bit_vector::rank_1_type r_r(&(reference_csa.wavelet_tree.bv));
    bit_vector::select_1_type r_s1(&(reference_csa.wavelet_tree.bv));
    bit_vector::select_0_type r_s0(&(reference_csa.wavelet_tree.bv));
    rlz.compress(reference_csa.wavelet_tree.bv, r_r, r_s1, r_s0);
    printSize(rlz_csa, "RLZ");

    uint64_t errors = 0;
    for(uint64_t i = 0; i < CORRECTNESS_QUERIES; i++)
    {
      uint64_t pos = rng() % text.size();
      if(plain_csa[pos] != rlz_csa[pos]) { errors++; }
    }
    std::cout << "Completed " << CORRECTNESS_QUERIES << " locate queries with " << errors << " error(s)." << std::endl;

    std::cout << std::endl;
  }
}

void
compressBitvector(std::string name, const bit_vector& vec, const bit_vector& ref)
{
  std::cout << name; std::cout.flush();

  bit_vector::rank_1_type vec_rank(&vec);
  bit_vector::select_1_type vec_select_1(&vec);
  bit_vector::select_0_type vec_select_0(&vec);
  uint64_t plain_bytes = size_in_bytes(vec) + size_in_bytes(vec_rank) + size_in_bytes(vec_select_1) + size_in_bytes(vec_select_0);
  std::cout << "Plain " << inMegabytes(plain_bytes) << " MB"; std::cout.flush();

  rrr_vector<> rrr(vec);
  rrr_vector<>::rank_1_type rrr_rank(&rrr);
  rrr_vector<>::select_1_type rrr_select_1(&rrr);
  rrr_vector<>::select_0_type rrr_select_0(&rrr);
  uint64_t rrr_bytes = size_in_bytes(rrr) + size_in_bytes(rrr_rank) + size_in_bytes(rrr_select_1) + size_in_bytes(rrr_select_0);
  std::cout << ", RRR " << inMegabytes(rrr_bytes) << " MB"; std::cout.flush();

  sd_vector<> sd(vec);
  sd_vector<>::rank_1_type sd_rank(&sd);
  sd_vector<>::select_1_type sd_select_1(&sd);
  sd_vector<>::select_0_type sd_select_0(&sd);
  uint64_t sd_bytes = size_in_bytes(sd) + size_in_bytes(sd_rank) + size_in_bytes(sd_select_1) + size_in_bytes(sd_select_0);
  std::cout << ", SD " << inMegabytes(sd_bytes) << " MB"; std::cout.flush();

  bit_vector::rank_1_type ref_rank(&ref);
  bit_vector::select_1_type ref_select_1(&ref);
  bit_vector::select_0_type ref_select_0(&ref);
  rlz_vector rlz_vec(vec);
  rlz_vec.compress(ref, ref_rank, ref_select_1, ref_select_0);
  uint64_t rlz_bytes = size_in_bytes(rlz_vec);
  std::cout << ", RLZ " << inMegabytes(rlz_bytes) << " MB"; std::cout.flush();

  std::cout << std::endl;
}

typedef cst_sada<csa_wt<wt_huff<>, 32, 64, text_order_sa_sampling<bit_vector> > > cst_type;

void
testCST(int argc, char** argv)
{
  std::cout << "Compressed suffix tree tests" << std::endl;
  std::cout << std::endl;

  if(argc < 3)
  {
    std::cerr << "testCST(): Usage: " << argv[0] << " reference_file text_file" << std::endl;
    return;
  }

  std::cout << "Reference:        " << argv[1] << std::endl;
  cst_type reference_cst;
  construct(reference_cst, argv[1], 1);
  std::cout << "Reference CST:    " << inMegabytes(size_in_bytes(reference_cst)) << " MB" << std::endl;

  std::cout << "Text:             " << argv[2] << std::endl;
  cst_type text_cst;
  construct(text_cst, argv[2], 1);
  std::cout << "Text CST:         " << inMegabytes(size_in_bytes(reference_cst)) << " MB" << std::endl;

  std::cout << std::endl;

  compressBitvector("CSA               ", text_cst.csa.wavelet_tree.bv, reference_cst.csa.wavelet_tree.bv);
  compressBitvector("Sampled positions ", text_cst.csa.sa_sample.marked, reference_cst.csa.sa_sample.marked);
  compressBitvector("Tree              ", text_cst.bp, reference_cst.bp);

  std::ofstream out("temp.dat", std::ios_base::binary);
  text_cst.lcp.serialize(out); reference_cst.lcp.serialize(out);
  out.close();
  std::ifstream in("temp.dat", std::ios_base::binary);
  bit_vector text_lcp, ref_lcp;
  bit_vector::select_1_type sel;
  text_lcp.load(in); sel.load(in);
  ref_lcp.load(in);
  in.close();
  compressBitvector("LCP               ", text_lcp, ref_lcp);

  std::cout << std::endl;
}

//------------------------------------------------------------------------------

bit_vector
generateReference(std::mt19937_64& rng, uint64_t size)
{
  bit_vector reference(size);
  uint64_t* data = reference.data();
  for(uint64_t i = 0; i < reference.capacity() >> 6; i++) { data[i] = rng(); }
  return reference;
}

bit_vector
generateVariant(std::mt19937_64& rng, const bit_vector& reference, double mutation_rate)
{
  bit_vector text(reference.size() * (1.0 + 2 * mutation_rate) + EXTRA_SIZE);

  uint64_t j = 0;
  for(uint64_t i = 0; i < reference.size(); i++)
  {
    double mutation = probability(rng);
    if(mutation < mutation_rate)
    {
      mutation /= mutation_rate;
      if(mutation < 0.1)  // Insertion
      {
        do
        {
          text[j] = rng() & 1; j++;
        }
        while(probability(rng) < EXTENSION_PROB);
        text[j] = reference[i]; j++;
      }
      else if(mutation < 0.2) // Deletion.
      {
        while(i < reference.size() - 1 && probability(rng) < EXTENSION_PROB) { i++; }
      }
      else  // Mismatch.
      {
        text[j] = !reference[i]; j++;
      }
    }
    else  // Match.
    {
      text[j] = reference[i]; j++;
    }
  }

  text.resize(j);
  return text;
}

//------------------------------------------------------------------------------

inline uint64_t randomChar(std::mt19937_64& rng, uint8_t alphabet_size)
{
  return rng() % alphabet_size + 1;
}

int_vector<8>
generateReference(std::mt19937_64& rng, uint64_t size, uint8_t alphabet_size)
{
  if(alphabet_size == 0)
  {
    std::cerr << "generateReference(): Alphabet size cannot be 0!" << std::endl;
  }

  int_vector<8> reference(size);
  for(uint64_t i = 0; i < reference.size(); i++) { reference[i] = randomChar(rng, alphabet_size) + 1; }
  return reference;
}

int_vector<8>
generateVariant(std::mt19937_64& rng, const int_vector<8>& reference, double mutation_rate)
{
  uint8_t alphabet_size = *std::max_element(reference.begin(), reference.end());
  if(alphabet_size < 2)
  {
    std::cerr << "generateVariant(): Alphabet size must be at least 2!" << std::endl;
  }
  int_vector<8> text(reference.size() * (1.0 + 2 * mutation_rate) + EXTRA_SIZE);

  uint64_t j = 0;
  for(uint64_t i = 0; i < reference.size(); i++)
  {
    double mutation = probability(rng);
    if(mutation < mutation_rate)
    {
      mutation /= mutation_rate;
      if(mutation < 0.1)  // Insertion
      {
        do
        {
          text[j] = randomChar(rng, alphabet_size); j++;
        }
        while(probability(rng) < EXTENSION_PROB);
        text[j] = reference[i]; j++;
      }
      else if(mutation < 0.2) // Deletion.
      {
        while(i < reference.size() - 1 && probability(rng) < EXTENSION_PROB) { i++; }
      }
      else  // Mismatch.
      {
        uint64_t temp = (text[j] + randomChar(rng, alphabet_size - 1) - 1) % alphabet_size + 1;
        text[j] = temp; j++;
      }
    }
    else  // Match.
    {
      text[j] = reference[i]; j++;
    }
  }

  text.resize(j);
  return text;
}

//------------------------------------------------------------------------------
