#include <cstdlib>

#include <sdsl/rrr_vector.hpp>

#include "rlz_vector.h"


const uint64_t SIZE = 6000 * (uint64_t)1048576;
const uint64_t EXTRA_SIZE = 1024; // Some extra space to avoid problems in text generation.
const uint64_t QUERIES = 100000;
const uint64_t TIMING_QUERIES = 100000000;
const uint64_t TOTAL_TIME_TO_NANOSECS = 10;
const double   EXTENSION_PROB = 0.3;

int
main(int argc, char** argv)
{
  std::cout << "Relative Lempel-Ziv bitvector test" << std::endl;
  std::cout << std::endl;

  for(int arg = 1; arg < argc; arg++)
  {
    double prob = atof(argv[arg]);
    srand(0xDEADBEEF);
    bit_vector reference(SIZE);
    bit_vector text(SIZE * (1.0 + 2 * prob) + EXTRA_SIZE);
    uint64_t j = 0;
    for(uint64_t i = 0; i < SIZE; i++)
    {
      reference[i] = rand() & 1;
      if(rand() / (RAND_MAX + 1.0) < prob)
      {
        uint64_t choice = rand() % 10;
        if(choice == 0) // Insertion.
        {
          do
          {
            text[j] = rand() & 1; j++;
          }
          while(rand() / (RAND_MAX + 1.0) < EXTENSION_PROB);
          text[j] = reference[i]; j++;
        }
        else if(choice == 1)  // Deletion.
        {
          while(i < SIZE - 1 && rand() / (RAND_MAX + 1.0) < EXTENSION_PROB) { i++; }
        }
        else  // Mismatch
        {
          text[j] = !reference[i]; j++;
        }
      }
      else
      {
        text[j] = reference[i]; j++;
      }
    }
    text.resize(j);
    uint64_t onebits = util::cnt_one_bits(text), errors = 0;
    std::cout << "Reference length " << SIZE << ", text length " << text.size() << ", mutation probability " << prob << "." << std::endl;

    std::vector<uint64_t> starts, lengths;
    bit_vector mismatches;
    double start_time = readTimer();
//    relativeLZ(text, reference, starts, lengths, mismatches);
    relativeLZSuccinct(text, reference, starts, lengths, mismatches);
    double seconds = readTimer() - start_time;
    std::cout << "Parsed the text as " << starts.size() << " phrases." << std::endl;
    std::cout << "Parsing took " << seconds << " seconds, " << inMegabytes(memoryUsage()) << " MB." << std::endl;

    uint64_t text_pos = 0;
    for(uint64_t phrase = 0; phrase < starts.size(); phrase++)
    {
      for(uint64_t ref_pos = starts[phrase]; ref_pos < starts[phrase] + lengths[phrase] - 1; text_pos++, ref_pos++)
      {
        if(text[text_pos] != reference[ref_pos]) { errors++; }
      }
      if(text[text_pos] != mismatches[phrase]) { errors++; } text_pos++;
    }
    std::cout << "Decompressed the text with " << errors << " error(s)." << std::endl;

    /*bit_vector::rank_1_type ref_rank; util::init_support(ref_rank, &reference);
    bit_vector::select_1_type ref_select; util::init_support(ref_select, &reference);
    bit_vector::rank_1_type text_rank; util::init_support(text_rank, &text);
    bit_vector::select_1_type text_select; util::init_support(text_select, &text);
    RLZVector relative(text, reference, ref_rank, ref_select);
    std::cout << "Built the relative bitvector: " << inBPC(relative.reportSize(), text.size()) << " vs. " <<
      inBPC(size_in_bytes(text) + size_in_bytes(text_rank) + size_in_bytes(text_select), text.size()) << " bpc." << std::endl;

    errors = 0;
    for(uint64_t i = 0; i < QUERIES; i++)
    {
      uint64_t pos = rand() % text.size();
      if(text_rank(pos) != relative.rank(pos)) { errors++; }
    }
    std::cout << "Completed " << QUERIES << " rank queries with " << errors << " error(s)." << std::endl;

    errors = 0;
    for(uint64_t i = 0; i < QUERIES; i++)
    {
      uint64_t pos = rand() % onebits + 1;
      if(text_select(pos) != relative.select(pos)) { errors++; }
    }
    std::cout << "Completed " << QUERIES << " select queries with " << errors << " error(s)." << std::endl;

    errors = 0;
    for(uint64_t i = 0; i < QUERIES; i++)
    {
      uint64_t pos = rand() % text.size();
      if(text[pos] != relative[pos]) { errors++; }
    }
    std::cout << "Completed " << QUERIES << " access queries with " << errors << " error(s)." << std::endl;

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
      for(uint64_t i = 0; i < TIMING_QUERIES; i++) { checksum += relative.select(queries[i]); }
      double seconds = readTimer() - start;
      std::cout << "select(): " << (seconds * TOTAL_TIME_TO_NANOSECS) << " ns/query (checksum " << checksum << ")" << std::endl;
    }

    {
      std::vector<uint64_t> queries(TIMING_QUERIES);
      for(uint64_t i = 0; i < TIMING_QUERIES; i++) { queries[i] = rand() % text.size(); }
      double start = readTimer();
      uint64_t checksum = 0;
      for(uint64_t i = 0; i < TIMING_QUERIES; i++) { checksum += relative[queries[i]]; }
      double seconds = readTimer() - start;
      std::cout << "access(): " << (seconds * TOTAL_TIME_TO_NANOSECS) << " ns/query (checksum " << checksum << ")" << std::endl;
    }*/

    std::cout << std::endl;
  }

  return 0;
}
