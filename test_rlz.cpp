#include <cstdlib>

#include "rlz_vector.h"


const uint64_t SIZE = 128 * 1048576;
const uint64_t QUERIES = 100000;

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
    bit_vector text(SIZE);
    for(uint64_t i = 0; i < SIZE; i++)
    {
      reference[i] = rand() & 1;
      if(rand() / (RAND_MAX + 1.0) < prob) { text[i] = !reference[i]; }
      else { text[i] = reference[i]; }
    }
    uint64_t onebits = util::cnt_one_bits(text);
    std::cout << "Generated sequences of length " << SIZE << " with mutation probability " << prob << "." << std::endl;

    std::vector<range_type> phrases;
    relativeLZ(text, reference, phrases);
    std::cout << "Parsed the text as " << phrases.size() << " phrases." << std::endl;

    uint64_t pos = 0, errors = 0;
    for(auto phrase : phrases)
    {
      for(uint64_t i = phrase.first; i < phrase.first + phrase.second; i++, pos++)
      {
        if(text[pos] != reference[i]) { errors++; }
      }
    }
    std::cout << "Decompressed the text with " << errors << " error(s)." << std::endl;

    bit_vector::rank_1_type ref_rank; util::init_support(ref_rank, &reference);
    bit_vector::select_1_type ref_select; util::init_support(ref_select, &reference);
    bit_vector::rank_1_type text_rank; util::init_support(text_rank, &text);
    bit_vector::select_1_type text_select; util::init_support(text_select, &text);
    RLZVector relative(text, reference, ref_rank, ref_select);
    std::cout << "Built the relative bitvector; " <<inBPC(relative.reportSize(), text.size()) << " vs. " <<
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

    std::cout << std::endl;
  }

  return 0;
}
