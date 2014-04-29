#include <fstream>

#include <sdsl/construct.hpp>

#include "utils.h"


void build_bwt(bwt_type& bwt, alphabet_type& alpha, const std::string& filename);


int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: build_bwt input1 [input2 ...]" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  for(int i = 1; i < argc; i++)
  {
    bwt_type bwt;
    alphabet_type alpha;
    std::string base_name = argv[i];
    build_bwt(bwt, alpha, base_name);

    std::string bwt_name = base_name + ".bwt";
    std::ofstream output(bwt_name.c_str(), std::ios_base::binary);
    bwt.serialize(output);
    alpha.serialize(output);
    output.close();
  }

  return 0;
}


void
build_bwt(bwt_type& bwt, alphabet_type& alpha, const std::string& filename)
{
  std::cout << "File: " << filename << std::endl;

  cache_config config;

  {
    int_vector<8> text;
    load_vector_from_file(text, filename, 1);
    std::cout << "Text size: " << text.size() << std::endl;
    if(contains_no_zero_symbol(text, filename)) { append_zero_symbol(text); }
    store_to_cache(text, conf::KEY_TEXT, config);
  }

  {
    std::cout << "Building SA..." << std::endl;
    construct_sa<8>(config);
  }

  {
    std::cout << "Building BWT..." << std::endl;
    construct_bwt<8>(config);
    remove(cache_file_name(conf::KEY_TEXT, config));
    remove(cache_file_name(conf::KEY_SA, config));
  }

  {
    std::cout << "Building alphabet..." << std::endl;
    int_vector_buffer<8> bwt_buf(cache_file_name(conf::KEY_BWT, config));
    auto n = bwt_buf.size();
    alphabet_type temp(bwt_buf, n);
    alpha.swap(temp);
  }

  {
    std::cout << "Building WT..." << std::endl;
    int_vector_buffer<8> bwt_buf(cache_file_name(conf::KEY_BWT, config));
    auto n = bwt_buf.size();
    bwt_type temp(bwt_buf, n);
    bwt.swap(temp);
  }

  remove(cache_file_name(conf::KEY_BWT, config));
  std::cout << "BWT size: " << bwt.size() << std::endl;
  std::cout << std::endl;
}
