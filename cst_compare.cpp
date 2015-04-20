#include <sdsl/suffix_trees.hpp>
#include <sdsl/suffix_tree_algorithm.hpp>

#include "relative_cst.h"

using namespace relative;

//------------------------------------------------------------------------------

//#define USE_HASH

template<class CST>
void buildCST(CST& cst, const std::string& base_name, const std::string& type);

template<class CST>
void maximalMatches(const CST& cst, const int_vector<8>& seq, const std::string& name, uint64_t indent = 18);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 4)
  {
    std::cerr << "Usage: cst_compare reference target sequence" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  std::cout << "Finding maximal matches with a CST" << std::endl;
  std::cout << std::endl;

  std::string ref_name = argv[1];
  std::cout << "Reference: " << ref_name << std::endl;

  SimpleFM<> ref_fm(ref_name);
  uint64_t fm_bytes = ref_fm.reportSize();
  printSize("FM-index", fm_bytes, ref_fm.size());
  RelativeLCP::lcp_type ref_lcp;
  load_from_file(ref_lcp, ref_name + LCP_EXTENSION);
  uint64_t lcp_bytes = size_in_bytes(ref_lcp);
  printSize("LCP array", lcp_bytes, ref_lcp.size());
  printSize("Reference data", fm_bytes + lcp_bytes, ref_fm.size());
  std::cout << std::endl;

  std::string target_name = argv[2];
  std::cout << "Target: " << target_name << std::endl;
  std::string seq_name = argv[3];
  std::cout << "Sequence: " << seq_name << std::endl;
  int_vector<8> seq;
  {
    std::ifstream in(seq_name.c_str(), std::ios_base::binary);
    if(!in)
    {
      std::cerr << "cst_compare: Cannot open sequence file " << seq_name << std::endl;
      return 2;
    }
    uint64_t size = util::file_size(seq_name); seq.resize(size);
    in.read((char*)(seq.data()), size); in.close();
  }
  std::cout << std::endl;

/*  {
    RelativeFM<> rfm(ref_fm, seq_name);
    RelativeLCP rlcp(ref_lcp, seq_name);
    RelativeCST<> rcst(rfm, rlcp);
    printSize("Relative CST", rcst.reportSize(), rcst.size());
#ifdef USE_HASH
    traverseHash(rcst, "Relative CST");
#else
    traverse(rcst, "Relative CST");
#endif
    std::cout << std::endl;
  }*/

/*    {
      std::string name = "cst_sct3_dac";
      cst_sct3<> cst;
      buildCST(cst, seq_name, name);
#ifdef USE_HASH
      traverseHash(cst, name);
#else
      traverse(cst, name);
#endif
      std::cout << std::endl;
    }

    {
      std::string name = "cst_sct3_plcp";
      cst_sct3<csa_wt<>, lcp_support_sada<>> cst;
      buildCST(cst, seq_name, name);
#ifdef USE_HASH
      traverseHash(cst, name);
#else
      traverse(cst, name);
#endif
      std::cout << std::endl;
    }*/

  {
    std::string name = "cst_sada";
    cst_sada<> cst;
    buildCST(cst, target_name, name);
    maximalMatches(cst, seq, name);
    std::cout << std::endl;
  }

  std::cout << "Memory used: " << inMegabytes(memoryUsage()) << " MB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

bool
file_exists(const std::string& name)
{
  std::ifstream in(name.c_str(), std::ios_base::binary);
  if(!in) { return false; }
  in.close();
  return true;
}

template<class CST>
void
buildCST(CST& cst, const std::string& base_name, const std::string& type)
{
  std::string cst_file = base_name + "." + type;
  if(file_exists(cst_file))
  {
    load_from_file(cst, cst_file);
  }
  else
  {
    construct(cst, base_name, 1);
    store_to_file(cst, cst_file);
  }
  printSize(type, size_in_bytes(cst), cst.size());
}

//------------------------------------------------------------------------------

template<class CST>
void
maximalMatches(const CST& cst, const int_vector<8>& seq, const std::string& name, uint64_t indent)
{
  double start = readTimer();
  uint64_t depth = 0, total_length = 0;
  typename CST::size_type char_pos = 0;
  typename CST::node_type curr = cst.root(), old;
  for(uint64_t i = 0; i < seq.size(); i++)
  {
    old = curr;
    while(forward_search(cst, curr, depth, seq[i], char_pos) == 0 && depth > 0)
    {
      std::cout << "Match at depth " << depth << std::endl;
      total_length += depth; depth--;
      // FIXME how to update char_pos after following suffix link?
      curr = old = cst.sl(old); char_pos = get_char_pos(curr.i, depth, cst);
    }
    if(curr != cst.root()) { depth++; }
  }
  double seconds = readTimer() - start;

  std::string padding;
  if(name.length() + 1 < indent) { padding = std::string(indent - 1 - name.length(), ' '); }
  std::cout << name << ":" << padding << "Average maximal match: " << (total_length / (double)(cst.size()))
                    << " (" << seconds << " seconds)" << std::endl;
}

//------------------------------------------------------------------------------
