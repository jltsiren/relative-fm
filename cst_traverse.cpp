#include <sdsl/suffix_trees.hpp>

#include "relative_cst.h"

using namespace relative;

//------------------------------------------------------------------------------

//#define USE_HASH

template<class CST>
void buildCST(CST& cst, const std::string& base_name, const std::string& type);

template<class CST>
void traverseHash(CST& cst, const std::string& name, uint64_t indent = 18);

template<class CST>
void traverse(CST& cst, const std::string& name, uint64_t indent = 18);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 3)
  {
    std::cerr << "Usage: cst_traverse ref seq1 [seq2 ...]" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  std::cout << "DFS traversal in compressed suffix trees" << std::endl;
  std::cout << std::endl;

  std::string ref_name = argv[1];
  std::cout << "Reference: " << ref_name << std::endl;
  std::cout << std::endl;

  SimpleFM<> ref_fm(ref_name);
  uint64_t fm_bytes = ref_fm.reportSize();
  printSize("FM-index", fm_bytes, ref_fm.size());
  RelativeLCP::lcp_type ref_lcp;
  load_from_file(ref_lcp, ref_name + LCP_EXTENSION);
  uint64_t lcp_bytes = size_in_bytes(ref_lcp);
  printSize("LCP array", lcp_bytes, ref_lcp.size());
  printSize("Reference data", fm_bytes + lcp_bytes, ref_fm.size());
  std::cout << std::endl;

  std::cout << std::endl;

  for(int arg = 2; arg < argc; arg++)
  {
    std::string seq_name = argv[arg];
    std::cout << "Sequence: " << seq_name << std::endl;
    std::cout << std::endl;

    {
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
    }

    {
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
    }

    {
      std::string name = "cst_sada";
      cst_sada<> cst;
      buildCST(cst, seq_name, name);
      traverse(cst, name);
      std::cout << std::endl;
    }

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

inline uint64_t
hashNode(const rcst_node& node, uint64_t hash)
{
  hash = fnv1a_hash(node.sp, hash);
  return fnv1a_hash(node.ep, hash);
}

inline uint64_t
hashNode(const bp_interval<uint64_t>& node, uint64_t hash)
{
  hash = fnv1a_hash(node.i, hash);
  return fnv1a_hash(node.j, hash);
}

template<class CST>
void
traverseHash(CST& cst, const std::string& name, uint64_t indent)
{
  double start = readTimer();
  uint64_t nodes = 0, hash = FNV_OFFSET_BASIS;
  for(auto iter = cst.begin(); iter != cst.end(); ++iter)
  {
    if(iter.visit() == 1) { nodes++; }
    hash = hashNode(*iter, hash);
  }
  double seconds = readTimer() - start;

  std::string padding;
  if(name.length() + 1 < indent) { padding = std::string(indent - 1 - name.length(), ' '); }
  std::cout << name << ":" << padding << nodes << " nodes in " << seconds
            << " seconds (hash " << hash << ")" << std::endl;
}

template<class CST>
void
traverse(CST& cst, const std::string& name, uint64_t indent)
{
  double start = readTimer();
  uint64_t nodes = 0;
  for(auto iter = cst.begin(); iter != cst.end(); ++iter)
  {
    if(iter.visit() == 1) { nodes++; }
  }
  double seconds = readTimer() - start;

  std::string padding;
  if(name.length() + 1 < indent) { padding = std::string(indent - 1 - name.length(), ' '); }
  std::cout << name << ":" << padding << nodes << " nodes in " << seconds << " seconds" << std::endl;
}

//------------------------------------------------------------------------------
