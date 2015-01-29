#include <cstdlib>
#include <fstream>

#include <sys/resource.h>

#include "utils.h"

//------------------------------------------------------------------------------

const std::string BWT_EXTENSION = ".bwt";
const std::string NATIVE_BWT_EXTENSION = ".cbwt";
const std::string ALPHA_EXTENSION = ".alpha";
const std::string SIMPLE_FM_DEFAULT_ALPHABET("\0ACGNT", 6);

//------------------------------------------------------------------------------

void
printSize(const std::string& header, uint64_t bytes, uint64_t data_size, uint indent)
{
  std::string padding;
  if(header.length() + 1 < indent) { padding = std::string(indent - 1 - header.length(), ' '); }

  std::cout << header << ':' << padding << inMegabytes(bytes) << " MB (" << inBPC(bytes, data_size) << " bpc)" << std::endl;
}

void
printTime(const std::string& header, uint64_t found, uint64_t matches, uint64_t bytes, double seconds, uint indent)
{
  std::string padding;
  if(header.length() + 1 < indent) { padding = std::string(indent - 1 - header.length(), ' '); }

  std::cout << header << ':' << padding << "Found " << found << " patterns with " << matches << " occ in "
    << seconds << " seconds (" << (inMegabytes(bytes) / seconds) << " MB / s)" << std::endl;
}

//------------------------------------------------------------------------------

double
readTimer()
{
  return clock() / (double)CLOCKS_PER_SEC;
}

uint64_t
memoryUsage()
{
  rusage usage;
  getrusage(RUSAGE_SELF, &usage);
#ifdef RUSAGE_IN_BYTES
  return usage.ru_maxrss;
#else
  return ((uint64_t)1024) * usage.ru_maxrss;
#endif
}

//------------------------------------------------------------------------------

uint64_t
readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows)
{
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "readRows(): Cannot open input file " << filename << std::endl;
    return 0;
  }

  uint64_t chars = 0;
  while(input)
  {
    std::string buf;
    std::getline(input, buf);
    if(skip_empty_rows && buf.length() == 0) { continue; }
    rows.push_back(buf);
    chars += buf.length();
  }

  input.close();
  return chars;
}

//------------------------------------------------------------------------------

Alphabet::Alphabet() :
  char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma)
{
  this->m_sigma = 0;
}

Alphabet::Alphabet(const Alphabet& a) :
  char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma)
{
  this->copy(a);
}

Alphabet::Alphabet(Alphabet&& a) :
  char2comp(m_char2comp), comp2char(m_comp2char), C(m_C), sigma(m_sigma)
{
  *this = std::move(a);
}

Alphabet::~Alphabet()
{
}

void
Alphabet::copy(const Alphabet& a)
{
  this->m_char2comp = a.m_char2comp;
  this->m_comp2char = a.m_comp2char;
  this->m_C = a.m_C;
  this->m_sigma = a.m_sigma;
}

void
Alphabet::swap(Alphabet& a)
{
  if(this != &a)
  {
    this->m_char2comp.swap(a.m_char2comp);
    this->m_comp2char.swap(a.m_comp2char);
    this->m_C.swap(a.m_C);
    std::swap(this->m_sigma, a.m_sigma);
  }
}

Alphabet&
Alphabet::operator=(const Alphabet& a)
{
  if(this != &a) { this->copy(a); }
  return *this;
}

Alphabet&
Alphabet::operator=(Alphabet&& a)
{
  if(this != &a)
  {
    this->m_char2comp = std::move(a.m_char2comp);
    this->m_comp2char = std::move(a.m_comp2char);
    this->m_C = std::move(a.m_C);
    this->m_sigma = a.m_sigma;
  }
  return *this;
}

uint64_t
Alphabet::serialize(std::ostream& out, structure_tree_node* s, std::string name) const
{
  structure_tree_node* child = structure_tree::add_child(s, name, util::class_name(*this));
  uint64_t written_bytes = 0;
  written_bytes += this->m_char2comp.serialize(out, child, "char2comp");
  written_bytes += this->m_comp2char.serialize(out, child, "comp2char");
  written_bytes += this->m_C.serialize(out, child, "C");
  written_bytes += write_member(this->m_sigma, out, child, "sigma");
  structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
Alphabet::load(std::istream& in)
{
  this->m_char2comp.load(in);
  this->m_comp2char.load(in);
  this->m_C.load(in);
  read_member(this->m_sigma, in);
}

bool
Alphabet::assign(const std::string& alphabet_string)
{
  if(alphabet_string.length() != this->m_sigma)
  {
    std::cerr << "Alphabet::assign(): Alphabet string has length " << alphabet_string.length()
              << " (should be " << this->m_sigma << ")" << std::endl;
    return false;
  }

  for(size_type i = 1; i < alphabet_string.length(); i++)
  {
    if(alphabet_string[i] <= alphabet_string[i - 1])
    {
      std::cerr << "Alphabet::assign(): The alphabet string is not sorted" << std::endl;
      return false;
    }
  }

  util::assign(this->m_char2comp, int_vector<8>(MAX_SIGMA, 0));
  for(size_type i = 0; i < alphabet_string.length(); i++)
  {
    this->m_char2comp[(uint8_t)(alphabet_string[i])] = i;
    this->m_comp2char[i] = (uint8_t)(alphabet_string[i]);
  }
  return true;
}

//------------------------------------------------------------------------------
