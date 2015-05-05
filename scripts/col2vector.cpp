#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>

//------------------------------------------------------------------------------

const uint64_t MEGABYTE = 1048576;

void writeBytes(std::ostream& out, char* bytes, uint64_t n);

void buildAlphabet(uint64_t* counts, char* char2comp, const std::string& alpha_name);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 5)
  {
    std::cerr << "Usage: col2vector input column output alphabet" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  std::cout << "Converting a column into int_vector<8> with a packed alphabet" << std::endl;
  std::cout << std::endl;
  std::string input_name = argv[1];
  std::cout << "Input: " << input_name << std::endl;
  uint64_t col = std::stoul(argv[2]);
  std::cout << "Column: " << col << std::endl;
  std::string output_name = argv[3];
  std::cout << "Output: " << output_name << std::endl;
  std::string alpha_name = argv[4];
  std::cout << "Alphabet: " << alpha_name << std::endl;
  std::cout << std::endl;

  if(input_name == output_name)
  {
    std::cerr << "col2vector: Input and output files must be separate!" << std::endl;
    return 2;
  }

  std::ifstream in(input_name.c_str(), std::ios_base::binary);
  if(!in)
  {
    std::cerr << "col2vector: Cannot open input file " << input_name << std::endl;
    return 3;
  }
  std::ofstream out(output_name.c_str(), std::ios_base::binary);
  if(!out)
  {
    std::cerr << "col2vector: Cannot open output file " << output_name << std::endl;
    in.close();
    return 4;
  }

  // First pass: Build the alphabet.
  uint64_t bits = 0;
  uint64_t counts[257] = {};
  char     char2comp[256] = {};
  while(in)
  {
    std::string line;
    std::getline(in, line);
    if(line.length() <= col) { continue; }
    unsigned char c = line[col];
    if(isdigit(c)) { break; } // Ugly hack.
    counts[c]++; bits += 8;
  }
  std::cout << "Vector size: " << (bits / 8) << std::endl;
  out.write((char*)&bits, sizeof(bits));
  buildAlphabet(counts, char2comp, alpha_name);

  // Second pass: Write the BWT.
  char* buffer = new char[MEGABYTE];
  uint64_t bytes = 0;
  in.seekg(0);
  while(in)
  {
    std::string line;
    std::getline(in, line);
    if(line.length() <= col) { continue; }
    unsigned char c = line[col];
    if(isdigit(c)) { break; } // Ugly hack.
    buffer[bytes % MEGABYTE] = char2comp[c]; bytes++;
    if(bytes % MEGABYTE == 0) { writeBytes(out, buffer, MEGABYTE); }
  }
  if(bytes % MEGABYTE > 0) { writeBytes(out, buffer, bytes % MEGABYTE); }

  delete[] buffer; buffer = 0;
  in.close(); out.close();
  return 0;
}

//------------------------------------------------------------------------------

void
writeBytes(std::ostream& out, char* bytes, uint64_t n)
{
  uint64_t padding = 0;
  while((n + padding) % sizeof(uint64_t) != 0) { bytes[n + padding] = 0; padding++; }
  out.write(bytes, n + padding);
}

void
writeBuffer(std::ostream& out, char* bytes, uint64_t n)
{
  uint64_t bits = 8 * n;
  out.write((char*)&bits, sizeof(bits));
  writeBytes(out, bytes, n);
}

void
writeBuffer(std::ostream& out, uint64_t* data, uint64_t n)
{
  uint64_t bits = 64 * n;
  out.write((char*)&bits, sizeof(bits));
  out.write((char*)data, n * sizeof(uint64_t));
}

//------------------------------------------------------------------------------

void
buildAlphabet(uint64_t* counts, char* char2comp, const std::string& alpha_name)
{
  uint64_t sigma = 0;
  char     comp2char[256] = {};

  for(uint64_t i = 0; i < 256; i++)
  {
    if(counts[i] > 0)
    {
      char2comp[i] = sigma;
      comp2char[sigma] = i;
      counts[sigma] = counts[i];
      sigma++;
    }
  }

  for(uint64_t i = 0; i < sigma; i++)
  {
    std::cout << "  counts[" << i << "] = " << counts[i] << std::endl;
  }
  std::cout << std::endl;

  // Cumulative counts.
  for(uint64_t i = sigma; i > 0; i--) { counts[i] = counts[i - 1]; }
  counts[0] = 0;
  for(uint64_t i = 1; i <= sigma; i++) { counts[i] += counts[i - 1]; }

  std::ofstream out(alpha_name.c_str(), std::ios_base::binary);
  if(!out)
  {
    std::cerr << "col2vector: buildAlphabet(): Cannot open alphabet file " << alpha_name << std::endl;
    return;
  }
  writeBuffer(out, char2comp, 256);
  writeBuffer(out, comp2char, sigma);
  writeBuffer(out, counts, sigma + 1);
  out.write((char*)&sigma, sizeof(sigma));
  out.close();
}

//------------------------------------------------------------------------------
