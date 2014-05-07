#include <fstream>
#include <iostream>
#include <limits>
#include <vector>


typedef std::vector<uint8_t> byte_vector;

byte_vector readFile(const std::string& filename);
int greedyLCS(byte_vector& seq1, byte_vector& seq2);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 3)
  {
    std::cerr << "Usage: lcs seq1 seq2" << std::endl;
    std::cerr << std::endl;
    return 1;
  }

  std::cout << "Input: " << argv[1] << " and " << argv[2] << std::endl;
  byte_vector seq1 = readFile(argv[1]);
  byte_vector seq2 = readFile(argv[2]);
  std::cout << "Lengths: " << seq1.size() << " and " << seq2.size() << std::endl;
  
  int lcs = greedyLCS(seq1, seq2);
  std::cout << "LCS length: " << lcs << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

uint64_t
fileSize(std::ifstream& file)
{
  std::streamoff curr = file.tellg();

  file.seekg(0, std::ios::end);
  std::streamoff size = file.tellg();
  file.seekg(0, std::ios::beg);
  size -= file.tellg();

  file.seekg(curr, std::ios::beg);
  return size;
}

byte_vector
readFile(const std::string& filename)
{
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "readFile(): Cannot open input file " << filename << "!" << std::endl;
    return byte_vector();
  }

  uint64_t size = fileSize(input);
  char* buffer = new char[size];
  input.read(buffer, size);
  input.close();

  byte_vector res(buffer, buffer + size);
  delete[] buffer; buffer = 0;
  return res;
}

//------------------------------------------------------------------------------

inline uint
mapToUint(int val)
{
  if(val >= 0) { return 2 * val; }
  else         { return 2 * (-val) - 1; }
}

/*
  Eugene W. Myers: An O(ND) Difference Algorithm and Its Variations. Algorithmica, 1986.

  The implementation assumes that offsets fit into int.
*/
int
greedyLCS(byte_vector& seq1, byte_vector& seq2)
{
  if(seq1.size() == 0 || seq1.size() == 0) { return 0; }

  // v[k] stores how many characters of ref have been processed on diagonal k.
  std::vector<int> v(3, 0); // v[mapToUint(1)] = 0;
  int seq1_len = seq1.size(), seq2_len = seq2.size();
  for(int d = 0; true; d++)
  {
    for(int k = -d; k <= d; k += 2)
    {
      int x = 0;
      if(k == -d || (k != d && v[mapToUint(k - 1)] < v[mapToUint(k + 1)]))
      {
        x = v[mapToUint(k + 1)];
      }
      else
      {
        x = v[mapToUint(k - 1)] + 1;
      }
      int y = x - k;
      while(x < seq1_len && y < seq2_len && seq1[x] == seq2[y]) { x++; y++; }
      v[mapToUint(k)] = x;
      if(x >= seq1_len && y >= seq2_len) { return (seq1_len + seq2_len - d) / 2; }
    }
    v.resize(v.size() + 2, 0);
  }

  return 0;
}

//------------------------------------------------------------------------------
