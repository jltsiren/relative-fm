#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

#include <unistd.h>


uint64_t fileSize(std::ifstream& file);

const uint64_t MEGABYTE = 1048576;


int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: fasta2seq [-n] base_name" << std::endl;
    std::cerr << "  -n  Truncate runs of Ns into single Ns" << std::endl;
    return 1;
  }

  bool truncate_Ns = false;
  int c = 0;
  while((c = getopt(argc, argv, "n")) != -1)
  {
    if(c == 'n') { truncate_Ns = true; }
    else { return 2; }
  }

  std::string base_name = argv[optind];
  std::string fasta_name = base_name + ".fa";
  std::cout << "Extracting sequences from FASTA file " << fasta_name << std::endl;
  if(truncate_Ns)
  {
    std::cout << "Truncating runs of Ns" << std::endl;
  }

  std::ifstream fasta_file(fasta_name.c_str(), std::ios_base::binary);
  if(!fasta_file)
  {
    std::cerr << "fasta2seq: Cannot open FASTA file " << fasta_name << std::endl;
    return 2;
  }
  uint64_t bytes = fileSize(fasta_file);
  char* buffer = new char[bytes];
  fasta_file.read(buffer, bytes);
  fasta_file.close();
  std::cout << "Bytes: " << bytes << std::endl;
  std::cout << std::endl;

  std::ofstream seq_file(base_name.c_str(), std::ios_base::binary);
  if(!seq_file)
  {
    std::cerr << "fasta2seq: Cannot open sequence file " << base_name << std::endl;
    return 3;
  }

  char* write_buffer = new char[1048576];
  uint64_t sequences = 0, size = 0;
  std::vector<uint64_t> counts(256, 0);
  bool header_line = false, in_run = false;
  for(uint64_t i = 0; i < bytes; i++)
  {
    if(header_line)
    {
      if(buffer[i] == '\n') { header_line = false; }
      continue;
    }
    if(buffer[i] == '>') { header_line = true; sequences++; continue; }
    if(buffer[i] == '\n') { continue; }
    unsigned char temp = toupper(buffer[i]);
    if(truncate_Ns)
    {
      if(temp == 'N')
      {
        if(in_run) { continue; }
        in_run = true;
      }
      else { in_run = false; }
    }
    counts[temp]++;
    write_buffer[size % MEGABYTE] = temp;
    size++;
    if(size % MEGABYTE == 0)
    {
      seq_file.write(write_buffer, MEGABYTE);
    }
  }
  if(size % MEGABYTE != 0) { seq_file.write(write_buffer, size % MEGABYTE); }
  seq_file.close();
  delete[] write_buffer; write_buffer = 0;
  delete[] buffer; buffer = 0;

  std::cout << "Sequences: " << sequences << std::endl;
  std::cout << "Bytes: " << size << std::endl;
  for(uint64_t i = 0; i < 256; i++)
  {
    if(counts[i] > 0) { std::cout << "counts[" << i << "] = " << counts[i] << std::endl; }
  }
  std::cout << std::endl;

  return 0;
}

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
