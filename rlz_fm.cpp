/*
  Copyright (c) 2015 Genome Research Ltd.
  Copyright (c) 2014 Jouni Siren

  Author: Jouni Siren <jouni.siren@iki.fi>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#include "rlz_fm.h"

namespace relative
{

//------------------------------------------------------------------------------

const std::string RLZFM::EXTENSION = ".rlzfm";

RLZFM::RLZFM(const SimpleFM<>& ref, const SimpleFM<>& seq, const csa_wt<>* csa) :
  reference(ref), alpha(seq.alpha)
{
  this->block_rank = 0;

  // Extract the BWT.
  int_vector<8> seq_buffer(seq.bwt.size()); seq.extractBWT(seq_buffer);

  // Parse the BWT.
  std::vector<uint64_t> phrase_starts, phrase_lengths;
  if(csa != 0)
  {
    relativeLZ(seq_buffer, *csa, phrase_starts, phrase_lengths, this->mismatches);
  }
  else
  {
    int_vector<8> ref_buffer(ref.bwt.size()); ref.extractBWT(ref_buffer);
    relativeLZ(seq_buffer, ref_buffer, phrase_starts, phrase_lengths, this->mismatches);
  }
  util::clear(seq_buffer);

  // Initialize phrases and blocks.
  this->phrases.init(phrase_starts, phrase_lengths);
  util::assign(this->blocks, CumulativeArray(phrase_lengths));

  // Initialize rank structures.
  this->block_rank = new CumulativeArray[this->alpha.sigma];
  std::vector<uint64_t>* rank_buffers = new std::vector<uint64_t>[this->alpha.sigma];
  for(uint64_t phrase = 0; phrase < phrase_starts.size(); phrase++)
  {
    for(uint64_t c = 0; c < this->alpha.sigma; c++) // Determine rank information.
    {
      uint64_t val = this->countOf(phrase_starts[phrase], phrase_lengths[phrase] - 1, this->alpha.comp2char[c]);
      if(this->mismatches[phrase] == this->alpha.comp2char[c]) { val++; }
      if(phrase % BLOCK_SIZE == 0)  // Add a new block.
      {
        rank_buffers[c].push_back(val);
      }
      else
      {
        *(rank_buffers[c].rbegin()) += val;
      }
    }
  }
  for(uint64_t c = 0; c < this->alpha.sigma; c++)
  {
    util::assign(this->block_rank[c], CumulativeArray(rank_buffers[c]));
  }
  delete[] rank_buffers; rank_buffers = 0;
}

RLZFM::RLZFM(const SimpleFM<>& ref, const std::string& base_name) :
  reference(ref)
{
  this->block_rank = 0;
  std::string filename = base_name + EXTENSION;
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "RLZFM::RLZFM(): Cannot open input file " << filename << std::endl;
    return;
  }
  this->loadFrom(input);
  input.close();
}

RLZFM::RLZFM(const SimpleFM<>& ref, std::istream& input) :
  reference(ref)
{
  this->block_rank = 0;
  this->loadFrom(input);
}

RLZFM::~RLZFM()
{
  delete[] this->block_rank; this->block_rank = 0;
}

//------------------------------------------------------------------------------

uint64_t
RLZFM::reportSize(bool print) const
{
  uint64_t phrase_bytes = this->phrases.reportSize();
  uint64_t block_bytes = size_in_bytes(this->blocks);
  uint64_t mismatch_bytes = size_in_bytes(this->mismatches);

  uint64_t rlz_bytes = phrase_bytes + block_bytes + mismatch_bytes;
  uint64_t rank_bytes = 0;
  for(uint64_t c = 0; c < this->alpha.sigma; c++) { rank_bytes += size_in_bytes(this->block_rank[c]); }
  uint64_t bytes = size_in_bytes(this->alpha) + rlz_bytes + rank_bytes;

  if(print)
  {
#ifdef VERBOSE_OUTPUT
    printSize("Phrases", phrase_bytes, this->size());
    printSize("Blocks", block_bytes, this->size());
    printSize("Mismatches", mismatch_bytes, this->size());
#else
    printSize("RLZ", rlz_bytes, this->size());
#endif
    printSize("Rank", rank_bytes, this->size());
    printSize("RLZ FM", bytes, this->size());
    std::cout << std::endl;
  }

  return bytes;
}

void
RLZFM::writeTo(const std::string& base_name) const
{
  std::string filename = base_name + EXTENSION;
  std::ofstream output(filename.c_str(), std::ios_base::binary);
  if(!output)
  {
    std::cerr << "RLZFM::writeTo(): Cannot open output file " << filename << std::endl;
    return;
  }
  this->writeTo(output);
  output.close();
}

void
RLZFM::writeTo(std::ostream& output) const
{
  this->alpha.serialize(output);
  this->phrases.serialize(output);
  this->blocks.serialize(output);
  for(uint64_t c = 0; c < this->alpha.sigma; c++) { this->block_rank[c].serialize(output); }
  this->mismatches.serialize(output);
}

void
RLZFM::loadFrom(std::istream& input)
{
  this->alpha.load(input);
  this->phrases.load(input);
  this->blocks.load(input);
  delete[] this->block_rank; this->block_rank = new CumulativeArray[this->alpha.sigma];
  for(uint64_t c = 0; c < this->alpha.sigma; c++) { this->block_rank[c].load(input); }
  this->mismatches.load(input);
}

//------------------------------------------------------------------------------

uint64_t
RLZFM::rank(uint64_t i, uint8_t c) const
{
  if(!(hasChar(this->alpha, c))) { return 0; }
  uint64_t comp = this->alpha.char2comp[c];
  if(i >= this->size()) { return this->alpha.C[comp + 1] - this->alpha.C[comp]; }

  uint64_t phrase = this->blocks.inverse(i); // Phrase is 0-based.
  if(phrase == 0) { return this->countOf(this->phrases.decode(0, 0), i, c); }

  uint64_t base_phrase = (phrase / BLOCK_SIZE) * BLOCK_SIZE;
  uint64_t text_pos = 0;  // Starting position of the phrase in the text.
  uint64_t base_rank = 0; // Rank so far.
  if(phrase >= BLOCK_SIZE)
  {
    text_pos = this->blocks.sum(base_phrase);
    base_rank = this->block_rank[comp].sum(base_phrase / BLOCK_SIZE);
  }
  for(uint64_t cur = base_phrase; cur < phrase; cur++)
  {
    uint64_t block_size = this->blocks[cur];
    if(block_size > 1) { base_rank += this->countOf(this->phrases.decode(cur, text_pos), block_size - 1, c); }
    if(this->mismatches[cur] == c) { base_rank++; }
    text_pos += block_size;
  }

  return base_rank + this->countOf(this->phrases.decode(phrase, text_pos), i - text_pos, c);
}

//------------------------------------------------------------------------------

} // namespace relative
