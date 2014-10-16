#include "rlz_fm.h"

//------------------------------------------------------------------------------

const std::string RLZFM::EXTENSION = ".rlzfm";

RLZFM::RLZFM(const SimpleFM<>& ref, const SimpleFM<>& seq, const csa_wt<>* csa) :
  reference(ref), alpha(seq.alpha)
{
  this->block_rank = 0;

  // Extract the BWT.
  int_vector<8> seq_buffer(seq.bwt.size());
  for(uint64_t i = 0; i < seq_buffer.size(); i++) { seq_buffer[i] = seq.bwt[i]; }

  // Parse the BWT.
  std::vector<uint64_t> phrase_starts, phrase_lengths;
  if(csa != 0)
  {
    relativeLZ(seq_buffer, *csa, phrase_starts, phrase_lengths, this->mismatches);
  }
  else
  {
    int_vector<8> ref_buffer(ref.bwt.size());
    for(uint64_t i = 0; i < ref_buffer.size(); i++) { ref_buffer[i] = ref.bwt[i]; }
    relativeLZ(seq_buffer, ref_buffer, phrase_starts, phrase_lengths, this->mismatches);
  }
  util::clear(seq_buffer);

  // Initialize phrases and blocks.
  this->phrases.init(phrase_starts, phrase_lengths);
  this->blocks.init(phrase_lengths);

  // Initialize rank structures.
  this->block_rank = new rlz_helper[this->alpha.sigma - 1];
  std::vector<uint64_t>* rank_buffers = new std::vector<uint64_t>[this->alpha.sigma - 1];
  for(uint64_t phrase = 0; phrase < phrase_starts.size(); phrase++)
  {
    for(uint64_t c = 1; c < this->alpha.sigma; c++) // Determine rank information.
    {
      uint64_t val = this->countOf(phrase_starts[phrase], phrase_lengths[phrase] - 1, this->alpha.comp2char[c]);
      if(this->mismatches[phrase] == this->alpha.comp2char[c]) { val++; }
      if(phrase % BLOCK_SIZE == 0)  // Add a new block.
      {
        rank_buffers[c - 1].push_back(val);
      }
      else
      {
        *(rank_buffers[c - 1].rbegin()) += val;
      }
    }
  }
  for(uint64_t c = 1; c < this->alpha.sigma; c++) { this->block_rank[c - 1].init(rank_buffers[c - 1]); }
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
  uint64_t block_bytes = this->blocks.reportSize();
  uint64_t mismatch_bytes = size_in_bytes(this->mismatches);

  uint64_t rlz_bytes = phrase_bytes + block_bytes + mismatch_bytes;
  uint64_t rank_bytes = 0;
  for(uint64_t c = 1; c < this->alpha.sigma; c++) { rank_bytes += this->block_rank[c - 1].reportSize(); }
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
  for(uint64_t c = 1; c < this->alpha.sigma; c++) { this->block_rank[c - 1].serialize(output); }
  this->mismatches.serialize(output);
}

void
RLZFM::loadFrom(std::istream& input)
{
  this->alpha.load(input);
  this->phrases.load(input);
  this->blocks.load(input);
  delete[] this->block_rank; this->block_rank = new rlz_helper[this->alpha.sigma - 1];
  for(uint64_t c = 1; c < this->alpha.sigma; c++) { this->block_rank[c - 1].load(input); }
  this->mismatches.load(input);
}

//------------------------------------------------------------------------------

uint64_t
RLZFM::rank(uint64_t i, uint8_t c) const
{
  if(c == 0 || !(hasChar(this->alpha, c))) { return 0; }
  uint64_t comp = this->alpha.char2comp[c];
  if(i >= this->size()) { return this->alpha.C[comp + 1] - this->alpha.C[comp]; }

  uint64_t phrase = this->blocks.blockFor(i); // Phrase is 0-based.
  if(phrase == 0) { return this->countOf(this->phrases.decode(0, 0), i, c); }

  uint64_t base_phrase = (phrase / BLOCK_SIZE) * BLOCK_SIZE;
  uint64_t text_pos = 0;  // Starting position of the phrase in the text.
  uint64_t base_rank = 0; // Rank so far.
  if(phrase >= BLOCK_SIZE)
  {
    text_pos = this->blocks.itemsAfter(base_phrase - 1);
    base_rank = this->block_rank[comp - 1].itemsAfter(base_phrase / BLOCK_SIZE - 1);
  }
  for(uint64_t cur = base_phrase; cur < phrase; cur++)
  {
    uint64_t block_size = this->blocks.blockSize(cur);
    if(block_size > 1) { base_rank += this->countOf(this->phrases.decode(cur, text_pos), block_size - 1, c); }
    if(this->mismatches[cur] == c) { base_rank++; }
    text_pos += block_size;
  }

  return base_rank + this->countOf(this->phrases.decode(phrase, text_pos), i - text_pos, c);
}

//------------------------------------------------------------------------------
