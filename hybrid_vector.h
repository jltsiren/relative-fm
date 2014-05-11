#ifndef _RELATIVE_FM_HYBRID_VECTOR
#define _RELATIVE_FM_HYBRID_VECTOR


#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <sdsl/bits.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/iterators.hpp>
#include <sdsl/util.hpp>


namespace sdsl
{

//------------------------------------------------------------------------------

template<uint8_t t_b = 1, uint16_t sb_rate = 16> class rank_support_hybrid;
template<uint8_t t_b = 1, uint16_t sb_rate = 16> class select_support_hybrid;

template<uint16_t sb_rate = 16>
class hybrid_vector
{
public:
  typedef bit_vector::size_type size_type;
  typedef bit_vector::value_type value_type;
  typedef bit_vector::difference_type difference_type;
  typedef random_access_const_iterator<hybrid_vector> iterator;
  typedef rank_support_hybrid<1, sb_rate>   rank_1_type;
  typedef rank_support_hybrid<0, sb_rate>   rank_0_type;
  typedef select_support_hybrid<1, sb_rate> select_1_type;
  typedef select_support_hybrid<0, sb_rate> select_0_type;
  friend class rank_support_hybrid<1, sb_rate>;
  friend class rank_support_hybrid<0, sb_rate>;
  friend class select_support_hybrid<1, sb_rate>;
  friend class select_support_hybrid<0, sb_rate>;

  const static size_type k_word_bits = CHAR_BIT * sizeof(size_type);

  const static size_type k_block_size  = 255;
  const static size_type k_block_bytes = 32;
  const static size_type k_sb_rate     = sb_rate;
  const static size_type k_sb_header   = 2 * sizeof(size_type);
  const static size_type k_header_size = k_sb_header + 2 * k_sb_rate;
  const static size_type k_sblock_size = k_block_size * k_sb_rate;

  const static uint8_t k_plain_encoding = 0x00;
  const static uint8_t k_gap_encoding   = 0x40;
  const static uint8_t k_run0_encoding  = 0x80;
  const static uint8_t k_run1_encoding  = 0xC0;
  const static uint8_t k_b_size_mask    = 0x3F;
  const static uint8_t k_encoding_mask  = 0xC0;

  const static size_type k_uniform0_sb = ((size_type)1) << (k_word_bits - 2);
  const static size_type k_uniform1_sb = ((size_type)1) << (k_word_bits - 1);
  const static size_type k_sb_ptr_mask = k_uniform0_sb - 1;

//------------------------------------------------------------------------------

  hybrid_vector(const bit_vector& bv) :
    m_length(bv.size()), n_ones(0),
    n_block((m_length + k_block_size - 1) / k_block_size),
    n_sblock((n_block + k_sb_rate - 1) / k_sb_rate),
    trunk_size(0),
    m_trunk(0),
    m_header(new uint8_t[n_sblock * k_header_size])
  {
    std::fill(this->m_header, this->m_header + (this->n_sblock + 1) * k_header_size, 0);
    this->compute_trunk_size(bv);
    this->m_trunk = new uint8_t[this->trunk_size];

    size_type sb_ones = 0, trunk_ptr = 0;
    for(size_type beg = 0, block_id = 0; beg < this->m_length; beg += k_block_size, block_id++)
    {
      size_type end = std::min(beg + k_block_size, this->m_length);
      size_type sb_id = block_id / k_sb_rate;
      uint8_t* block_header =
        this->m_header + (sb_id * k_header_size + k_sb_header + (block_id % k_sb_rate) * 2));
      int_vector<1> buffer(k_block_size, 0);

      // Update superblock information.
      if((block_id % k_sb_rate) == 0)
      {
        size_type* sb_header = (size_type*)(this->m_header + sb_id * k_header_size);
        sb_header[0] = trunk_ptr; sb_header[1] = this->n_ones;
        if(sb_id > 0 && (sb_ones == 0 || sb_ones == k_sblock_size))
        {
          sb_header = (size_type*)(this->m_header + (sb_id - 1) * k_header_size);
          sb_header[0] |= (sb_ones == 0 ? k_uniform0_sb : k_uniform1_sb);
        }
        sb_ones = 0;
      }

      // Fill the buffer and count the statistics.
      size_type ones = 0, runs = 1;
      for(size_type i = 0; i < k_block_size; i++)
      {
        buffer[i - beg] = bv[i];
        if(bv[i] == 1) { ones++; }
        if(i + 1 < end && bv[i + 1] != bv[i]) { runs++; }
      }
      block_header[1] = ones;

      // Encode the blocks.
      size_type rle_size = (runs > 2 ? runs - 2 : 0);
      size_type gap_size = std::min(ones, k_block_size - ones);
      if(rle_size >= k_block_bytes && gap_size >= k_block_bytes)  // Use plain encoding.
      {
        size_type* block_data = (size_type*)(this->m_trunk + trunk_ptr);
        block_header[0] = k_block_bytes | k_plain_encoding;
        for(size_type i = 0; i < k_block_size; i += k_word_bits);
        {
          *block_data = 0;
          for(size_type j = 0; j < k_word_bits; j++)
          {
            if(i + j < k_block_size && buffer[i + j] == 1)
            {
              *block_data |= ((size_type)1) << j;
            }
          }
          ++block_data; trunk_ptr += sizeof(size_type);
        }
      }
      else if(gap_size > rle_size) // Use run-length encoding.
      {
        block_header[0] = rle_size | (buffer[0] == 1 ? k_run1_encoding : k_run0_encoding);
        for(size_type i = 0, run = 0; run < rle_size; i++)
        {
          if(buffer[i + 1] != buffer[i])
          {
            this->m_trunk[trunk_ptr] = i;
            run++; trunk_ptr++;
          }
        }
      }
      else  // Use gap encoding for the minority bits. Trunk size may be 0.
      {
        block_header[0] = (ones < k_block_size - ones ? ones : k_block_size - ones) | k_gap_encoding;
        value_type val = (ones < k_block_size - ones ? 1 : 0);
        for(size_type i = 0; i < k_block_size; i++)
        {
          if(buffer[i] == val) { this->m_trunk[trunk_ptr] = i; trunk_ptr++; }
        }
      }

      // Update global ranks.
      sb_ones += ones; this->n_ones += ones;
    }
  }

  hybrid_vector() :
    m_length(0), n_ones(0), n_block(0), n_sblock(0), trunk_size(0),
    m_trunk(0), m_header(0)
  {
  }

  hybrid_vector(const hybrid_vector& another)
  {
    this->copy(another);
  }

  hybrid_vector(hybrid_vector&& another)
  {
    this->swap(another);
  }

  hybrid_vector(std::istream& in) :
    m_length(0), n_ones(0), n_block(0), n_sblock(0), trunk_size(0),
    m_trunk(0), m_header(0)
  {
    this->load(in);
  }

  ~hybrid_vector()
  {
    this->clear();
  }

  hybrid_vector& operator=(const hybrid_vector& another)
  {
    this->copy(another);
    return *this;
  }

  hybrid_vector& operator=(hybrid_vector&& another)
  {
    this->swap(another);
    return *this;
  }

//------------------------------------------------------------------------------

  size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const
  {
    size_t bytes = sizeof(this->m_length) + sizeof(this->n_ones) + sizeof(this->n_block) + sizeof(this->n_sblock) + sizeof(this->trunk_size);
    out.write((char*)&(this->m_length), sizeof(this->m_length));
    out.write((char*)&(this->n_ones), sizeof(this->n_ones));
    out.write((char*)&(this->n_block), sizeof(this->n_block));
    out.write((char*)&(this->n_sblock), sizeof(this->n_sblock));
    out.write((char*)&(this->trunk_size), sizeof(this->trunk_size));

    if(this->m_trunk != 0)
    {
      bytes += this->trunk_size;
      out.write((char*)(this->m_trunk), this->trunk_size);
    }
    if(this->m_header != 0)
    {
      bytes += this->n_sblock * k_header_size;
      out.write((char*)(this->m_header), this->n_sblock * k_header_size);
    }

    return bytes;
  }

  void load(std::istream& in)
  {
    this->clear();

    in.read((char*)&(this->m_length), sizeof(this->m_length));
    in.read((char*)&(this->n_ones), sizeof(this->n_ones));
    in.read((char*)&(this->n_block), sizeof(this->n_block));
    in.read((char*)&(this->n_sblock), sizeof(this->n_sblock));
    in.read((char*)&(this->trunk_size), sizeof(this->trunk_size));

    this->m_trunk = new uint8_t[this->trunk_size];
    in.read((char*)(this->m_trunk), this->trunk_size);
    this->m_header = new uint8_t[this->n_sblock * k_header_size];
    in.read((char*)(this->m_header), this->n_sblock * k_header_size);
  }

  void copy(const hybrid_vector& another)
  {
    if(this == &another) { return; }
    this->clear();

    this->m_length = another.m_length;
    this->n_ones = another.n_ones;
    this->n_block = another.n_block;
    this->n_sblock = another.n_sblock;
    this->trunk_size = another.trunk_size;

    if(another.m_trunk != 0)
    {
      this->m_trunk = new uint8_t[this->trunk_size];
      std::copy(another.m_trunk, another.m_trunk + another.trunk_size, this->m_trunk);
    }
    if(another.m_header != 0)
    {
      this->m_header = new uint8_t[this->n_sblock * k_header_size];
      std::copy(another.m_header, another.m_header + another.n_sblock * k_header_size, this->m_header);
    }
  }

  void swap(hybrid_vector& another)
  {
    if(this == &another) { return; }
    std::swap(m_length, another.m_length);
    std::swap(n_ones, another.n_ones);
    std::swap(n_block, another.n_block);
    std::swap(n_sblock, another.n_sblock);
    std::swap(trunk_size, another.trunk_size);
    std::swap(m_trunk, another.m_trunk);
    std::swap(m_header, another.m_header);
  }

//------------------------------------------------------------------------------

  size_type size() const { return this->m_length; }
  size_type items() const { return this->n_ones; }
  iterator begin() const { return iterator(this, 0); }
  iterator end() const { return iterator(this, this->size()); }

  value_type operator[](size_type i) const
  {
    if(i < 0 || i >= this->size()) { return 0; }

    // Handle the superblock.
    size_type trunk_ptr = 0, header_ptr = 0, ones = 0, uniform = 0;
    this->vec->superblock_for_index(i, trunk_ptr, header_ptr, ones, uniform);
    if(uniform == k_uniform0_sb) { return 0; }  // A superblock of 0-bits.
    if(uniform == k_uniform1_sb) { return 1; }  // A superblock of 1-bits.

    // Find the correct block.
    size_type block_bytes = 0, block_ones = 0; uint8_t block_encoding = 0;
    this->vec->block_for_index(i, trunk_ptr, header_ptr, ones, block_bytes, block_ones, block_encoding);
    uint8_t* block_data = this->vec->m_trunk + trunk_ptr;

    // Decode the block.
    size_type offset = i % bit_vector_type::k_block_size;
    value_type val = 0;
    if(block_encoding == k_plain_encoding)
    {
      size_type* block_words = (size_type*)block_data; block_words += offset / k_word_bits;
      offset %= k_word_bits;
      return (*block_words) & (((size_type)1) << (k_word_bits - (offset + 1)));
    }
    else if(block_encoding == bit_vector_type::k_gap_encoding)
    {
      if(block_ones == bit_vector_type::k_block_size) { return 1; } // A block of 1-bits.
      else if(block_ones > 0)
      {
        bool mark_ones = (block_ones < k_block_size - block_ones);
        for(size_type i = 0; i < block_bytes; i++)
        {
          if(block_data[i] == offset) { return mark_ones; }
          if(block_data[i] > offset) { return !mark_ones; }
        }
        return !mark_ones;  // Are the remaining bits 1-bits?
      }
      else { return 0; }  // A block of 0-bits.
    }
    else  // Run-length encoding.
    {
      size_type run_start = 0, ones_found = 0;
      bool next_ones = (block_encoding == k_run1_encoding);
      for(size_type i = 0; i < block_bytes; i++)
      {
        if(block_data[i] >= offset) { return next_ones; }
        if(next_ones) { ones_found += block_data[i] + 1 - run_start; }
        run_start = block_data[i] + 1; next_ones = !next_ones;
      }
      // The last two runs are implicit.
      if(!next_ones) { return (offset >= k_block_size - (block_ones - ones_found)); }
      return (offset < run_start + block_ones - ones_found);
    }

    return 0;
  }

//------------------------------------------------------------------------------

private:

  void compute_trunk_size(const bit_vector& bv)
  {
    this->trunk_size = 0;

    for(size_type beg = 0; beg < this->m_length; beg += k_block_size)
    {
      size_type end = std::min(beg + k_block_size, this->m_length);
      size_type ones = 0, runs = 1;
      for(size_type i = beg; i < end; i++)
      {
        if(bv[i] == 1) { ones++; }
        if(i + 1 < end && bv[i + 1] != bv[i]) { runs++; }
      }

      size_type rle_size = (runs > 2 ? runs - 2 : 0);
      size_type gap_size = std::min(ones, k_block_size - ones);
      this->trunk_size += std::min(k_block_bytes, std::min(rle_size, gap_size));
    }
  }

  void clear()
  {
    delete[] this->m_trunk; this->m_trunk = 0;
    delete[] this->m_header; this->m_header = 0;
  }

//------------------------------------------------------------------------------

  void superblock_for_index(size_type i,
    size_type& trunk_ptr, size_type& header_ptr, size_type& ones, size_type& uniform) const
  {
    size_type superblock = (i / k_block_size) / k_sb_rate;
    size_type* sb_header = (size_type*)(this->m_header + superblock * k_header_size);

    trunk_ptr = sb_header[0] & k_sb_ptr_mask;
    header_ptr = superblock * k_header_size + k_sb_header;
    ones = sb_header[1];
    uniform = sb_header[0] ^ trunk_ptr;
  }

  void block_for_index(size_type i, size_type& trunk_ptr, size_type& header_ptr, size_type& ones,
    size_type& block_bytes, size_type&block_ones, uint8_t& block_encoding)
  {
    size_type block = i / k_block_size;
    for(size_type b = (block / k_sb_rate) * k_sb_rate; b != block; b++)
    {
      trunk_ptr += this->m_header[header_ptr] & k_b_size_mask;
      ones += this->m_header[header_ptr + 1];
      header_ptr += 2;
    }

    block_bytes = this->m_header[header_ptr] & k_b_size_mask;
    block_ones = this->m_header[header_ptr + 1];
    block_encoding = this->m_header[header_ptr] & k_encoding_mask;
  }

  void superblock_for_rank(size_type r, bool bit,
    size_type& trunk_ptr, size_type& header_ptr, size_type& bits, size_type& ones, size_type& uniform) const
  {
    size_type* sb_header = 0;
    size_type low = 0, high = this->n_sblock - 1;
    while(low < high) // Find the last sb with at most r - 1 correct bits before it.
    {
      size_type mid = low + (high - low + 1) / 2;
      sb_header = (size_type*)(this->m_header + mid * k_header_size);
      size_type count = (bit ? sb_header[1] : mid * k_sb_rate * k_block_size - sb_header[1]);
      if(count >= r) { high = mid - 1; }
      else           { low = mid; }
    }

    sb_header = (size_type*)(this->m_header + high * k_header_size);
    trunk_ptr = sb_header[0] & k_sb_ptr_mask;
    header_ptr = superblock * k_header_size + k_sb_header;
    bits = high * k_sb_rate * k_block_size;
    ones = sb_header[1];
    uniform = sb_header[0] ^ trunk_ptr;
  }

  void block_for_rank(size_type r, bool bit,
    size_type& trunk_ptr, size_type& header_ptr, size_type& bits, size_type& ones,
    size_type& block_bytes, size_type&block_ones, uint8_t& block_encoding)
  {
    // Find the first block after which there are >= r correct bits.
    while(true)
    {
      size_type count = ones + this->m_header[header_ptr + 1];
      if(!bit) { count = bits + k_block_size - count; }
      if(count >= r) { break; }
      trunk_ptr += this->m_header[header_ptr] & k_b_size_mask;
      bits += k_block_size;
      ones += this->m_header[header_ptr + 1];
      header_ptr += 2;
    }

    block_bytes = this->m_header[header_ptr] & k_b_size_mask;
    block_ones = this->m_header[header_ptr + 1];
    block_encoding = this->m_header[header_ptr] & k_encoding_mask;
  }

  uint8_t*  m_trunk;
  uint8_t*  m_header;
  size_type m_length;          // bitvector length (in bits)
  size_type n_block, n_sblock; // number of blocks
  size_type trunk_size;        // size of proper encoding
};

//------------------------------------------------------------------------------

/////////////////////////////// RANK ///////////////////////////////////////////

template<uint8_t t_b, uint16_t sb_rate>
class rank_support_hybrid
{
public:
  typedef hybrid_vector<sb_rate> bit_vector_type;
  typedef typename bit_vector_type::size_type size_type;

  explicit rank_support_hybrid(const bit_vector_type* v = nullptr)
  {
    this->set_vector(v);
  }

  const size_type size() const
  {
    if(this->vec == 0) { return 0; }
    return this->vec->size();
  }

  void load(std::istream&, const bit_vector_type* v = nullptr)
  {
    this->set_vector(v);
  }

  size_type serialize(std::ostream&, structure_tree_node* v = nullptr, std::string name = "") const
  {
    return 0;
  }
  
  rank_support_hybrid& operator=(const rank_support_hybrid& another)
  {
    this->set_vector(another.vec);
    return *this;
  }

  void swap(rank_support_hybrid& another)
  {
    std::swap(this->vec, another.vec);
  }

  void set_vector(const bit_vector_type* v = nullptr)
  {
    this->vec = v;
  }

//------------------------------------------------------------------------------

  const size_type operator()(size_type i) const { return this->rank(i); }

  const size_type rank(size_type i) const
  {
    if(i <= 0 || this->vec == 0) { return 0; }
    if(i >= this->size()) { i = this->size(); }
    i--;  // From now on, we can count the number of 1-bits/0-bits up to i.

    // Handle the superblock.
    size_type trunk_ptr = 0, header_ptr = 0, ones = 0, uniform = 0;
    this->vec->superblock_for_index(i, trunk_ptr, header_ptr, ones, uniform);
    if(uniform == bit_vector_type::k_uniform0_sb) // A superblock of 0-bits.
    {
      return this->real_rank(i, ones);
    }
    if(uniform == bit_vector_type::k_uniform1_sb) // A superblock of 1-bits.
    {
      ones += 1 + (i % (bit_vector_type::k_block_size * bit_vector_type::k_sb_rate));
      return this->real_rank(i, ones);
    }

    // Find the correct block.
    size_type block_bytes = 0, block_ones = 0; uint8_t block_encoding = 0;
    this->vec->block_for_index(i, trunk_ptr, header_ptr, ones, block_bytes, block_ones, block_encoding);
    uint8_t* block_data = this->vec->m_trunk + trunk_ptr;

    // Decode the block.
    size_type offset = i % bit_vector_type::k_block_size;
    if(block_encoding == bit_vector_type::k_plain_encoding)
    {
      size_type* block_words = (size_type*)block_data;
      while(offset >= bit_vector_type::k_word_bits)
      {
        ones += bits::cnt(*block_words);
        ++block_words; offset -= bit_vector_type::k_word_bits;
      }
      ones += bits::cnt(*block_words & bits::lo_set(offset + 1));
    }
    else if(block_encoding == bit_vector_type::k_gap_encoding)
    {
      if(block_ones == bit_vector_type::k_block_size) { ones += offset + 1; } // A block of 1-bits.
      else if(block_ones > 0)
      {
        size_type bits = 0;
        for(size_type i = 0; i < block_bytes; i++)
        {
          if(block_data[i] <= offset) { bits++; }
          if(block_data[i] >= offset) { break; }
        }
        ones += (block_ones < bit_vector_type::k_block_size - block_ones ? bits : offset + 1 - bits);
      }
    }
    else  // Run-length encoding.
    {
      size_type run_start = 0, ones_found = 0;
      bool next_ones = (block_encoding == bit_vector_type::k_run1_encoding);
      for(size_type i = 0; i < block_bytes; i++)
      {
        size_type end = std::min(block_data[i], offset);
        if(next_ones) { ones_found += end + 1 - run_start; }
        run_start = end + 1; next_ones = !next_ones;
        if(end >= offset) { break; }
      }
      if(run_start <= offset) // The last two runs are implicit.
      {
        if(!next_ones) { run_start = bit_vector_type::k_block_size - (block_ones - ones_found); }
        if(run_start <= offset) { ones_found += std::min(block_ones - ones_found, offset + 1 - run_start); }
      }
      ones += ones_found;
    }

    return this->real_rank(i, ones);
  }

private:
  const bit_vector_type* vec;

  size_type real_rank(size_type i, size_type ones) const
  {
    return (t_b == 0 ? i + 1 - ones : ones);
  }
};

//------------------------------------------------------------------------------

///////////////////////////////// SELECT ///////////////////////////////////////

template<uint8_t t_b, uint16_t sb_rate>
class select_support_hybrid
{
public:
  typedef hybrid_vector<sb_rate> bit_vector_type;
  typedef typename bit_vector_type::size_type size_type;

  explicit select_support_hybrid(const bit_vector_type* v = nullptr)
  {
    this->set_vector(v);
  }

  const size_type size() const
  {
    if(this->vec == 0) { return 0; }
    return this->vec->size();
  }

  void load(std::istream&, const bit_vector_type* v = nullptr)
  {
    this->set_vector(v);
  }

  size_type serialize(std::ostream&, structure_tree_node* v = nullptr, std::string name = "") const
  {
    return 0;
  }
  
  select_support_hybrid& operator=(const select_support_hybrid& another)
  {
    this->set_vector(another.vec);
    return *this;
  }

  void swap(select_support_hybrid& another)
  {
    std::swap(this->vec, another.vec);
  }

  void set_vector(const bit_vector_type* v = nullptr)
  {
    this->vec = v;
  }

//------------------------------------------------------------------------------

  const size_type operator()(size_type i) const { return this->select(i); }

  const size_type select(size_type r) const
  {
    // We return the position of r'th 1-bit/0-bit.
    if(r <= 0 || this->vec == 0) { return 0; }
    if(t_b == 1 && r > this->vec->items()) { return this->vec->size(); }
    if(t_b == 0 && r > this->vec->size() - this->vec->items()) { return this->vec->size(); }

    // Handle the superblock.
    size_type trunk_ptr = 0, header_ptr = 0, bits = 0, ones = 0, uniform = 0;
    this->vec->superblock_for_rank(r, t_b, trunk_ptr, header_ptr, bits, ones, uniform);
    if(uniform != 0) // A superblock of t_b.
    {
      if(t_b) { return bits + r - ones - 1; }
      else    { return r + ones - 1; }
    }

    // Find the correct block.
    size_type block_bytes = 0, block_ones = 0; uint8_t block_encoding = 0;
    this->vec->block_for_rank(r, t_b, trunk_ptr, header_ptr, ones, block_bytes, block_ones, block_encoding);
    uint8_t* block_data = this->vec->m_trunk + trunk_ptr;

    // Decode the block.
    if(block_encoding == bit_vector_type::k_plain_encoding)
    {
      size_type* block_words = (size_type*)block_data;
      while(true)
      {
        size_type temp = ones + bits::cnt(*block_words);
        size_type count = (t_b ? temp : bits + k_word_bits - temp); 
        if(count >= r) { break; }
        ++block_words; ones = temp; bits += k_word_bits;
      }
      return bits + (t_b ? bits::sel(*block_words, r - ones) : bits::sel(~(*block_words), r - (bits - ones)));
    }
    else if(block_encoding == bit_vector_type::k_gap_encoding)
    {
      // FIXME implement
    }
    else  // Run-length encoding.
    {
      // FIXME implement
    }

    return 0;
  }

private:
  const bit_vector_type* vec;
};

//------------------------------------------------------------------------------

} // namespace sdsl

#endif // _RELATIVE_FM_HYBRID_VECTOR
