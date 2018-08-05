#include "bit_vector.h"

#include <cmath>

#include <iostream>

#include <x86intrin.h>

namespace succient_bv {

using std::vector;

template<class T>
void BitVector::InitVector(const T &v) {
  uint64_t n = v.size();
  b_.resize((n + 32 - 1) / 32, 0);

  for (uint64_t i = 0; i < n; ++i)
    if (v[i]) b_[i / 32] |= 1u << (32 - 1 - (i % 32));
}

void BitVector::InitRankIndex() {
  // the number of w^2 bits blocks is [n/w^2]+1.
  // every w^2 bits block contains 2w small (1/2 w bits) blocks.
  r1_.reserve(b_.size() / (2 * 64) + 1);
  // the number of 1/2 w bits blocks is [2n/w]+1.
  r2_.reserve(b_.size());

  uint64_t r1_sum = 0;
  uint64_t r2_sum = 0;

  for (uint64_t i = 0, m = b_.size(); i < m; ++i) {
    if (i % (2 * 64) == 0) {
      r1_.push_back(r1_sum);
      r2_sum = 0;
    }

    r2_.push_back(r2_sum);
    uint64_t count = _mm_popcnt_u32(b_[i]);
    r1_sum += count;
    r2_sum += count;
  }
}

uint64_t BitVector::Rank(uint64_t x) const {
  size_t r2_index = x / 32;
  uint32_t bits = b_[r2_index] >> (32 - 1 - (x % 32));

  // popcnt instruction is used instead of the pattern table for 1/2 w bits.
  return r1_[x / (64 * 64)] + r2_[r2_index] + _mm_popcnt_u32(bits);
}

void BitVector::InitSelectIndex() {
  InitSelectTable();

  vector<uint64_t> s;
  vector<uint64_t> next_s;

  for (uint64_t i = 0, m = b_.size(); i < m; ++i) {
    uint8_t count = _mm_popcnt_u32(b_[i]);

    for (uint8_t j = 0; j < count; ++j)
      s.push_back(i * 32 + static_cast<uint64_t>(SelectOn32bits(b_[i], j)));

    if (s.size() > (64 * 64)) {
      for (size_t j = 64 * 64; j < s.size(); ++j)
        next_s.push_back(s[j]);

      s.resize(64 * 64);
    }

    // a block contains w^2 ones.
    if ((s.size() == (64 * 64)) || (i == m - 1)) {

      // a block is sparse if the size of block > w^4 bits.
      if ((s.back() - s.front() + 1) > 64 * 64 * 64 * 64) {
        is_sparse_.push_back(1);
        vector_indices_.push_back(sparse_s_.size());
        sparse_s_.push_back(s);
      } else {
        is_sparse_.push_back(1);
        vector_indices_.push_back(sparse_s_.size());
        sparse_s_.push_back(s);
      }

      s.clear();
      s.swap(next_s);
    }
  }
}

void BitVector::InitSelectTable() {
  select_table_.resize(8 * 256, 8);

  for (size_t i = 0; i < 256; ++i) {
    uint16_t pattern = static_cast<uint16_t>(i);
    uint16_t index = 0;

    for (uint8_t j = 0; j < 8; ++j) {
      if (pattern & (1u << (8 - 1 - j))) {
        select_table_[(index << 8) + pattern] = j;
        ++index;
      }
    }
  }
}

uint8_t BitVector::SelectOn32bits(uint32_t bits, uint8_t i) const {
  uint8_t j = 0;

  while (j < 4) {
    uint8_t count = _mm_popcnt_u32((bits >> ((3 - j) * 8)) & 0xffU);
    if (i < count) break;
    ++j;
    i -= count;
  }

  return j * 8 + select_table_[(i << 8) + ((bits >> ((3 - j) * 8)) & 0xffU)];
}

uint64_t BitVector::Select(uint64_t i) const {
  uint64_t block_index = i / (64 * 64);

  if (is_sparse_[block_index])
    return sparse_s_[vector_indices_[block_index]][i % (64 * 64)];

  return sparse_s_[vector_indices_[block_index]][i % (64 * 64)];
}

size_t BitVector::n_bytes() const {
  size_t n = b_.capacity() * sizeof(uint32_t);
  n += r1_.capacity() * sizeof(uint64_t) + r2_.capacity() * sizeof(uint16_t);
  n += is_sparse_.capacity() * sizeof(uint8_t);
  n += vector_indices_.capacity() * sizeof(size_t);

  for (auto &v : sparse_s_)
    n += v.capacity() * sizeof(uint64_t);

  n += select_table_.capacity() * sizeof(uint8_t);

  return n;
}

template void BitVector::InitVector<std::deque<bool> >(
    const std::deque<bool> &v);

template void BitVector::InitVector<std::vector<bool> >(
    const std::vector<bool> &v);

} // namespace succient_bv
