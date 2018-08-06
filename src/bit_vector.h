#ifndef BIT_VECTOR_H_
#define BIT_VECTOR_H_

#include <cstddef>
#include <cstdint>

#include <deque>
#include <memory>
#include <vector>
#include <iostream>

#include <x86intrin.h>

namespace succient_bv {

class BitVector {
 public:
  explicit BitVector(const std::deque<bool> &v) : b_(nullptr) { Init(v); }

  explicit BitVector(const std::vector<bool> &v) : b_(nullptr) { Init(v); }

  ~BitVector() { free(b_); }

  uint64_t Rank(uint64_t x) const {
    size_t r2_index = x / 32;
    uint32_t bits = b_[r2_index] >> (32 - 1 - (x % 32));

    // popcnt instruction is used instead of the pattern table for 1/2 w bits.
    return r1_[x / (64 * 64)] + r2_[r2_index] + _mm_popcnt_u32(bits);
  }

  uint64_t Select(uint64_t i) const {
    return s_[i / (64 * 64)]->Select(this, i % (64 * 64));
  }

  size_t n_bytes() const;

 private:
  template<class T> void Init(const T &v) {
    if (v.empty())
      throw std::runtime_error("Given container is empty.");

    InitVector(v);
    InitRankIndex();
    InitSelectIndex();
  }

  template<class T> void InitVector(const T &v);

  void InitRankIndex();

  void InitSelectIndex();

  void InitSelectTable();

  uint8_t SelectOn32bits(uint32_t bits, uint8_t i) const;

  class SelectIndex {
   public:
    virtual ~SelectIndex() = 0;

    virtual uint64_t Select(const BitVector *b, int16_t i) const = 0;

    virtual size_t n_bytes() const = 0;
  };

  class SelectIndexArray : public SelectIndex {
   public:
    SelectIndexArray(const BitVector *b, const std::vector<uint64_t> &s)
      : s_(s) {}

    ~SelectIndexArray() {}

    uint64_t Select(const BitVector *b, int16_t i) const override {
      return s_[i];
    }

    size_t n_bytes() const override {
      return s_.capacity() * sizeof(uint64_t);
    }

   private:
    std::vector<uint64_t> s_;
  };

  class SelectIndexTree : public SelectIndex {
   public:
    SelectIndexTree(const BitVector *b, const std::vector<uint64_t> &s)
      : cumsums_(nullptr) { Init(b, s); }

    ~SelectIndexTree() { if (cumsums_ != nullptr) free(cumsums_); }

    uint64_t Select(const BitVector *b, int16_t i) const override;

    size_t n_bytes() const override {
      return 8 * n_inner_ * sizeof(uint16_t);
    }

    void Dump() const;

   private:
    void Init(const BitVector *b, const std::vector<uint64_t> &s);

    int height_;
    unsigned int n_inner_;
    unsigned int first_block_index_;
    uint64_t first_block_offset_;
    /**
     This array stores cumsum of #ones in subtree for each node.
     Cumsum of #ones is at most w^2 and each node has sqrt(w) children.
     However, each node has 8 children in this implementation for efficiency.
     128 bits is needed per node because each node uses 8 int16_t.
     */
    int16_t *cumsums_;
  };

  uint64_t n_b_;
  // bit vector storing every 1/2 w bits.
  uint32_t *b_;
  // store rank at i * w^2 in the bit vector. rank is at most w bits.
  std::vector<uint64_t> r1_;
  // store rank at i * 1/2 w in a w^2 bits block. rank is at most 2lg(w) bits.
  std::vector<uint16_t> r2_;
  std::vector<std::shared_ptr<SelectIndex> > s_;
  std::vector<uint8_t> select_table_;
};

} // namespace succient_bv

#endif // BIT_VECTOR_H_
