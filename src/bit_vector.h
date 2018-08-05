#ifndef BIT_VECTOR_H_
#define BIT_VECTOR_H_

#include <cstddef>
#include <cstdint>

#include <deque>
#include <memory>
#include <vector>

namespace succient_bv {

class BitVector {
 public:
  BitVector(const std::deque<bool> &v) { Init(v); }

  BitVector(const std::vector<bool> &v) { Init(v); }

  ~BitVector() {}

  uint64_t Rank(uint64_t x) const;

  uint64_t  Select(uint64_t i) const;

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

  // aligned bit vector storing every 1/2 w bits.
  std::vector<uint32_t> b_;
  // store rank at i * w^2 in the bit vector. rank is at most w bits.
  std::vector<uint64_t> r1_;
  // store rank at i * 1/2 w in a w^2 bits block. rank is at most 2lg(w) bits.
  std::vector<uint16_t> r2_;
  std::vector<uint8_t> is_sparse_;
  std::vector<size_t> vector_indices_;
  std::vector<std::vector<uint64_t> > sparse_s_;
  std::vector<uint8_t> select_table_;
};

} // namespace succient_bv

#endif // BIT_VECTOR_H_
