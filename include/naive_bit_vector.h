#ifndef NAIVE_BIT_VECTOR_H_
#define NAIVE_BIT_VECTOR_H_

#include <cmath>

#include <deque>
#include <vector>

namespace succinct_bv {

class NaiveBitVector {
 public:
  explicit NaiveBitVector(const std::deque<bool> &v) { Init(v); }

  explicit NaiveBitVector(const std::vector<bool> &v) { Init(v); }

  ~NaiveBitVector() {}

  uint64_t Rank(uint64_t x) const { return rank_[x]; }

  uint64_t Select(uint64_t i) const { return select_[i]; }

  size_t n_bytes() const {
    return (rank_.capacity() + select_.capacity()) * sizeof(uint64_t);
  }

 private:
  template<class T> void Init(const T &v);

  std::vector<uint64_t> rank_;
  std::vector<uint64_t> select_;
};

} // namespace succinct_bv

#endif // NAIVE_BIT_VECTOR_H_
