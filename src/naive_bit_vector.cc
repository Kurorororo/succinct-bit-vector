#include "naive_bit_vector.h"

#include <cstdint>

namespace succient_bv {

template<class T>
void NaiveBitVector::Init(const T &v) {
  if (v.size() == 0)
    throw std::runtime_error("Given container is empty.");

  rank_.reserve(v.size());

  uint64_t count = 0;

  for (auto c : v) {
    if (c) ++count;
    rank_.push_back(count);
  }

  select_.reserve(count);

  for (uint64_t i = 0, n = v.size(); i < n; ++i)
    if (v[i]) select_.push_back(i);
}

template void NaiveBitVector::Init<std::deque<bool> >(
    const std::deque<bool> &v);

template void NaiveBitVector::Init<std::vector<bool> >(
    const std::vector<bool> &v);

} // namespace succient_bv
