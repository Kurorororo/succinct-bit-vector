# Succient Bit Vector

Succient bit vector supports following queries in constant time using n + O(nlglgn/lgn) space.

- B: 0,1 vector of length n B[0]B[1]...B[n-1]
- rank(x)
  - number of ones in B[0..x]=B[0]B[1]...B[x].
- select(i)
  - position of i-th 1 from the head.

## Usage
Please build library, include `bit_vector.h` and link library.

`BitVector` can be instantiated from only `deque<bool>` and `vector<bool>`.

`uint64_t Rank(uint64_t x)` and `uint64_t Select(uint64_t i)` are supported.

### Example
```c++
#include <vector>

#include <bit_vector.h>

int main() {
  vector<bool> v(100, true)
  v[10] = true;
  v[50] = true;
  BitVector bv(v);
  // r = 1
  int r = bv.Rank(30);
  // s = 50
  int s = bv.Select(1);
}
```

```
$ g++ -std=c++17 -I./succient_bv/include -L./succient_bv/lib -lsuccient_bv example.cc
```

## Build
SSE4.2 and POPCNT is needed.

```
$ mkdir build
$ cd build
$ cmake ..
$ make
$ ctest
```

## References
R. Raman, V. Raman, and S. S. Rao. Succinct Indexable Dictionaries with Applications to Encoding k-ary Trees and Multisets, ACM Transactions on Algorithms (TALG) , Vol. 3, Issue 4, 2007.
