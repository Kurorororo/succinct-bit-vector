#include "bit_vector.h"

#include <vector>

#include "gtest/gtest.h"

#include "naive_bit_vector.h"

namespace succient_bv {

class BitVectorTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    v1_.resize(8, false);
    v1_[0] = true;
    v1_[2] = true;
    v1_[3] = true;
    v1_[7] = true; // 10110001

    v2_.resize(10000, false);
    v2_[111] = true;
    v2_[831] = true;
    v2_[5215] = true;

    v3_.resize(1000000, false);
    n_true_ = 0;

    for (int i = 0; i < v3_.size(); ++i) {
      if (rand() % 2 == 0) {
        ++n_true_;
        v3_[i] = true;
      }
    }

    v4_.resize(1000000, false);
    n_sparse_true_ = 0;

    for (int i = 0; i < v4_.size(); ++i) {
      if (rand() % 1000 == 0) {
        ++n_sparse_true_;
        v4_[i] = true;
      }
    }
  }

  int n_true_;
  int n_sparse_true_;
  std::vector<bool> v1_;
  std::vector<bool> v2_;
  std::vector<bool> v3_;
  std::vector<bool> v4_;
};

TEST_F(BitVectorTest, RankWorks) {
  BitVector bv1(v1_);

  EXPECT_EQ(1u, bv1.Rank(0));
  EXPECT_EQ(1u, bv1.Rank(1));
  EXPECT_EQ(2u, bv1.Rank(2));
  EXPECT_EQ(3u, bv1.Rank(3));
  EXPECT_EQ(3u, bv1.Rank(4));
  EXPECT_EQ(3u, bv1.Rank(5));
  EXPECT_EQ(3u, bv1.Rank(6));
  EXPECT_EQ(4u, bv1.Rank(7));

  BitVector bv2(v2_);

  EXPECT_EQ(0u, bv2.Rank(110));
  EXPECT_EQ(1u, bv2.Rank(111));
  EXPECT_EQ(1u, bv2.Rank(830));
  EXPECT_EQ(2u, bv2.Rank(831));
  EXPECT_EQ(2u, bv2.Rank(5214));
  EXPECT_EQ(3u, bv2.Rank(5215));
  EXPECT_EQ(3u, bv2.Rank(9999));

  BitVector bv3(v3_);
  NaiveBitVector nbv3(v3_);

  for (int i = 0; i < v3_.size(); ++i)
    EXPECT_EQ(nbv3.Rank(i), bv3.Rank(i));

  BitVector bv4(v4_);
  NaiveBitVector nbv4(v4_);

  for (int i = 0; i < v4_.size(); ++i)
    EXPECT_EQ(nbv4.Rank(i), bv4.Rank(i));
}

TEST_F(BitVectorTest, SelectWorks) {
  BitVector bv1(v1_);

  EXPECT_EQ(0, bv1.Select(0));
  EXPECT_EQ(2, bv1.Select(1));
  EXPECT_EQ(3, bv1.Select(2));
  EXPECT_EQ(7, bv1.Select(3));

  BitVector bv2(v2_);

  EXPECT_EQ(111u, bv2.Select(0));
  EXPECT_EQ(831u, bv2.Select(1));
  EXPECT_EQ(5215u, bv2.Select(2));

  BitVector bv3(v3_);
  NaiveBitVector nbv3(v3_);

  for (int i = 0; i < n_true_; ++i)
    EXPECT_EQ(nbv3.Select(i), bv3.Select(i));

  BitVector bv4(v4_);
  NaiveBitVector nbv4(v4_);

  for (int i = 0; i < n_sparse_true_; ++i)
    EXPECT_EQ(nbv4.Select(i), bv4.Select(i));
}

} // namespace succient_bv
