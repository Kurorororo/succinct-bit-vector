#include "naive_bit_vector.h"

#include <vector>

#include "gtest/gtest.h"

namespace succient_bv {

class NaiveBitVectorTest : public ::testing::Test {
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
  }

  std::vector<bool> v1_;
  std::vector<bool> v2_;
};

TEST_F(NaiveBitVectorTest, RankWorks) {
  NaiveBitVector bv1(v1_);

  EXPECT_EQ(1u, bv1.Rank(0));
  EXPECT_EQ(1u, bv1.Rank(1));
  EXPECT_EQ(2u, bv1.Rank(2));
  EXPECT_EQ(3u, bv1.Rank(3));
  EXPECT_EQ(3u, bv1.Rank(4));
  EXPECT_EQ(3u, bv1.Rank(5));
  EXPECT_EQ(3u, bv1.Rank(6));
  EXPECT_EQ(4u, bv1.Rank(7));

  NaiveBitVector bv2(v2_);

  EXPECT_EQ(0u, bv2.Rank(110));
  EXPECT_EQ(1u, bv2.Rank(111));
  EXPECT_EQ(1u, bv2.Rank(830));
  EXPECT_EQ(2u, bv2.Rank(831));
  EXPECT_EQ(2u, bv2.Rank(5214));
  EXPECT_EQ(3u, bv2.Rank(5215));
  EXPECT_EQ(3u, bv2.Rank(9999));
}

TEST_F(NaiveBitVectorTest, SelectWorks) {
  NaiveBitVector bv1(v1_);

  EXPECT_EQ(0, bv1.Select(0));
  EXPECT_EQ(2, bv1.Select(1));
  EXPECT_EQ(3, bv1.Select(2));
  EXPECT_EQ(7, bv1.Select(3));

  NaiveBitVector bv2(v2_);

  EXPECT_EQ(111u, bv2.Select(0));
  EXPECT_EQ(831u, bv2.Select(1));
  EXPECT_EQ(5215u, bv2.Select(2));
}

} // namespace succient_bv
