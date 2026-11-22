// Comparison operators across sign permutations and magnitudes.

#include "unit_test_framework.h"

#include "biginteger/BigInteger.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Comparator.h"
#include "biginteger/ops/Comparison.h"

using namespace BigMath;

static void CheckEquality(const BigInteger &a, const BigInteger &b, bool expected)
{
  ASSERT_EQ(a.CompareTo(b) == 0, expected);
  ASSERT_EQ(a == b, expected);
  ASSERT_EQ(b == a, expected);
  ASSERT_EQ(a != b, !expected);
  ASSERT_EQ(b != a, !expected);
}

REGISTER_TEST(Compare, EqualityAcrossLimbPositions)
{
  for (SizeT size : {SizeT{1}, SizeT{2}, SizeT{3}, BigInteger::InlineLimbs,
                     BigInteger::InlineLimbs + 1, SizeT{8}, SizeT{9}, SizeT{17}, SizeT{64}})
  {
    std::vector<DataT> limbs(size, 0x12345678);
    for (bool negative : {false, true})
    {
      const BigInteger a(limbs, negative);
      CheckEquality(a, a, true);
      CheckEquality(a, BigInteger(limbs, negative), true);
      CheckEquality(a, BigInteger(limbs, !negative), false);
      for (SizeT i = 0; i < size; ++i)
      {
        limbs[i] ^= DataT{1} << (LimbBits - 1);
        CheckEquality(a, BigInteger(limbs, negative), false);
        limbs[i] ^= DataT{1} << (LimbBits - 1);
      }
      auto longer = limbs;
      longer.push_back(1);
      CheckEquality(a, BigInteger(longer, negative), false);
    }
  }
}

REGISTER_TEST(Compare, EqualityWithPaddedZeros)
{
  const BigInteger zero;
  const BigInteger inlineZero(BigInteger::InlineLimbs, true);
  const BigInteger heapZero(BigInteger::InlineLimbs + 1, false);
  const BigInteger one(std::vector<DataT>{1}, false);
  CheckEquality(zero, inlineZero, true);
  CheckEquality(zero, heapZero, true);
  CheckEquality(inlineZero, heapZero, true);
  CheckEquality(heapZero, BigInteger(heapZero), true);
  CheckEquality(one, heapZero, false);
}

REGISTER_TEST(Compare, EqualPositive)
{
  BigInteger a = BigIntegerBuilder::From("123456789");
  BigInteger b = BigIntegerBuilder::From("123456789");
  ASSERT_EQ(a.CompareTo(b), 0);
  ASSERT_TRUE(a == b);
  ASSERT_TRUE(a <= b);
  ASSERT_TRUE(a >= b);
  ASSERT_FALSE(a < b);
  ASSERT_FALSE(a > b);
  ASSERT_FALSE(a != b);
}

REGISTER_TEST(Compare, EqualNegative)
{
  BigInteger a = BigIntegerBuilder::From("-9999999999");
  BigInteger b = BigIntegerBuilder::From("-9999999999");
  ASSERT_TRUE(a == b);
  ASSERT_EQ(a.CompareTo(b), 0);
}

REGISTER_TEST(Compare, PositiveGreaterThanNegative)
{
  BigInteger a = BigIntegerBuilder::From("1");
  BigInteger b = BigIntegerBuilder::From("-1000000000000");
  ASSERT_GT(a, b);
  ASSERT_LT(b, a);
  ASSERT_EQ(a.CompareTo(b), 1);
  ASSERT_EQ(b.CompareTo(a), -1);
}

REGISTER_TEST(Compare, NegativeOrderReversed)
{
  // |-5| > |-3|  ⇒  -5 < -3
  BigInteger a = BigIntegerBuilder::From("-5");
  BigInteger b = BigIntegerBuilder::From("-3");
  ASSERT_LT(a, b);
  ASSERT_GT(b, a);
}

REGISTER_TEST(Compare, ZeroVsNegativeZero)
{
  BigInteger pos;
  BigInteger neg(std::vector<DataT>{0, 0}, true);
  ASSERT_EQ(pos.CompareTo(neg), 0);
  ASSERT_TRUE(pos == neg);
}

REGISTER_TEST(Compare, MultiLimbBoundary)
{
  // Difference in a high limb of a 4-limb value.
  std::vector<DataT> av{1, 2, 3, 4};
  std::vector<DataT> bv{1, 2, 3, 5};
  BigInteger a(av, false), b(bv, false);
  ASSERT_LT(a, b);
  // Same high limbs but a longer tail (b is bigger).
  std::vector<DataT> cv{1, 2, 3};
  std::vector<DataT> dv{1, 2, 3, 1};
  BigInteger c(cv, false), d(dv, false);
  ASSERT_LT(c, d);
}

REGISTER_TEST(Compare, RawVectorCompare)
{
  std::vector<DataT> a{1, 2, 3};
  std::vector<DataT> b{1, 2, 3};
  ASSERT_EQ(Compare(a, b), 0);

  std::vector<DataT> c{1, 2, 4};
  ASSERT_LT(Compare(a, c), 0);
  ASSERT_GT(Compare(c, a), 0);

  std::vector<DataT> d{1, 2, 3, 0};   // canonical equality with trailing-zero pad
  ASSERT_EQ(Compare(a, d), 0);
}

REGISTER_TEST(Compare, NumericOrderAcrossBlocks)
{
  const auto check = [](BigInteger const &a, BigInteger const &b, Int order) {
    ASSERT_EQ(a.CompareTo(b), order);
    ASSERT_EQ(b.CompareTo(a), -order);
    ASSERT_EQ(a > b, order > 0);
    ASSERT_EQ(a >= b, order >= 0);
    ASSERT_EQ(a < b, order < 0);
    ASSERT_EQ(a <= b, order <= 0);
  };
  for (SizeT size : {17u, 18u, 19u, 31u, 32u, 33u, 65u})
  {
    const std::vector<DataT> a(size, 0x100);
    check(BigInteger(a, false), BigInteger(a, false), 0);
    for (SizeT index = 1; index < size; ++index)
    {
      auto b = a;
      b.front() = LimbMask;
      b[index] = 0xff;
      // The higher difference decides even though bytewise comparison on
      // little-endian hardware and the low limb both favor the other side.
      ASSERT_EQ(Compare(a, b), 1);
      check(BigInteger(a, false), BigInteger(b, false), 1);
      check(BigInteger(a, true), BigInteger(b, true), -1);
    }
  }

  check(BigInteger(), BigInteger(65), 0);
  check(BigIntegerBuilder::From("-1"), BigInteger(65), -1);
}

REGISTER_TEST(Compare, RawPaddingAndScalar)
{
  const std::vector<DataT> empty;
  const std::vector<DataT> zeros(70, 0);
  std::vector<DataT> padded(75, 0);
  padded.front() = 5;
  const std::vector<DataT> five{5};
  ASSERT_EQ(Compare(empty, zeros), 0);
  ASSERT_EQ(Compare(zeros, empty), 0);
  ASSERT_EQ(Compare(empty, DataT{0}), 0);
  ASSERT_EQ(Compare(padded, five), 0);
  ASSERT_EQ(Compare(padded, DataT{5}), 0);
  ASSERT_EQ(Compare(padded, DataT{6}), -1);
  ASSERT_EQ(Compare(padded, DataT{4}), 1);
  for (SizeT i : {1u, 16u, 17u, 32u, 73u, 74u})
  {
    padded[i] = 1;
    ASSERT_EQ(Compare(padded, five), 1);
    ASSERT_EQ(Compare(five, padded), -1);
    ASSERT_EQ(Compare(padded, LimbMask), 1);
    padded[i] = 0;
  }
}
