// Comparison operators across sign permutations and magnitudes.

#include "unit_test_framework.h"

#include "biginteger/BigInteger.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Comparator.h"
#include "biginteger/ops/Comparison.h"

using namespace BigMath;

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
