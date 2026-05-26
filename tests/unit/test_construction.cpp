// Construction, accessors, sign normalization, builder.

#include "unit_test_framework.h"

#include <limits>
#include <string>

#include "biginteger/BigInteger.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Parser.h"
#include "biginteger/ops/Comparison.h"

using namespace BigMath;

REGISTER_TEST(Construction, DefaultIsZero)
{
  BigInteger x;
  ASSERT_TRUE(x.Zero());
  ASSERT_FALSE(x.IsNegative());
  ASSERT_EQ(x.size(), 1u);
  ASSERT_EQ(x[0], 0u);
}

REGISTER_TEST(Construction, SizedDefaultZero)
{
  BigInteger x(4, false);
  ASSERT_TRUE(x.Zero());
  ASSERT_EQ(x.size(), 4u);
  for (SizeT i = 0; i < x.size(); ++i)
    ASSERT_EQ(x[i], 0u);
}

REGISTER_TEST(Construction, NegativeZeroNormalized)
{
  // Constructing a zero BigInteger with negative=true must canonicalize to non-negative.
  BigInteger x(1, true);
  ASSERT_TRUE(x.Zero());
  ASSERT_FALSE(x.IsNegative());

  std::vector<DataT> zero{0, 0, 0};
  BigInteger y(zero, true);
  ASSERT_TRUE(y.Zero());
  ASSERT_FALSE(y.IsNegative());
}

REGISTER_TEST(Construction, FromVectorTrimmed)
{
  // Leading-zero limbs are trimmed; value preserved.
  std::vector<DataT> v{42, 0, 0, 0};
  BigInteger x(v, false);
  // Trimmed to one limb of value 42.
  ASSERT_EQ(x.size(), 1u);
  ASSERT_EQ(x[0], 42u);

  std::vector<DataT> w{1, 2, 0};
  BigInteger y(w, true);
  ASSERT_EQ(y.size(), 2u);
  ASSERT_EQ(y[0], 1u);
  ASSERT_EQ(y[1], 2u);
  ASSERT_TRUE(y.IsNegative());
}

REGISTER_TEST(Construction, SizedFill)
{
  BigInteger x(3, false, 7);
  ASSERT_FALSE(x.Zero());
  ASSERT_EQ(x.size(), 3u);
  ASSERT_EQ(x[0], 7u);
  ASSERT_EQ(x[1], 7u);
  ASSERT_EQ(x[2], 7u);
}

REGISTER_TEST(Construction, CopyAndAssign)
{
  BigInteger a = BigIntegerBuilder::From("-12345678901234567890");
  BigInteger b(a);
  ASSERT_EQ(a, b);
  ASSERT_EQ(a.IsNegative(), b.IsNegative());

  BigInteger c;
  c = a;
  ASSERT_EQ(a, c);

  // Self-assign safe.
  c = c;
  ASSERT_EQ(a, c);
}

REGISTER_TEST(Construction, UnaryNegationFlipsSign)
{
  BigInteger a = BigIntegerBuilder::From("123456789012345678");
  ASSERT_FALSE(a.IsNegative());
  -a;
  ASSERT_TRUE(a.IsNegative());
  -a;
  ASSERT_FALSE(a.IsNegative());
}

REGISTER_TEST(Construction, UnaryNegationZeroStaysPositive)
{
  BigInteger z;
  -z;
  ASSERT_FALSE(z.IsNegative());
  ASSERT_TRUE(z.Zero());
}

REGISTER_TEST(Construction, SetSignZeroIgnored)
{
  BigInteger z;
  z.SetSign(true);
  ASSERT_FALSE(z.IsNegative());
}

REGISTER_TEST(Construction, SetSignNonZero)
{
  BigInteger a = BigIntegerBuilder::From("42");
  a.SetSign(true);
  ASSERT_TRUE(a.IsNegative());
  a.SetSign(false);
  ASSERT_FALSE(a.IsNegative());
}

REGISTER_TEST(Builder, FromULongZero)
{
  BigInteger x = BigIntegerBuilder::From((ULong)0);
  ASSERT_TRUE(x.Zero());
  ASSERT_FALSE(x.IsNegative());
}

REGISTER_TEST(Builder, FromULongMax)
{
  ULong v = std::numeric_limits<ULong>::max();
  BigInteger x = BigIntegerBuilder::From(v);
  ASSERT_FALSE(x.IsNegative());
  ASSERT_EQ(ToString(x), std::to_string(v));
}

REGISTER_TEST(Builder, FromLongNegative)
{
  Long v = -123456789012345LL;
  BigInteger x = BigIntegerBuilder::From(v);
  ASSERT_TRUE(x.IsNegative());
  ASSERT_EQ(ToString(x), std::to_string(v));
}

REGISTER_TEST(Builder, FromString)
{
  BigInteger x = BigIntegerBuilder::From(std::string("100000000000000000000"));
  ASSERT_FALSE(x.IsNegative());
  ASSERT_EQ(ToString(x), "100000000000000000000");
}

REGISTER_TEST(Builder, FromCStrLeadingPlus)
{
  BigInteger x = BigIntegerBuilder::From("+999");
  ASSERT_FALSE(x.IsNegative());
  ASSERT_EQ(ToString(x), "999");
}

REGISTER_TEST(Builder, FromCStrLeadingZeros)
{
  BigInteger x = BigIntegerBuilder::From("00012345");
  ASSERT_EQ(ToString(x), "12345");
}

REGISTER_TEST(Builder, FromLongDoubleZero)
{
  BigInteger x = BigIntegerBuilder::From((long double)0.0L);
  ASSERT_TRUE(x.Zero());
  ASSERT_FALSE(x.IsNegative());
}

REGISTER_TEST(Builder, FromLongDoublePositive)
{
  BigInteger x = BigIntegerBuilder::From((long double)4294967296.0L);  // = 2^32
  ASSERT_EQ(ToString(x), "4294967296");
}

REGISTER_TEST(Builder, FromLongDoubleNegative)
{
  BigInteger x = BigIntegerBuilder::From((long double)-1000000.0L);
  ASSERT_TRUE(x.IsNegative());
  ASSERT_EQ(ToString(x), "-1000000");
}

#if !BIGMATH_LIMB_64
// Pins the Base2_32 limb encoding of 2^32. Under LIMB_64=1 the same value
// fits in one 64-bit limb (v = {2^32}); skip rather than assert wrong layout.
REGISTER_TEST(Builder, VectorFromULong)
{
  std::vector<DataT> v = BigIntegerBuilder::VectorFrom((ULong)Base2_32);  // 2^32
  // 2^32 = 1 * 2^32 + 0  →  limbs {0, 1}
  ASSERT_EQ(v.size(), 2u);
  ASSERT_EQ(v[0], 0u);
  ASSERT_EQ(v[1], 1u);
}
#endif
