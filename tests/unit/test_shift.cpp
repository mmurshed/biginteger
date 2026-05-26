// BigInteger << / >> operate at LIMB granularity (base-2^32), not bit.
// So  a << k  ≡  a · (2^32)^k  and  a >> k ≡ a / (2^32)^k.

#include "unit_test_framework.h"

#include "biginteger/BigInteger.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Parser.h"
#include "biginteger/ops/Comparison.h"
#include "biginteger/ops/Multiplication.h"
#include "biginteger/ops/Shift.h"

using namespace BigMath;

REGISTER_TEST(Shift, LeftByZeroIsIdentity)
{
  BigInteger a = BigIntegerBuilder::From("12345678901234567890");
  ASSERT_EQ(a << 0, a);
}

REGISTER_TEST(Shift, RightByZeroIsIdentity)
{
  BigInteger a = BigIntegerBuilder::From("12345678901234567890");
  ASSERT_EQ(a >> 0, a);
}

REGISTER_TEST(Shift, ZeroShiftedIsZero)
{
  BigInteger z;
  ASSERT_TRUE((z << 5).Zero());
  ASSERT_TRUE((z >> 5).Zero());
}

#if BIGMATH_LIMB_64
REGISTER_TEST(Shift, LeftOneLimbEquals2Pow64)
{
  BigInteger one = BigIntegerBuilder::From("1");
  BigInteger r = one << 1;          // 1 * 2^64
  ASSERT_EQ(ToString(r), "18446744073709551616");
}

REGISTER_TEST(Shift, LeftTwoLimbsEquals2Pow128)
{
  BigInteger one = BigIntegerBuilder::From("1");
  BigInteger r = one << 2;          // 1 * 2^128
  ASSERT_EQ(ToString(r), "340282366920938463463374607431768211456");
}

REGISTER_TEST(Shift, LeftMatchesMultiplyByPowOfTwo64)
{
  // a << 3  ==  a * (2^64)^3 = a * 2^192
  BigInteger a = BigIntegerBuilder::From("123456789012345");
  BigInteger pow = BigIntegerBuilder::From(
      "6277101735386680763835789423207666416102355444464034512896");  // 2^192
  ASSERT_EQ(a << 3, a * pow);
}
#else
REGISTER_TEST(Shift, LeftOneLimbEquals2Pow32)
{
  BigInteger one = BigIntegerBuilder::From("1");
  BigInteger r = one << 1;          // 1 * 2^32
  ASSERT_EQ(ToString(r), "4294967296");
}

REGISTER_TEST(Shift, LeftTwoLimbsEquals2Pow64)
{
  BigInteger one = BigIntegerBuilder::From("1");
  BigInteger r = one << 2;          // 1 * 2^64
  ASSERT_EQ(ToString(r), "18446744073709551616");
}

REGISTER_TEST(Shift, LeftMatchesMultiplyByPowOfTwo32)
{
  // a << 3  ==  a * (2^32)^3
  BigInteger a = BigIntegerBuilder::From("123456789012345");
  BigInteger pow = BigIntegerBuilder::From(
      "6277101735386680763835789423207666416102355444464034512896");  // 2^192
  ASSERT_EQ(a << 6, a * pow);
}
#endif

REGISTER_TEST(Shift, RightInvertsLeftWhenMultipleOfLimb)
{
  BigInteger a = BigIntegerBuilder::From("123456789012345678901234567890");
  ASSERT_EQ((a << 4) >> 4, a);
}

REGISTER_TEST(Shift, RightOverflowProducesZero)
{
  // Right-shift by more limbs than the value holds → 0.
  BigInteger a = BigIntegerBuilder::From("1");   // 1 limb
  ASSERT_TRUE((a >> 5).Zero());
}

REGISTER_TEST(Shift, RightTruncatesLowLimbs)
{
  // (a · 2^64 + r) >> 2 == a    when r < 2^64
  BigInteger a = BigIntegerBuilder::From("999999999999");
  BigInteger r = (a << 2) >> 2;
  ASSERT_EQ(r, a);
}

REGISTER_TEST(Shift, SignPreserved)
{
  BigInteger a = BigIntegerBuilder::From("-123456789");
  BigInteger s = a << 2;
  ASSERT_TRUE(s.IsNegative());
  ASSERT_EQ(s >> 2, a);
}
