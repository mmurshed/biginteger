// Sign-aware Add / Subtract on BigInteger.
// Magnitudes range from 0 up to ~600 limbs so per-limb carries and borrows
// across multiple limb boundaries are exercised.

#include "unit_test_framework.h"

#include <random>
#include <string>

#include "biginteger/BigInteger.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Parser.h"
#include "biginteger/ops/Addition.h"
#include "biginteger/ops/Subtraction.h"
#include "biginteger/ops/Comparison.h"

using namespace BigMath;

static std::string RandomDigits(int digits, std::mt19937 &gen)
{
  std::uniform_int_distribution<int> first(1, 9);
  std::uniform_int_distribution<int> rest(0, 9);
  std::string v;
  v.reserve(digits);
  v.push_back((char)('0' + first(gen)));
  for (int i = 1; i < digits; ++i)
    v.push_back((char)('0' + rest(gen)));
  return v;
}

// ─── Add ─────────────────────────────────────────────────────────────────────

REGISTER_TEST(Add, ZeroIdentity)
{
  BigInteger a = BigIntegerBuilder::From("12345678901234567890");
  BigInteger z;
  ASSERT_EQ(Add(a, z), a);
  ASSERT_EQ(Add(z, a), a);
  ASSERT_EQ(a + z, a);
  ASSERT_EQ(z + a, a);
}

REGISTER_TEST(Add, BothZero)
{
  BigInteger z;
  BigInteger r = z + z;
  ASSERT_TRUE(r.Zero());
  ASSERT_FALSE(r.IsNegative());
}

REGISTER_TEST(Add, SmallPositive)
{
  BigInteger a = BigIntegerBuilder::From("999");
  BigInteger b = BigIntegerBuilder::From("1");
  ASSERT_EQ(ToString(a + b), "1000");
}

REGISTER_TEST(Add, MultiLimbCarry)
{
  // 2^32-1 + 1 → 2^32 (carries across limb).
  BigInteger a = BigIntegerBuilder::From("4294967295");
  BigInteger b = BigIntegerBuilder::From("1");
  ASSERT_EQ(ToString(a + b), "4294967296");
}

REGISTER_TEST(Add, CarryChain)
{
  // (2^256 − 1) + 1 = 2^256, a full carry chain across many limbs.
  BigInteger a = BigIntegerBuilder::From(
      "115792089237316195423570985008687907853269984665640564039457584007913129639935");
  BigInteger b = BigIntegerBuilder::From("1");
  BigInteger r = a + b;
  ASSERT_EQ(ToString(r),
            "115792089237316195423570985008687907853269984665640564039457584007913129639936");
}

REGISTER_TEST(Add, Commutative)
{
  std::mt19937 gen(0xA11);
  for (int trial = 0; trial < 20; ++trial)
  {
    int da = 1 + (gen() % 500);
    int db = 1 + (gen() % 500);
    BigInteger a = BigIntegerBuilder::From(RandomDigits(da, gen));
    BigInteger b = BigIntegerBuilder::From(RandomDigits(db, gen));
    ASSERT_EQ(a + b, b + a);
  }
}

REGISTER_TEST(Add, Associative)
{
  std::mt19937 gen(0xB22);
  for (int trial = 0; trial < 10; ++trial)
  {
    int da = 1 + (gen() % 300);
    int db = 1 + (gen() % 300);
    int dc = 1 + (gen() % 300);
    BigInteger a = BigIntegerBuilder::From(RandomDigits(da, gen));
    BigInteger b = BigIntegerBuilder::From(RandomDigits(db, gen));
    BigInteger c = BigIntegerBuilder::From(RandomDigits(dc, gen));
    ASSERT_EQ((a + b) + c, a + (b + c));
  }
}

REGISTER_TEST(Add, NegativePlusPositiveSmaller)
{
  // -10 + 3 = -7
  BigInteger a = BigIntegerBuilder::From("-10");
  BigInteger b = BigIntegerBuilder::From("3");
  BigInteger r = a + b;
  ASSERT_EQ(ToString(r), "-7");
}

REGISTER_TEST(Add, NegativePlusPositiveLarger)
{
  // -3 + 10 = 7
  BigInteger a = BigIntegerBuilder::From("-3");
  BigInteger b = BigIntegerBuilder::From("10");
  BigInteger r = a + b;
  ASSERT_EQ(ToString(r), "7");
}

REGISTER_TEST(Add, NegativePlusNegative)
{
  BigInteger a = BigIntegerBuilder::From("-100");
  BigInteger b = BigIntegerBuilder::From("-250");
  ASSERT_EQ(ToString(a + b), "-350");
}

REGISTER_TEST(Add, CancellationToZero)
{
  BigInteger a = BigIntegerBuilder::From("12345678901234567890");
  BigInteger b = BigIntegerBuilder::From("-12345678901234567890");
  BigInteger r = a + b;
  ASSERT_TRUE(r.Zero());
  ASSERT_FALSE(r.IsNegative());
}

// ─── Subtract ────────────────────────────────────────────────────────────────

REGISTER_TEST(Sub, ZeroMinusZero)
{
  BigInteger z;
  BigInteger r = z - z;
  ASSERT_TRUE(r.Zero());
  ASSERT_FALSE(r.IsNegative());
}

REGISTER_TEST(Sub, AnyMinusZero)
{
  BigInteger a = BigIntegerBuilder::From("987654321987654321");
  BigInteger z;
  ASSERT_EQ(a - z, a);
}

REGISTER_TEST(Sub, ZeroMinusPositive)
{
  BigInteger z;
  BigInteger a = BigIntegerBuilder::From("42");
  ASSERT_EQ(ToString(z - a), "-42");
}

REGISTER_TEST(Sub, ZeroMinusNegative)
{
  BigInteger z;
  BigInteger a = BigIntegerBuilder::From("-42");
  ASSERT_EQ(ToString(z - a), "42");
}

REGISTER_TEST(Sub, EqualToZero)
{
  BigInteger a = BigIntegerBuilder::From("999999999999999");
  BigInteger r = a - a;
  ASSERT_TRUE(r.Zero());
  ASSERT_FALSE(r.IsNegative());
}

REGISTER_TEST(Sub, PositiveMinusSmallerPositive)
{
  BigInteger a = BigIntegerBuilder::From("1000000000000");
  BigInteger b = BigIntegerBuilder::From("1");
  ASSERT_EQ(ToString(a - b), "999999999999");
}

REGISTER_TEST(Sub, PositiveMinusLargerPositive)
{
  BigInteger a = BigIntegerBuilder::From("5");
  BigInteger b = BigIntegerBuilder::From("12");
  ASSERT_EQ(ToString(a - b), "-7");
}

REGISTER_TEST(Sub, BorrowChain)
{
  // 2^256 − 1: full borrow chain.
  BigInteger a = BigIntegerBuilder::From(
      "115792089237316195423570985008687907853269984665640564039457584007913129639936");
  BigInteger b = BigIntegerBuilder::From("1");
  BigInteger r = a - b;
  ASSERT_EQ(ToString(r),
            "115792089237316195423570985008687907853269984665640564039457584007913129639935");
}

REGISTER_TEST(Sub, OppositeSigns)
{
  // a - b where signs differ → magnitudes add.
  BigInteger a = BigIntegerBuilder::From("100");
  BigInteger b = BigIntegerBuilder::From("-50");
  ASSERT_EQ(ToString(a - b), "150");

  BigInteger c = BigIntegerBuilder::From("-100");
  BigInteger d = BigIntegerBuilder::From("50");
  ASSERT_EQ(ToString(c - d), "-150");
}

REGISTER_TEST(Sub, BothNegative)
{
  // -5 - (-3) = -2
  BigInteger a = BigIntegerBuilder::From("-5");
  BigInteger b = BigIntegerBuilder::From("-3");
  ASSERT_EQ(ToString(a - b), "-2");

  // -3 - (-5) = 2
  BigInteger c = BigIntegerBuilder::From("-3");
  BigInteger d = BigIntegerBuilder::From("-5");
  ASSERT_EQ(ToString(c - d), "2");
}

// ─── round-trip Add/Sub ──────────────────────────────────────────────────────

REGISTER_TEST(AddSub, AddThenSubtractRecoversOriginal)
{
  std::mt19937 gen(0xC33);
  for (int trial = 0; trial < 30; ++trial)
  {
    int da = 1 + (gen() % 400);
    int db = 1 + (gen() % 400);
    BigInteger a = BigIntegerBuilder::From(RandomDigits(da, gen));
    BigInteger b = BigIntegerBuilder::From(RandomDigits(db, gen));
    if (trial % 3 == 0) a.SetSign(true);
    if (trial % 4 == 0) b.SetSign(true);
    BigInteger sum = a + b;
    BigInteger back = sum - b;
    ASSERT_EQ(back, a);
  }
}
