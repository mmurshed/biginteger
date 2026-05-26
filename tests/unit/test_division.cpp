// Division and modulo: identity q*b+r==a, |r|<|b|, sign rules, scalar.
// Sizes span the FastDivision / BurnikelZiegler / Newton dispatch bands.

#include "unit_test_framework.h"

#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "biginteger/BigInteger.h"
#include "biginteger/algorithms/Division.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Comparator.h"
#include "biginteger/common/Parser.h"
#include "biginteger/ops/Addition.h"
#include "biginteger/ops/Comparison.h"
#include "biginteger/ops/Division.h"
#include "biginteger/ops/Multiplication.h"
#include "biginteger/ops/ScalarDivision.h"

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

// Magnitude-level identity check using the canonical limb-level dispatcher.
static void CheckRawIdentity(SizeT aLimbs, SizeT bLimbs, uint32_t seed)
{
  std::mt19937_64 gen(seed);
  std::uniform_int_distribution<uint64_t> digit(0, 0xFFFFFFFFULL);
  std::vector<DataT> a(aLimbs), b(bLimbs);
  for (auto &x : a) x = digit(gen);
  for (auto &x : b) x = digit(gen);
  if (a.back() == 0) a.back() = 1;
  if (b.back() == 0) b.back() = 1;

  auto qr = DivideAndRemainder(a, b, BigInteger::Base());
  std::vector<DataT> q = qr.first;
  std::vector<DataT> r = qr.second;

  // q*b + r == a
  std::vector<DataT> qb = Multiply(q, b, BigInteger::Base());
  std::vector<DataT> back = Add(qb, r, BigInteger::Base());
  ASSERT_EQ(Compare(back, a), 0);

  // r < b
  ASSERT_LT(Compare(r, b), 0);
}

// ─── single-limb divisor (scalar path) ───────────────────────────────────────

REGISTER_TEST(DivScalar, ByOne)
{
  BigInteger a = BigIntegerBuilder::From("12345678901234567890");
  ASSERT_EQ(a / (DataT)1, a);
  ASSERT_TRUE((a % (DataT)1).Zero());
}

REGISTER_TEST(DivScalar, ZeroOverAny)
{
  BigInteger z;
  ASSERT_TRUE((z / (DataT)7).Zero());
  ASSERT_TRUE((z % (DataT)7).Zero());
}

REGISTER_TEST(DivScalar, KnownValue)
{
  BigInteger a = BigIntegerBuilder::From("1000000007");
  BigInteger q = a / (DataT)13;
  BigInteger r = a % (DataT)13;
  ASSERT_EQ(ToString(q), "76923077");
  ASSERT_EQ(ToString(r), "6");
}

REGISTER_TEST(DivScalar, ByZeroThrows)
{
  BigInteger a = BigIntegerBuilder::From("42");
  bool threw = false;
  try { (void)(a / (DataT)0); }
  catch (const std::invalid_argument &) { threw = true; }
  ASSERT_TRUE(threw);
}

REGISTER_TEST(DivScalar, NegativeNumerator)
{
  BigInteger a = BigIntegerBuilder::From("-100");
  auto [q, r] = DivideAndRemainder(a, (DataT)7);
  // Truncate toward zero: q = -14, r = -2.  Library makes both negative.
  ASSERT_TRUE(q.IsNegative());
  ASSERT_TRUE(r.IsNegative());
  ASSERT_EQ(ToString(q), "-14");
  ASSERT_EQ(ToString(r), "-2");
}

// ─── BigInt / BigInt: small known values ─────────────────────────────────────

REGISTER_TEST(Div, ByOne)
{
  BigInteger a = BigIntegerBuilder::From("99999999999999999999");
  BigInteger one = BigIntegerBuilder::From("1");
  ASSERT_EQ(a / one, a);
  ASSERT_TRUE((a % one).Zero());
}

REGISTER_TEST(Div, BySelfIsOne)
{
  BigInteger a = BigIntegerBuilder::From("123456789012345678901234567890");
  ASSERT_EQ(ToString(a / a), "1");
  ASSERT_TRUE((a % a).Zero());
}

REGISTER_TEST(Div, SmallerOverLarger)
{
  BigInteger a = BigIntegerBuilder::From("3");
  BigInteger b = BigIntegerBuilder::From("1000000000");
  ASSERT_TRUE((a / b).Zero());
  ASSERT_EQ(a % b, a);
}

REGISTER_TEST(Div, KnownValue)
{
  // (10^30 + 7) / 1234567
  BigInteger a = BigIntegerBuilder::From("1000000000000000000000000000007");
  BigInteger b = BigIntegerBuilder::From("1234567");
  auto [q, r] = DivideAndRemainder(a, b);
  // Verify identity.
  BigInteger back = q * b + r;
  ASSERT_EQ(back, a);
  // Remainder bound: 0 ≤ r < b
  ASSERT_LT(r, b);
  ASSERT_FALSE(r.IsNegative());
}

REGISTER_TEST(Div, SignRulesQuotient)
{
  // q sign = XOR(a, b)
  BigInteger a = BigIntegerBuilder::From("100");
  BigInteger b = BigIntegerBuilder::From("7");
  ASSERT_FALSE((a / b).IsNegative());

  BigInteger na = BigIntegerBuilder::From("-100");
  ASSERT_TRUE((na / b).IsNegative());

  BigInteger nb = BigIntegerBuilder::From("-7");
  ASSERT_TRUE((a / nb).IsNegative());
  ASSERT_FALSE((na / nb).IsNegative());
}

REGISTER_TEST(Div, ByZeroThrows)
{
  BigInteger a = BigIntegerBuilder::From("42");
  BigInteger z;
  bool threw = false;
  try { (void)(a / z); }
  catch (...) { threw = true; }
  ASSERT_TRUE(threw);
}

// ─── magnitude-level dispatch bands ──────────────────────────────────────────

REGISTER_TEST(DivDispatch, FastBand_Small)        { CheckRawIdentity(32,   8,   0x10); }
REGISTER_TEST(DivDispatch, FastBand_Balanced_100) { CheckRawIdentity(100,  100, 0x11); }
REGISTER_TEST(DivDispatch, BZBand_Balanced_1024)  { CheckRawIdentity(1024, 1024, 0x12); }
REGISTER_TEST(DivDispatch, NewtonBand_Skewed)
{
  // a much larger than b, with b above NEWTON_MEDIUM_B → Newton band.
  CheckRawIdentity(4096, 1100, 0x13);
}
REGISTER_TEST(DivDispatch, NewtonBand_VeryLarge)  { CheckRawIdentity(8192, 2048, 0x14); }
REGISTER_TEST(DivDispatch, BZBand_AlmostBalanced) { CheckRawIdentity(1200, 1024, 0x15); }

// ─── multiplicative round-trip on the public API ─────────────────────────────

REGISTER_TEST(Div, MulThenDivRecovers)
{
  std::mt19937 gen(0x999);
  for (int trial = 0; trial < 20; ++trial)
  {
    int da = 1 + (gen() % 300);
    int db = 1 + (gen() % 200);
    BigInteger a = BigIntegerBuilder::From(RandomDigits(da, gen));
    BigInteger b = BigIntegerBuilder::From(RandomDigits(db, gen));
    BigInteger prod = a * b;
    ASSERT_EQ(prod / b, a);
    ASSERT_TRUE((prod % b).Zero());
  }
}

REGISTER_TEST(Div, IdentityHolds)
{
  // q*b + r == a  and  0 ≤ r < b  for positive operands.
  std::mt19937 gen(0xAAA);
  for (int trial = 0; trial < 15; ++trial)
  {
    int da = 50 + (gen() % 300);
    int db = 10 + (gen() % 100);
    BigInteger a = BigIntegerBuilder::From(RandomDigits(da, gen));
    BigInteger b = BigIntegerBuilder::From(RandomDigits(db, gen));
    auto [q, r] = DivideAndRemainder(a, b);
    ASSERT_EQ(q * b + r, a);
    ASSERT_LT(r, b);
    ASSERT_FALSE(r.IsNegative());
  }
}
