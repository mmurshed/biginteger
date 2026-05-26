// BigInteger multiplication: scalar and BigInt × BigInt, hitting each
// dispatch band (Classic / Karatsuba / NTT). Sizes are chosen relative to
// the runtime thresholds exposed in algorithms/Multiplication.h.

#include "unit_test_framework.h"

#include <random>
#include <string>
#include <vector>

#include "biginteger/BigInteger.h"
#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Parser.h"
#include "biginteger/ops/Addition.h"
#include "biginteger/ops/Comparison.h"
#include "biginteger/ops/Multiplication.h"
#include "biginteger/ops/ScalarMultiplication.h"
#include "biginteger/ops/Subtraction.h"

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

// ─── scalar multiplication ───────────────────────────────────────────────────

REGISTER_TEST(MulScalar, ZeroTimesAny)
{
  BigInteger a = BigIntegerBuilder::From("123456789012345678901234567890");
  BigInteger r = a * (DataT)0;
  ASSERT_TRUE(r.Zero());
}

REGISTER_TEST(MulScalar, OneTimesAny)
{
  BigInteger a = BigIntegerBuilder::From("123456789012345678901234567890");
  BigInteger r = a * (DataT)1;
  ASSERT_EQ(r, a);
}

REGISTER_TEST(MulScalar, NineByNine)
{
  BigInteger a = BigIntegerBuilder::From("999999999999999999");
  BigInteger r = a * (DataT)9;
  ASSERT_EQ(ToString(r), "8999999999999999991");
}

REGISTER_TEST(MulScalar, KeepsSign)
{
  BigInteger a = BigIntegerBuilder::From("-7");
  BigInteger r = a * (DataT)6;
  ASSERT_TRUE(r.IsNegative());
  ASSERT_EQ(ToString(r), "-42");
}

// ─── BigInt × BigInt: identities ─────────────────────────────────────────────

REGISTER_TEST(Mul, ZeroTimesBig)
{
  BigInteger z;
  BigInteger a = BigIntegerBuilder::From("123456789012345678");
  ASSERT_TRUE((z * a).Zero());
  ASSERT_TRUE((a * z).Zero());
}

REGISTER_TEST(Mul, OneTimesBig)
{
  BigInteger one = BigIntegerBuilder::From("1");
  BigInteger a = BigIntegerBuilder::From("123456789012345678");
  ASSERT_EQ(one * a, a);
  ASSERT_EQ(a * one, a);
}

REGISTER_TEST(Mul, NegativeOneFlipsSign)
{
  BigInteger m1 = BigIntegerBuilder::From("-1");
  BigInteger a = BigIntegerBuilder::From("424242424242");
  BigInteger r = m1 * a;
  ASSERT_TRUE(r.IsNegative());
  ASSERT_EQ(ToString(r), "-424242424242");
}

REGISTER_TEST(Mul, SignRules)
{
  BigInteger a = BigIntegerBuilder::From("7");
  BigInteger b = BigIntegerBuilder::From("11");
  ASSERT_EQ(ToString(a * b), "77");
  ASSERT_EQ(ToString((-a) * b), "-77");
  -a; // restore a
  ASSERT_EQ(ToString(a * (-b)), "-77");
  -b; // restore b
  ASSERT_EQ(ToString((-a) * (-b)), "77");
}

REGISTER_TEST(Mul, ProductSquare)
{
  // 1234567890123456789 × 1234567890123456789
  BigInteger a = BigIntegerBuilder::From("1234567890123456789");
  BigInteger r = a * a;
  ASSERT_EQ(ToString(r), "1524157875323883675019051998750190521");
}

REGISTER_TEST(Mul, Commutative)
{
  std::mt19937 gen(0xD44);
  for (int trial = 0; trial < 15; ++trial)
  {
    int da = 1 + (gen() % 400);
    int db = 1 + (gen() % 400);
    BigInteger a = BigIntegerBuilder::From(RandomDigits(da, gen));
    BigInteger b = BigIntegerBuilder::From(RandomDigits(db, gen));
    ASSERT_EQ(a * b, b * a);
  }
}

REGISTER_TEST(Mul, DistributesOverAdd)
{
  std::mt19937 gen(0xE55);
  for (int trial = 0; trial < 10; ++trial)
  {
    int da = 1 + (gen() % 250);
    int db = 1 + (gen() % 250);
    int dc = 1 + (gen() % 250);
    BigInteger a = BigIntegerBuilder::From(RandomDigits(da, gen));
    BigInteger b = BigIntegerBuilder::From(RandomDigits(db, gen));
    BigInteger c = BigIntegerBuilder::From(RandomDigits(dc, gen));
    ASSERT_EQ(a * (b + c), a * b + a * c);
  }
}

// ─── dispatch bands ──────────────────────────────────────────────────────────
// The bands are configured at link time; we exercise small/mid/large operand
// sizes so that whichever algorithm is selected, the limb-level dispatcher
// agrees with ClassicMultiplication on the same input.

static void CrossCheckAgainstClassic(int digits)
{
  std::mt19937 gen(0xF66 + digits);
  BigInteger a = BigIntegerBuilder::From(RandomDigits(digits, gen));
  BigInteger b = BigIntegerBuilder::From(RandomDigits(digits, gen));
  std::vector<DataT> dispatched = Multiply(a.GetInteger(), b.GetInteger(), BigInteger::Base());
  std::vector<DataT> classic = ClassicMultiplication::Multiply(
      a.GetInteger(), b.GetInteger(), BigInteger::Base());
  ASSERT_EQ(Compare(dispatched, classic), 0);
}

REGISTER_TEST(MulDispatch, Small32Digits)        { CrossCheckAgainstClassic(32);    }
REGISTER_TEST(MulDispatch, ClassicBand_100)       { CrossCheckAgainstClassic(100);   }
REGISTER_TEST(MulDispatch, KaratsubaBand_2000)    { CrossCheckAgainstClassic(2000);  }
REGISTER_TEST(MulDispatch, KaratsubaBand_5000)    { CrossCheckAgainstClassic(5000);  }
REGISTER_TEST(MulDispatch, NTTBand_20000)         { CrossCheckAgainstClassic(20000); }

REGISTER_TEST(MulDispatch, AsymmetricSizes)
{
  std::mt19937 gen(0x123);
  BigInteger a = BigIntegerBuilder::From(RandomDigits(8000, gen));
  BigInteger b = BigIntegerBuilder::From(RandomDigits(50, gen));
  std::vector<DataT> dispatched = Multiply(a.GetInteger(), b.GetInteger(), BigInteger::Base());
  std::vector<DataT> classic = ClassicMultiplication::Multiply(
      a.GetInteger(), b.GetInteger(), BigInteger::Base());
  ASSERT_EQ(Compare(dispatched, classic), 0);
}

REGISTER_TEST(MulDispatch, ScalarFastPath)
{
  // Single-limb b operand should route through the scalar dispatcher path.
  std::vector<DataT> a{0xDEADBEEF, 0x12345678, 0xCAFEBABE, 0x9};
  std::vector<DataT> b{3};
  std::vector<DataT> dispatched = Multiply(a, b, BigInteger::Base());
  std::vector<DataT> classic = ClassicMultiplication::Multiply(a, b, BigInteger::Base());
  ASSERT_EQ(Compare(dispatched, classic), 0);
}
