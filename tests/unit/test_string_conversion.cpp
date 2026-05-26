// Parsing + formatting. Round-trips across linear and divide-conquer thresholds
// (DecimalDcThreshold = 8192 for parse, ToStringDcThreshold = 2048 for format).

#include "unit_test_framework.h"

#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "biginteger/BigInteger.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Parser.h"
#include "biginteger/ops/Comparison.h"
#include "biginteger/ops/IO.h"

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

// ─── parse basics ────────────────────────────────────────────────────────────

REGISTER_TEST(Parse, EmptyAndNullSafe)
{
  BigInteger e = Parse("");
  ASSERT_TRUE(e.Zero());
  BigInteger n = Parse(nullptr);
  ASSERT_TRUE(n.Zero());
}

REGISTER_TEST(Parse, JustZero)
{
  BigInteger x = Parse("0");
  ASSERT_TRUE(x.Zero());
  ASSERT_FALSE(x.IsNegative());
}

REGISTER_TEST(Parse, LeadingPlus)
{
  BigInteger x = Parse("+12345");
  ASSERT_FALSE(x.IsNegative());
  ASSERT_EQ(ToString(x), "12345");
}

REGISTER_TEST(Parse, LeadingMinus)
{
  BigInteger x = Parse("-12345");
  ASSERT_TRUE(x.IsNegative());
  ASSERT_EQ(ToString(x), "-12345");
}

REGISTER_TEST(Parse, LeadingZeros)
{
  BigInteger x = Parse("000000000042");
  ASSERT_EQ(ToString(x), "42");
}

REGISTER_TEST(Parse, AllZeros)
{
  BigInteger x = Parse("0000000");
  ASSERT_TRUE(x.Zero());
  ASSERT_FALSE(x.IsNegative());
}

REGISTER_TEST(Parse, NegativeZeroNormalized)
{
  BigInteger x = Parse("-0");
  ASSERT_TRUE(x.Zero());
  ASSERT_FALSE(x.IsNegative());
}

REGISTER_TEST(Parse, StopsAtNonDigit)
{
  Int processed = 0;
  BigInteger x = Parse("123abc", &processed);
  ASSERT_EQ(ToString(x), "123");
  ASSERT_EQ(processed, 3);
}

REGISTER_TEST(Parse, StopsAtNonDigitWithSign)
{
  Int processed = 0;
  BigInteger x = Parse("-7xyz", &processed);
  ASSERT_EQ(ToString(x), "-7");
  ASSERT_EQ(processed, 2);
}

REGISTER_TEST(Parse, EighteenDigitChunkBoundary)
{
  // ParseUnsignedLinear groups by 10^18. Exercise lengths around 18.
  for (int digits : {17, 18, 19, 36, 37}) {
    std::string s(digits, '9');
    BigInteger x = Parse(s.c_str());
    ASSERT_EQ(ToString(x), s);
  }
}

// ─── ToString basics ─────────────────────────────────────────────────────────

REGISTER_TEST(ToString, ZeroIsZero)
{
  BigInteger z;
  ASSERT_EQ(ToString(z), "0");
}

REGISTER_TEST(ToString, PadDoesNotApplyOnPublicApi)
{
  BigInteger x = BigIntegerBuilder::From("42");
  ASSERT_EQ(ToString(x), "42");
}

REGISTER_TEST(ToString, NegativeFormatted)
{
  BigInteger x = BigIntegerBuilder::From("-99");
  ASSERT_EQ(ToString(x), "-99");
}

// ─── round-trip across thresholds ────────────────────────────────────────────

static void RoundTripDigits(int digits, uint32_t seed)
{
  std::mt19937 gen(seed);
  std::string s = RandomDigits(digits, gen);
  BigInteger a = Parse(s.c_str());
  ASSERT_EQ(ToString(a), s);
}

REGISTER_TEST(RoundTrip, Small_1)         { RoundTripDigits(1,     0x1); }
REGISTER_TEST(RoundTrip, Small_18)        { RoundTripDigits(18,    0x2); }
REGISTER_TEST(RoundTrip, Small_100)       { RoundTripDigits(100,   0x3); }

// Below format DC threshold (2048).
REGISTER_TEST(RoundTrip, Mid_1000)        { RoundTripDigits(1000,  0x4); }
// Crosses format DC threshold (2048) but below parse DC threshold (8192).
REGISTER_TEST(RoundTrip, AboveFormatDc_3000)  { RoundTripDigits(3000,  0x5); }
// Just below parse DC threshold.
REGISTER_TEST(RoundTrip, NearParseDc_8000)    { RoundTripDigits(8000,  0x6); }
// Above parse DC threshold (8192).
REGISTER_TEST(RoundTrip, AboveParseDc_12000)  { RoundTripDigits(12000, 0x7); }
// Substantially above both.
REGISTER_TEST(RoundTrip, Large_25000)         { RoundTripDigits(25000, 0x8); }

// ─── ostream / istream ───────────────────────────────────────────────────────

REGISTER_TEST(IO, OstreamMatchesToString)
{
  BigInteger x = BigIntegerBuilder::From("-12345678901234567890");
  std::ostringstream os;
  os << x;
  ASSERT_EQ(os.str(), ToString(x));
}

REGISTER_TEST(IO, IstreamRoundTrip)
{
  BigInteger original = BigIntegerBuilder::From("99887766554433221100");
  std::ostringstream os;
  os << original;
  std::istringstream is(os.str());
  BigInteger parsed;
  is >> parsed;
  ASSERT_EQ(parsed, original);
}

// ─── Pow10 sanity ────────────────────────────────────────────────────────────

REGISTER_TEST(Pow10, SmallPowers)
{
  // Pow10 returns a base-2^32 limb vector; round-trip via ToString.
  ASSERT_EQ(ToString(Pow10(0)),  "1");
  ASSERT_EQ(ToString(Pow10(1)),  "10");
  ASSERT_EQ(ToString(Pow10(18)), "1000000000000000000");
  ASSERT_EQ(ToString(Pow10(19)), "10000000000000000000");
}

REGISTER_TEST(Pow10, LargePowerMatchesString)
{
  // 10^100 has 101 chars: '1' followed by 100 zeros.
  std::string expected = "1" + std::string(100, '0');
  ASSERT_EQ(ToString(Pow10(100)), expected);

  std::string expected1000 = "1" + std::string(1000, '0');
  ASSERT_EQ(ToString(Pow10(1000)), expected1000);
}
