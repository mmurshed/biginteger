// BigDecimal: construction, parse, arithmetic, division with rounding modes,
// scale manipulation, comparison, formatting.

#include "unit_test_framework.h"

#include <string>

#include "bigdecimal/BigDecimal.h"
#include "bigdecimal/RoundingMode.h"

using namespace BigMath;

namespace
{
  inline BigDecimal D(char const *s) { return BigDecimal::FromString(s); }
  inline std::string S(BigDecimal const &x) { return x.ToPlainString(); }
}

// ─── parse + format ──────────────────────────────────────────────────────────

REGISTER_TEST(BigDecimal, ParsePlain)
{
  ASSERT_EQ(S(D("0")),         "0");
  ASSERT_EQ(S(D("1")),         "1");
  ASSERT_EQ(S(D("-1")),        "-1");
  ASSERT_EQ(S(D("123.456")),   "123.456");
  ASSERT_EQ(S(D("-123.456")),  "-123.456");
  ASSERT_EQ(S(D("0.001")),     "0.001");
  ASSERT_EQ(S(D("-0.001")),    "-0.001");
  ASSERT_EQ(S(D("100")),       "100");
  ASSERT_EQ(S(D(".5")),        "0.5");
  ASSERT_EQ(S(D("5.")),        "5");
  ASSERT_EQ(S(D("+1.25")),     "1.25");
}

REGISTER_TEST(BigDecimal, ParseScientific)
{
  ASSERT_EQ(S(D("1e3")),       "1000");
  ASSERT_EQ(S(D("1.5e2")),     "150");
  ASSERT_EQ(S(D("1.5e-2")),    "0.015");
  ASSERT_EQ(S(D("123.456e0")), "123.456");
  ASSERT_EQ(S(D("-1.0E3")),    "-1000");
  ASSERT_EQ(S(D("2e-3")),      "0.002");
}

REGISTER_TEST(BigDecimal, ParseLeadingZeros)
{
  ASSERT_EQ(S(D("00123")),     "123");
  ASSERT_EQ(S(D("000.5")),     "0.5");
  ASSERT_EQ(S(D("00")),        "0");
}

REGISTER_TEST(BigDecimal, ParseErrors)
{
  bool threw = false;
  try { (void)D(""); } catch (std::invalid_argument const &) { threw = true; }
  ASSERT_TRUE(threw);

  threw = false;
  try { (void)D("abc"); } catch (std::invalid_argument const &) { threw = true; }
  ASSERT_TRUE(threw);

  threw = false;
  try { (void)D("1.2.3"); } catch (std::invalid_argument const &) { threw = true; }
  ASSERT_TRUE(threw);

  threw = false;
  try { (void)D("1e"); } catch (std::invalid_argument const &) { threw = true; }
  ASSERT_TRUE(threw);

  threw = false;
  try { (void)D("1.2x"); } catch (std::invalid_argument const &) { threw = true; }
  ASSERT_TRUE(threw);
}

REGISTER_TEST(BigDecimal, ScaleAccessor)
{
  ASSERT_EQ(D("1.23").Scale(), 2);
  ASSERT_EQ(D("100").Scale(), 0);
  ASSERT_EQ(D("1.5e2").Scale(), -1);     // 150 = unscaled 15 · 10^1
  ASSERT_EQ(D("0.001").Scale(), 3);
}

// ─── addition / subtraction ──────────────────────────────────────────────────

REGISTER_TEST(BigDecimal, Addition)
{
  ASSERT_EQ(S(D("1.5")  + D("2.5")),    "4.0");
  ASSERT_EQ(S(D("0.1")  + D("0.2")),    "0.3");
  ASSERT_EQ(S(D("100")  + D("0.01")),   "100.01");
  ASSERT_EQ(S(D("-5")   + D("3")),      "-2");
  ASSERT_EQ(S(D("0")    + D("0")),      "0");
  ASSERT_EQ(S(D("999999999999999999999999.99") + D("0.01")),
            "1000000000000000000000000.00");
}

REGISTER_TEST(BigDecimal, Subtraction)
{
  ASSERT_EQ(S(D("3.0")  - D("1.5")),    "1.5");
  ASSERT_EQ(S(D("0")    - D("0.001")),  "-0.001");
  ASSERT_EQ(S(D("100")  - D("99.99")),  "0.01");
  ASSERT_EQ(S(D("-3")   - D("-5")),     "2");      // both-negative path
  ASSERT_EQ(S(D("-3")   - D("5")),      "-8");
}

REGISTER_TEST(BigDecimal, AddSubInverse)
{
  BigDecimal a = D("3.14159265358979323846");
  BigDecimal b = D("2.71828182845904523536");
  ASSERT_TRUE((a + b - b).ValueEquals(a));
  ASSERT_TRUE((a - b + b).ValueEquals(a));
}

// ─── multiplication ──────────────────────────────────────────────────────────

REGISTER_TEST(BigDecimal, Multiplication)
{
  ASSERT_EQ(S(D("2")     * D("3")),       "6");
  ASSERT_EQ(S(D("1.5")   * D("2")),       "3.0");
  ASSERT_EQ(S(D("0.1")   * D("0.2")),     "0.02");
  ASSERT_EQ(S(D("-1.5")  * D("2")),       "-3.0");
  ASSERT_EQ(S(D("-1.5")  * D("-2")),      "3.0");
  // Java-style: scale of zero result is sum of input scales (0 + 2)
  ASSERT_EQ(S(D("0")     * D("123.45")),  "0.00");
}

REGISTER_TEST(BigDecimal, MultiplicationLarge)
{
  // 10! computed via repeated multiplication
  BigDecimal r = D("1");
  for (int i = 2; i <= 10; ++i)
    r = r * D(std::to_string(i).c_str());
  ASSERT_EQ(S(r), "3628800");
}

REGISTER_TEST(BigDecimal, ScaleAddsOnMultiply)
{
  ASSERT_EQ((D("1.23") * D("4.5")).Scale(), 3);   // 2 + 1
  ASSERT_EQ((D("1.00") * D("1.00")).Scale(), 4);
}

// ─── division ────────────────────────────────────────────────────────────────

REGISTER_TEST(BigDecimal, DivideExact)
{
  ASSERT_EQ(S(D("10").Divide(D("4"),  2, RoundingMode::HALF_UP)),    "2.50");
  ASSERT_EQ(S(D("10").Divide(D("4"),  0, RoundingMode::HALF_UP)),    "3");
  ASSERT_EQ(S(D("10").Divide(D("4"),  0, RoundingMode::HALF_DOWN)),  "2");
  ASSERT_EQ(S(D("10").Divide(D("4"),  0, RoundingMode::DOWN)),       "2");
  ASSERT_EQ(S(D("10").Divide(D("4"),  0, RoundingMode::UP)),         "3");
}

REGISTER_TEST(BigDecimal, DivideHalfUpHalfDown)
{
  // 1/3 to 10 places (HALF_UP, last digit "...3333333333")
  std::string s = S(D("1").Divide(D("3"), 10, RoundingMode::HALF_UP));
  ASSERT_EQ(s, "0.3333333333");

  // 2/3 to 10 places (rounds up at last digit)
  s = S(D("2").Divide(D("3"), 10, RoundingMode::HALF_UP));
  ASSERT_EQ(s, "0.6666666667");

  // tie cases: 0.5
  ASSERT_EQ(S(D("1").Divide(D("2"), 0, RoundingMode::HALF_UP)),    "1");
  ASSERT_EQ(S(D("1").Divide(D("2"), 0, RoundingMode::HALF_DOWN)),  "0");
  // 1.5 → tie
  ASSERT_EQ(S(D("3").Divide(D("2"), 0, RoundingMode::HALF_UP)),    "2");
  ASSERT_EQ(S(D("3").Divide(D("2"), 0, RoundingMode::HALF_DOWN)),  "1");
}

REGISTER_TEST(BigDecimal, DivideHalfEven)
{
  // 0.5 → 0 (round to even)
  ASSERT_EQ(S(D("1").Divide(D("2"), 0, RoundingMode::HALF_EVEN)),  "0");
  // 1.5 → 2 (round to even)
  ASSERT_EQ(S(D("3").Divide(D("2"), 0, RoundingMode::HALF_EVEN)),  "2");
  // 2.5 → 2 (round to even)
  ASSERT_EQ(S(D("5").Divide(D("2"), 0, RoundingMode::HALF_EVEN)),  "2");
  // 3.5 → 4 (round to even)
  ASSERT_EQ(S(D("7").Divide(D("2"), 0, RoundingMode::HALF_EVEN)),  "4");
  // non-tie: 0.6 → 1
  ASSERT_EQ(S(D("3").Divide(D("5"), 0, RoundingMode::HALF_EVEN)),  "1");
}

REGISTER_TEST(BigDecimal, DivideCeilingFloor)
{
  // 1/3 ≈ 0.333; CEILING → 0.34, FLOOR → 0.33
  ASSERT_EQ(S(D("1").Divide(D("3"),   2, RoundingMode::CEILING)),  "0.34");
  ASSERT_EQ(S(D("1").Divide(D("3"),   2, RoundingMode::FLOOR)),    "0.33");
  // -1/3; CEILING toward +inf → -0.33; FLOOR toward -inf → -0.34
  ASSERT_EQ(S(D("-1").Divide(D("3"),  2, RoundingMode::CEILING)),  "-0.33");
  ASSERT_EQ(S(D("-1").Divide(D("3"),  2, RoundingMode::FLOOR)),    "-0.34");
}

REGISTER_TEST(BigDecimal, DivideSignCombos)
{
  ASSERT_EQ(S(D("-10").Divide(D("4"),  2, RoundingMode::HALF_UP)), "-2.50");
  ASSERT_EQ(S(D("10").Divide(D("-4"),  2, RoundingMode::HALF_UP)), "-2.50");
  ASSERT_EQ(S(D("-10").Divide(D("-4"), 2, RoundingMode::HALF_UP)), "2.50");
}

REGISTER_TEST(BigDecimal, DivideByZeroThrows)
{
  bool threw = false;
  try { (void)D("1").Divide(D("0"), 5, RoundingMode::HALF_UP); }
  catch (std::invalid_argument const &) { threw = true; }
  ASSERT_TRUE(threw);
}

REGISTER_TEST(BigDecimal, DivideUnnecessaryThrowsWhenNeeded)
{
  bool threw = false;
  try { (void)D("1").Divide(D("3"), 5, RoundingMode::UNNECESSARY); }
  catch (std::runtime_error const &) { threw = true; }
  ASSERT_TRUE(threw);

  // exact: 10 / 4 = 2.5 at scale 1 — fits without rounding
  BigDecimal exact = D("10").Divide(D("4"), 1, RoundingMode::UNNECESSARY);
  ASSERT_EQ(S(exact), "2.5");
}

REGISTER_TEST(BigDecimal, DivideRoundTrip)
{
  // (a / b) * b should equal `a` up to rounding at sufficient precision
  BigDecimal a = D("1");
  BigDecimal b = D("7");
  BigDecimal q = a.Divide(b, 30, RoundingMode::HALF_EVEN);
  BigDecimal back = q * b;
  // back should be very close to 1; difference < 10^-29
  BigDecimal diff = (back - a);
  if (diff.IsNegative()) diff = diff.Negate();
  BigDecimal tol = D("1e-29");
  ASSERT_TRUE(diff <= tol);
}

// ─── scale manipulation ─────────────────────────────────────────────────────

REGISTER_TEST(BigDecimal, SetScaleUp)
{
  ASSERT_EQ(S(D("1.5").SetScale(4, RoundingMode::UNNECESSARY)), "1.5000");
  ASSERT_EQ(S(D("100").SetScale(2, RoundingMode::UNNECESSARY)), "100.00");
}

REGISTER_TEST(BigDecimal, SetScaleDown)
{
  ASSERT_EQ(S(D("1.2345").SetScale(2, RoundingMode::HALF_UP)),   "1.23");
  ASSERT_EQ(S(D("1.2355").SetScale(3, RoundingMode::HALF_UP)),   "1.236");
  ASSERT_EQ(S(D("1.2355").SetScale(3, RoundingMode::HALF_DOWN)), "1.235");
  ASSERT_EQ(S(D("1.2355").SetScale(3, RoundingMode::HALF_EVEN)), "1.236");
  ASSERT_EQ(S(D("1.2345").SetScale(3, RoundingMode::HALF_EVEN)), "1.234");
}

REGISTER_TEST(BigDecimal, StripTrailingZeros)
{
  ASSERT_EQ(S(D("1.2000").StripTrailingZeros()),  "1.2");
  ASSERT_EQ(S(D("100").StripTrailingZeros()),     "100");
  ASSERT_EQ(S(D("0").StripTrailingZeros()),       "0");
  ASSERT_EQ(S(D("1.5").StripTrailingZeros()),     "1.5");
}

// ─── comparison ──────────────────────────────────────────────────────────────

REGISTER_TEST(BigDecimal, CompareTo)
{
  ASSERT_TRUE(D("1.5") <  D("2"));
  ASSERT_TRUE(D("2")   >  D("1.5"));
  ASSERT_TRUE(D("1.5") == D("1.5"));
  ASSERT_TRUE(D("1.50") == D("1.5"));          // value equality across scales
  ASSERT_TRUE(D("0") == D("0.0"));
  ASSERT_TRUE(D("-1") <  D("0"));
  ASSERT_TRUE(D("-1.5") < D("-1.0"));
  ASSERT_FALSE(D("1") < D("1"));
}

REGISTER_TEST(BigDecimal, EqualsStrict)
{
  // value-equal but scale-differ ⇒ Equals() is strict false
  ASSERT_FALSE(D("1.50").Equals(D("1.5")));
  ASSERT_TRUE(D("1.50").Equals(D("1.50")));
  // value-equal with same scale ⇒ true
  ASSERT_TRUE(D("0.0").Equals(D("0.0")));
  ASSERT_FALSE(D("0").Equals(D("0.0")));     // different scales: 0 vs 1
}

REGISTER_TEST(BigDecimal, SignumAndZero)
{
  ASSERT_EQ(D("0").Signum(), 0);
  ASSERT_EQ(D("0.000").Signum(), 0);
  ASSERT_EQ(D("3.14").Signum(), 1);
  ASSERT_EQ(D("-3.14").Signum(), -1);
  ASSERT_TRUE(D("0").IsZero());
  ASSERT_TRUE(D("0.0000").IsZero());
  ASSERT_FALSE(D("0.0001").IsZero());
}

REGISTER_TEST(BigDecimal, NegateAbs)
{
  ASSERT_EQ(S(D("3.14").Negate()),  "-3.14");
  ASSERT_EQ(S(D("-3.14").Negate()), "3.14");
  ASSERT_EQ(S(D("0").Negate()),     "0");
  ASSERT_EQ(S(D("-3.14").Abs()),    "3.14");
  ASSERT_EQ(S(D("3.14").Abs()),     "3.14");
}

REGISTER_TEST(BigDecimal, OperatorOverloads)
{
  BigDecimal a = D("1.5");
  BigDecimal b = D("2.5");
  ASSERT_EQ(S(a + b), "4.0");
  ASSERT_EQ(S(b - a), "1.0");
  ASSERT_EQ(S(a * b), "3.75");
  ASSERT_EQ(S(-a),    "-1.5");
}

REGISTER_TEST(BigDecimal, Precision)
{
  ASSERT_EQ(D("0").Precision(), 1);
  ASSERT_EQ(D("9").Precision(), 1);
  ASSERT_EQ(D("10").Precision(), 2);
  ASSERT_EQ(D("123.456").Precision(), 6);
  ASSERT_EQ(D("-123").Precision(), 3);
}

// ─── large-scale sanity ──────────────────────────────────────────────────────

REGISTER_TEST(BigDecimal, BigPrecisionDivide)
{
  // pi-like: 22/7 to 50 digits, HALF_EVEN
  BigDecimal q = D("22").Divide(D("7"), 50, RoundingMode::HALF_EVEN);
  // 22/7 = 3.14285714285714285714... (period 6)
  std::string s = S(q);
  ASSERT_TRUE(s.find("3.142857142857142857142857142857142857142857142857") != std::string::npos);
  ASSERT_EQ((int)s.size(), 52);   // "3." + 50 digits
}

REGISTER_TEST(BigDecimal, LargeNumberArithmetic)
{
  BigDecimal a = D("123456789012345678901234567890.12345678901234567890");
  BigDecimal b = D("987654321098765432109876543210.98765432109876543210");
  BigDecimal sum = a + b;
  ASSERT_EQ(S(sum), "1111111110111111111011111111101.11111111011111111100");
  BigDecimal diff = b - a;
  ASSERT_EQ(S(diff), "864197532086419753208641975320.86419753208641975320");
}
