#include "unit_test_framework.h"

#include <stdexcept>

#include "biginteger/common/Builder.h"
#include "biginteger/ops/Operations.h"

using namespace BigMath;

REGISTER_TEST(NumberTheory, PowMod)
{
  const auto value = [](const char *s) { return BigIntegerBuilder::From(s); };
  ASSERT_TRUE(PowMod(value("2"), value("1000000"), value("1000000007")) == value("235042059"));
  ASSERT_TRUE(PowMod(value("-3"), value("7"), value("100")) == value("13"));
  ASSERT_TRUE(PowMod(value("0"), value("0"), value("17")) == value("1"));
  ASSERT_TRUE(PowMod(value("0"), value("0"), value("1")).Zero());
  ASSERT_TRUE(PowMod(value("-30"), value("9"), value("15")).Zero());
  ASSERT_TRUE(PowMod(value("3"), value("340282366920938463463374607431768211457"),
                     value("1000000007")) == value("800740594"));
  ASSERT_TRUE(PowMod(value("2"), BigInteger(12), value("17")) == value("1"));

  for (const auto &[exponent, modulus] : {std::pair{"-1", "5"}, {"1", "0"}, {"1", "-5"}})
  {
    bool threw = false;
    try { (void)PowMod(value("2"), value(exponent), value(modulus)); }
    catch (std::invalid_argument const &) { threw = true; }
    ASSERT_TRUE(threw);
  }
}

REGISTER_TEST(NumberTheory, GreatestCommonDivisors)
{
  const auto value = [](const char *s) { return BigIntegerBuilder::From(s); };
  const auto check = [](BigInteger const &a, BigInteger const &b, BigInteger const &expected) {
    ASSERT_TRUE(Gcd(a, b) == expected);
    const auto [divisor, x, y] = ExtendedGcd(a, b);
    ASSERT_TRUE(divisor == expected);
    ASSERT_TRUE(a * x + b * y == divisor);
  };

  check(BigInteger(12), BigInteger(), BigInteger());
  check(value("-18"), value("0"), value("18"));
  check(value("0"), value("-18"), value("18"));
  check(value("19"), value("35"), value("1"));

  const auto factor = value("340282366920938463463374607431768211457");
  const auto left = factor * value("462");
  const auto right = factor * value("1071");
  const auto expected = factor * value("21");
  for (bool negativeLeft : {false, true})
    for (bool negativeRight : {false, true})
      check(BigInteger(left).SetSign(negativeLeft), BigInteger(right).SetSign(negativeRight), expected);
  check(right, left, expected);
}

REGISTER_TEST(NumberTheory, SquareRootAndRemainder)
{
  const auto value = [](const char *s) { return BigIntegerBuilder::From(s); };
  const auto check = [](BigInteger const &n, BigInteger const &expectedRoot,
                        BigInteger const &expectedRemainder) {
    const auto [root, remainder] = SqrtRem(n);
    ASSERT_TRUE(root == expectedRoot);
    ASSERT_TRUE(remainder == expectedRemainder);
  };
  check(BigInteger(12), value("0"), value("0"));
  check(value("1"), value("1"), value("0"));
  check(value("3"), value("1"), value("2"));
  check(value("18446744073709551615"), value("4294967295"), value("8589934590"));

  const auto one = value("1");
  const auto root = value("340282366920938463463374607431768211456");
  const auto square = root * root;
  check(square - one, root - one, root + root - value("2"));
  check(square, root, value("0"));
  check(square + one, root, one);

  bool threw = false;
  try { (void)SqrtRem(value("-1")); }
  catch (std::invalid_argument const &) { threw = true; }
  ASSERT_TRUE(threw);
}
