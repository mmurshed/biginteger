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
