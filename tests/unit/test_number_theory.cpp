#include "unit_test_framework.h"

#include <stdexcept>
#include <thread>

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

  const auto base = value("-12345678901234567890123456789012345");
  const auto sparse = value("340282366920938463481821351505477828609");
  ASSERT_TRUE(PowMod(base, sparse, value("340282366920938463463374607431768211455")) ==
              value("293408806814659110257804787314180246610"));
  ASSERT_TRUE(PowMod(base, sparse, value("170141183460469231731687303715884105824")) ==
              value("47847993867134196767325167157259808199"));
  ASSERT_TRUE(PowMod(base, sparse, value("42535295865117308763025309245900849152")) ==
              value("28939521772930582829908416378337304711"));
  ASSERT_TRUE(PowMod(base, sparse, value("680564733841876926926749214863536422912")) ==
              value("481524715727764675644108224933147254919"));

  for (const auto &[exponent, modulus] : {std::pair{"-1", "5"}, {"1", "0"}, {"1", "-5"}})
  {
    bool threw = false;
    try { (void)PowMod(value("2"), value(exponent), value(modulus)); }
    catch (std::invalid_argument const &) { threw = true; }
    ASSERT_TRUE(threw);
  }
}

REGISTER_TEST(NumberTheory, Primality)
{
  const auto value = [](const char *s) { return BigIntegerBuilder::From(s); };
  for (const char *prime : {"2", "3", "53", "4294967291", "18446744073709551557",
                            "18446744073709551629", "461168601842738790401",
                            "170141183460469231731687303715884105727"})
    ASSERT_TRUE(IsProbablePrime(value(prime)));
  for (const char *composite : {"-3", "0", "1", "4", "561", "2047",
                                "3215031751", "341550071728321", "3825123056546413051",
                                "18446744073709551615", "18446744073709551617",
                                "318665857834031151167461"})
    ASSERT_TRUE(!IsProbablePrime(value(composite)));
  const auto largePrime = value("170141183460469231731687303715884105727");
  ASSERT_TRUE(!IsProbablePrime(largePrime * largePrime));
  ASSERT_TRUE(!IsProbablePrime(BigInteger(12)));

  bool threw = false;
  try { (void)IsProbablePrime(value("7"), 0); }
  catch (std::invalid_argument const &) { threw = true; }
  ASSERT_TRUE(threw);
}

REGISTER_TEST(NumberTheory, MultiBatchAndConcurrentPrimality)
{
  // The Mersenne prime 2^2203 - 1 reaches the threaded witness path in
  // both limb layouts. Thirty-three rounds cross the batch boundary.
  constexpr DataT high = DataT{1} << (2203 % LimbBits);
  const auto one = BigIntegerBuilder::From("1");
  const auto prime = (BigInteger::FromLimbs({&high, 1}, false) << (2203 / LimbBits)) - one;
  ASSERT_TRUE(IsProbablePrime(prime, 33));

  bool first = false, second = false;
  std::thread workerA([&] { first = IsProbablePrime(prime, 8); });
  std::thread workerB([&] { second = IsProbablePrime(prime, 8); });
  workerA.join();
  workerB.join();
  ASSERT_TRUE(first && second);
}

REGISTER_TEST(NumberTheory, NextProbablePrime)
{
  const auto value = [](const char *s) { return BigIntegerBuilder::From(s); };
  for (const auto &[input, expected] : {
           std::pair{"-100", "2"}, {"0", "2"}, {"1", "2"}, {"2", "3"},
           {"100", "101"}, {"101", "103"}, {"4093", "4099"},
           {"18446744073709551557", "18446744073709551629"},
           {"340282366920938463463374607431768211455",
            "340282366920938463463374607431768211507"}})
    ASSERT_TRUE(NextProbablePrime(value(input)) == value(expected));

  bool threw = false;
  try { (void)NextProbablePrime(value("0"), 0); }
  catch (std::invalid_argument const &) { threw = true; }
  ASSERT_TRUE(threw);

  // In this search the first sieve survivor is composite. The primes at
  // offsets 1493 and 1583 both fit in the first batch on an eight-thread build.
  constexpr DataT high = DataT{1} << (1023 % LimbBits);
  const auto base = BigInteger::FromLimbs({&high, 1}, false) << (1023 / LimbBits);
  const auto start = base + value("1445");
  const auto lower = base + value("1493");
  const auto higher = base + value("1583");
  ASSERT_TRUE(NextProbablePrime(start) == lower);
  ASSERT_TRUE(NextProbablePrime(lower) == higher);
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

  // Fibonacci pairs give long runs of quotient one. A large shared shift
  // puts the nonzero common divisor well above the low limbs.
  auto previous = value("1");
  auto current = value("1");
  for (unsigned i = 0; i < 350; ++i)
  {
    auto next = previous + current;
    previous = std::move(current);
    current = std::move(next);
  }
  const auto shared = value("1") << 5;
  check(-(current << 5), previous << 5, shared);

  // Equal limb counts can still yield a quotient too large for a batch.
  const auto divisor = (value("1") << 4) + value("1");
  check(divisor * value("2147483648") + value("1"), divisor, value("1"));
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

  // Exercise the native-word boundary and the first divide-and-conquer
  // splits. The inputs just below a square require an exact correction.
  for (unsigned bits : {32, 63, 64, 65, 95, 96, 126, 127, 128, 129, 192, 256})
  {
    const DataT bit = DataT{1} << (bits % LimbBits);
    const auto boundary = BigInteger::FromLimbs({&bit, 1}, false) << (bits / LimbBits);
    const auto exact = boundary * boundary;
    check(exact - one, boundary - one, boundary + boundary - value("2"));
    check(exact, boundary, value("0"));
    check(exact + one, boundary, one);
  }

  bool threw = false;
  try { (void)SqrtRem(value("-1")); }
  catch (std::invalid_argument const &) { threw = true; }
  ASSERT_TRUE(threw);
}
