#include "unit_test_framework.h"

#include <random>
#include <vector>

#include "biginteger/BigInteger.h"
#include "biginteger/common/Comparator.h"
#include "biginteger/ops/Division.h"
#include "biginteger/ops/Multiplication.h"

using namespace BigMath;

namespace
{
  std::vector<DataT> RandomLimbs(SizeT limbs, std::mt19937_64 &gen)
  {
    std::uniform_int_distribution<uint64_t> d(0, 0xFFFFFFFFULL);
    std::vector<DataT> v(limbs);
    for (auto &x : v)
      x = (DataT)d(gen);
    if (!v.empty() && v.back() == 0)
      v.back() = 1 + (d(gen) & 0xFFFF);
    return v;
  }

  BigInteger RandomInteger(SizeT limbs, std::mt19937_64 &gen, bool negative = false)
  {
    return BigInteger(RandomLimbs(limbs, gen), negative);
  }
}

REGISTER_TEST(PreparedAPI, ReusedMultiplyMatchesDispatcher)
{
  std::mt19937_64 gen(0xA111C0DE);
  BigInteger fixed = RandomInteger(1024, gen, false);
  PreparedMultiplication prepared(fixed, 1024);

  for (SizeT otherLimbs : {SizeT{128}, SizeT{256}, SizeT{512}, SizeT{1024}})
  {
    BigInteger other = RandomInteger(otherLimbs, gen, false);
    BigInteger expected = fixed * other;
    BigInteger actual = prepared * other;
    ASSERT_EQ(actual.CompareTo(expected), 0);
  }
}

REGISTER_TEST(PreparedAPI, ReusedMultiplyPreservesSign)
{
  std::mt19937_64 gen(0xA111C0DF);
  BigInteger fixed = RandomInteger(512, gen, true);
  PreparedMultiplication prepared(fixed, 512);

  for (bool negative : {false, true})
  {
    BigInteger other = RandomInteger(256, gen, negative);
    BigInteger expected = fixed * other;
    BigInteger actual = prepared.Multiply(other);
    ASSERT_EQ(actual.CompareTo(expected), 0);
  }
}

REGISTER_TEST(CachedAPI, ReusedDivisionMatchesDispatcher)
{
  std::mt19937_64 gen(0xBADC0DE);
  BigInteger divisor = RandomInteger(512, gen, false);
  CachedDivision cached(divisor);

  for (SizeT dividendLimbs : {SizeT{512}, SizeT{768}, SizeT{1024}, SizeT{2048}})
  {
    BigInteger dividend = RandomInteger(dividendLimbs, gen, false);
    auto expected = DivideAndRemainder(dividend, divisor);
    auto actual = cached.DivideAndRemainder(dividend);
    ASSERT_EQ(actual.first.CompareTo(expected.first), 0);
    ASSERT_EQ(actual.second.CompareTo(expected.second), 0);
  }
}

REGISTER_TEST(CachedAPI, ReusedDivisionPreservesSign)
{
  std::mt19937_64 gen(0xBADC0DF);
  BigInteger divisor = RandomInteger(256, gen, true);
  CachedDivision cached(divisor);

  for (bool negative : {false, true})
  {
    BigInteger dividend = RandomInteger(1024, gen, negative);
    auto expected = DivideAndRemainder(dividend, divisor);
    auto actual = cached.DivideAndRemainder(dividend);
    ASSERT_EQ(actual.first.CompareTo(expected.first), 0);
    ASSERT_EQ(actual.second.CompareTo(expected.second), 0);
  }
}

