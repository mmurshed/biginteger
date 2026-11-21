// BigInteger small-buffer storage and the small-operand operator paths.
//
// Values of up to BigInteger::InlineLimbs limbs live inside the object;
// longer ones on the heap. The operators combine small operands in a stack
// buffer (64-bit limbs only). These tests pin the representation at the
// boundary, copy/move behaviour, and check every operator against the
// vector-based dispatchers on operands straddling the boundary, for all sign
// combinations.

#include "unit_test_framework.h"

#include <climits>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#include "biginteger/BigInteger.h"
#include "biginteger/algorithms/Addition.h"
#include "biginteger/algorithms/Division.h"
#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/algorithms/Subtraction.h"
#include "biginteger/common/Builder.h"
#include "biginteger/ops/Operations.h"

using namespace BigMath;

namespace
{
  constexpr SizeT K = BigInteger::InlineLimbs;

  std::vector<DataT> Trimmed(std::vector<DataT> v)
  {
    while (v.size() > 1 && v.back() == 0)
      v.pop_back();
    return v;
  }

  DataT RandomLimb(std::mt19937_64 &gen)
  {
#if BIGMATH_LIMB_64
    return (DataT)gen();
#else
    return (DataT)(gen() & 0xFFFFFFFFULL);
#endif
  }

  DataT MaxLimb()
  {
#if BIGMATH_LIMB_64
    return ~(DataT)0;
#else
    return (DataT)0xFFFFFFFFULL;
#endif
  }

  // n limbs, top limb nonzero; kind 0 random, 1 all ones, 2 power of the base.
  std::vector<DataT> Magnitude(SizeT n, int kind, std::mt19937_64 &gen)
  {
    std::vector<DataT> v(n, 0);
    for (auto &x : v)
      x = kind == 1 ? MaxLimb() : (kind == 2 ? 0 : RandomLimb(gen));
    if (kind == 2 || v.back() == 0)
      v.back() = 1;
    return v;
  }

  // Signed reference arithmetic on (magnitude, sign) pairs via the vector
  // dispatchers.
  struct Ref
  {
    std::vector<DataT> mag;
    bool neg;
  };

  Ref RefAdd(Ref const &a, Ref const &b)
  {
    const BaseT base = BigInteger::Base();
    if (a.neg == b.neg)
      return {Trimmed(Add(a.mag, b.mag, base)), a.neg};
    Int c = Compare(a.mag, b.mag);
    if (c == 0)
      return {{0}, false};
    if (c > 0)
      return {Trimmed(Subtract(a.mag, b.mag, base)), a.neg};
    return {Trimmed(Subtract(b.mag, a.mag, base)), b.neg};
  }

  void ExpectEqual(BigInteger const &got, Ref const &want)
  {
    const bool zero = Trimmed(want.mag) == std::vector<DataT>{0};
    ASSERT_TRUE(got.GetInteger() == Trimmed(want.mag));
    ASSERT_EQ(got.IsNegative(), want.neg && !zero);
    // Representation follows the length.
    ASSERT_EQ(got.IsInline(), got.size() <= K);
  }
}

REGISTER_TEST(SmallBuffer, RepresentationFollowsLength)
{
  std::mt19937_64 gen(7);
  for (SizeT n = 1; n <= K + 3; ++n)
  {
    auto v = Magnitude(n, 0, gen);
    BigInteger a(v, false);
    ASSERT_EQ(a.size(), n);
    ASSERT_EQ(a.IsInline(), n <= K);
    ASSERT_TRUE(a.GetInteger() == v);
    auto l = a.Limbs();
    ASSERT_TRUE(std::vector<DataT>(l.begin(), l.end()) == v);
    ASSERT_TRUE(a.AsVector().get() == v);
    for (SizeT i = 0; i < n; ++i)
      ASSERT_EQ(a[i], v[i]);

    // Leading zeros are trimmed before choosing the representation.
    auto padded = v;
    padded.resize(K + 5, 0);
    BigInteger b(padded, true);
    ASSERT_EQ(b.IsInline(), n <= K);
    ASSERT_TRUE(b.GetInteger() == v);
    BigInteger c(std::move(padded), true);
    ASSERT_EQ(c.IsInline(), n <= K);
    ASSERT_TRUE(c.GetInteger() == v);
    BigInteger d = BigInteger::FromLimbs(std::span<const DataT>(b.Limbs()), false);
    ASSERT_EQ(d.IsInline(), n <= K);
    ASSERT_TRUE(d == BigInteger(v, false));
  }

  // The zero-of-size constructor keeps its length, inline or not.
  BigInteger z(K + 2);
  ASSERT_EQ(z.size(), K + 2);
  ASSERT_TRUE(!z.IsInline());
  ASSERT_TRUE(z.Zero());
  z.Trim();
  ASSERT_EQ(z.size(), 1);
  ASSERT_TRUE(z.IsInline());
  BigInteger negZero(3, true);
  ASSERT_TRUE(!negZero.IsNegative());
}

REGISTER_TEST(SmallBuffer, CopyMoveAndRelease)
{
  std::mt19937_64 gen(8);
  for (SizeT n : {(SizeT)1, K, K + 1, (SizeT)40})
  {
    auto v = Magnitude(n, 0, gen);
    BigInteger a(v, true);

    BigInteger copy(a);
    ASSERT_TRUE(copy == a && copy.IsNegative());
    BigInteger assigned;
    assigned = a;
    ASSERT_TRUE(assigned == a);

    BigInteger moved(std::move(copy));
    ASSERT_TRUE(moved == a);
    // Moved-from values read as (positive) zero and stay usable.
    ASSERT_TRUE(copy.Zero() && !copy.IsNegative() && copy.size() == 1);
    copy = moved;
    ASSERT_TRUE(copy == a);

    BigInteger moveAssigned;
    moveAssigned = std::move(assigned);
    ASSERT_TRUE(moveAssigned == a);
    ASSERT_TRUE(assigned.Zero() && !assigned.IsNegative());
    moveAssigned = std::move(moveAssigned); // self-move leaves the value alone
    ASSERT_TRUE(moveAssigned == a);

    // Assigning across representations in both directions.
    BigInteger small(std::vector<DataT>{5}, false), large(Magnitude(K + 4, 0, gen), false);
    BigInteger x = small;
    x = large;
    ASSERT_TRUE(x == large && !x.IsInline());
    x = small;
    ASSERT_TRUE(x == small && x.IsInline());
    x = std::move(large);
    ASSERT_TRUE(!x.IsInline() && large.Zero());

    BigInteger r(a);
    std::vector<DataT> released = r.ReleaseInteger();
    ASSERT_TRUE(released == v);
    ASSERT_TRUE(r.Zero() && !r.IsNegative());
  }
}

REGISTER_TEST(SmallBuffer, OperatorsMatchVectorDispatchers)
{
  std::mt19937_64 gen(0x5B0ULL);
  const BaseT base = BigInteger::Base();
  // Inline values, just past them, and both sides of the 16-limb stack path.
  std::vector<SizeT> lengths;
  for (SizeT n = 1; n <= K + 2; ++n)
    lengths.push_back(n);
  for (SizeT n : {(SizeT)15, (SizeT)16, (SizeT)17})
    lengths.push_back(n);

  for (SizeT na : lengths)
    for (SizeT nb : lengths)
      for (int kind = 0; kind < 3; ++kind)
        for (int signs = 0; signs < 4; ++signs)
        {
          const bool aNeg = signs & 1, bNeg = signs & 2;
          auto am = Magnitude(na, kind, gen);
          auto bm = Magnitude(nb, (kind + signs) % 3, gen);
          BigInteger a(am, aNeg), b(bm, bNeg);
          Ref ra{am, aNeg}, rb{bm, bNeg};

          ExpectEqual(a + b, RefAdd(ra, rb));
          ExpectEqual(a - b, RefAdd(ra, Ref{bm, !bNeg}));
          ExpectEqual(a * b, Ref{Trimmed(Multiply(am, bm, base)), aNeg != bNeg});
          ExpectEqual(a * a, Ref{Trimmed(Multiply(am, std::vector<DataT>(am), base)), false});

          auto qr = DivideAndRemainder(am, bm, base, true);
          Ref q{Trimmed(qr.first), aNeg != bNeg}, r{Trimmed(qr.second), aNeg};
          ExpectEqual(a / b, q);
          ExpectEqual(a % b, r);
          auto [q2, r2] = DivideAndRemainder(a, b);
          ExpectEqual(q2, q);
          ExpectEqual(r2, r);

          // x - x, x + (-x): zero results are positive and inline.
          ExpectEqual(a - a, Ref{{0}, false});
          ExpectEqual(a + (-a), Ref{{0}, false});
          ExpectEqual(a / a, Ref{{1}, false});
          ExpectEqual(a % a, Ref{{0}, false});
        }
}

REGISTER_TEST(SmallBuffer, DivisionEdgeCases)
{
  std::mt19937_64 gen(11);
  const BaseT base = BigInteger::Base();
  const DataT M = MaxLimb();
  // Quotient digits of B - 1, divisor top limbs 1 and B - 1, exact multiples,
  // and dividends just below a multiple.
  std::vector<std::vector<DataT>> divisors{{M}, {1, 1}, {M, M}, {0, 1}, {1, 0, 1}, {M, 0, M}, {3, M, 1}};
  for (auto const &bm : divisors)
    for (SizeT qn = 1; qn + bm.size() <= K + 1; ++qn)
    {
      auto qm = Magnitude(qn, 1, gen);
      auto prod = Trimmed(Multiply(qm, bm, base));
      std::vector<std::vector<DataT>> dividends{prod, Trimmed(Subtract(prod, std::vector<DataT>{1}, base)),
                                                Trimmed(Add(prod, Trimmed(Subtract(bm, std::vector<DataT>{1}, base)), base))};
      for (auto const &am : dividends)
      {
        BigInteger a(am, false), b(bm, false);
        auto qr = DivideAndRemainder(am, bm, base, true);
        ExpectEqual(a / b, Ref{Trimmed(qr.first), false});
        ExpectEqual(a % b, Ref{Trimmed(qr.second), false});
      }
    }

  bool threw = false;
  try
  {
    BigInteger one(std::vector<DataT>{1}, false);
    (void)(one / BigInteger());
  }
  catch (std::invalid_argument const &)
  {
    threw = true;
  }
  ASSERT_TRUE(threw);
}

REGISTER_TEST(SmallBuffer, ScalarOpsShiftsAndBuilder)
{
  std::mt19937_64 gen(12);
  const BaseT base = BigInteger::Base();
  for (SizeT n = 1; n <= K + 2; ++n)
    for (DataT d : {(DataT)1, (DataT)7, MaxLimb(), (DataT)(MaxLimb() >> 1) + 1})
      for (bool neg : {false, true})
      {
        auto am = Magnitude(n, 0, gen);
        BigInteger a(am, neg);
        std::vector<DataT> dv{d};
        ExpectEqual(a * d, Ref{Trimmed(Multiply(am, dv, base)), neg});
        auto qr = DivideAndRemainder(am, dv, base, true);
        ExpectEqual(a / d, Ref{Trimmed(qr.first), neg});
        ExpectEqual(a % d, Ref{Trimmed(qr.second), neg});
        auto [q, r] = DivideAndRemainder(a, d);
        ExpectEqual(q, Ref{Trimmed(qr.first), neg});
        ExpectEqual(r, Ref{Trimmed(qr.second), neg});

        for (SizeT s : {(SizeT)0, (SizeT)1, K - 1, K, K + 3})
        {
          std::vector<DataT> left(s, 0);
          left.insert(left.end(), am.begin(), am.end());
          ExpectEqual(a << s, Ref{left, neg});
          std::vector<DataT> right = s >= am.size() ? std::vector<DataT>{0}
                                                    : std::vector<DataT>(am.begin() + s, am.end());
          ExpectEqual(a >> s, Ref{right, neg});
        }
      }

  // Integers build directly, inline.
  for (Long v : {0LL, 1LL, -1LL, 42LL, -123456789012345LL, LLONG_MAX, LLONG_MIN})
  {
    BigInteger x = BigIntegerBuilder::From(v);
    ASSERT_TRUE(x.IsInline());
    ASSERT_TRUE(x == BigIntegerBuilder::From(std::to_string(v)));
  }
  for (ULong v : {0ULL, 1ULL, 0xFFFFFFFFULL, 0x100000000ULL, ULLONG_MAX})
  {
    BigInteger x = BigIntegerBuilder::From(v);
    ASSERT_TRUE(x.IsInline());
    ASSERT_TRUE(x == BigIntegerBuilder::From(std::to_string(v)));
  }
}
