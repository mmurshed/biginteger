// Quotient-only division (Divide / DivideAndRemainder(..., false)).
//
// When the quotient is shorter than the divisor the dispatcher computes it
// from the operands' top limbs and only falls back to a full-width check
// when the top remainder is too small to decide. These tests aim at exactly
// that boundary: dividends built as q*b + r with r = 0, 1, b - 1 and
// q*b - 1, where the truncated quotient is one too large or undecidable.

#include "unit_test_framework.h"

#include <cstdint>
#include <random>
#include <vector>

#include "biginteger/algorithms/Addition.h"
#include "biginteger/algorithms/Division.h"
#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/algorithms/Subtraction.h"
#include "biginteger/common/Comparator.h"

using namespace BigMath;

#if BIGMATH_LIMB_64

namespace
{
  const DataT M = ~(DataT)0;
  const std::vector<DataT> kOne{1};

  std::vector<DataT> Trimmed(std::vector<DataT> v)
  {
    while (v.size() > 1 && v.back() == 0)
      v.pop_back();
    return v;
  }

  std::vector<DataT> Random(SizeT n, std::mt19937_64 &gen)
  {
    std::vector<DataT> v(n);
    for (auto &x : v) x = gen();
    if (v.back() == 0) v.back() = 1;
    return v;
  }

  // Divisor patterns that stress the truncation: random, smallest and
  // largest top limb, a power of the base, base power plus one, all ones.
  std::vector<DataT> Divisor(SizeT n, int kind, std::mt19937_64 &gen)
  {
    std::vector<DataT> b = Random(n, gen);
    switch (kind)
    {
    case 1: b.back() = 1; break;
    case 2: b.back() = M; break;
    case 3: for (SizeT i = 0; i + 1 < n; ++i) b[i] = 0; b.back() = 1; break;
    case 4: for (SizeT i = 0; i + 1 < n; ++i) b[i] = 0; b.back() = 1; b[0] = 1; break;
    case 5: for (auto &x : b) x = M; break;
    default: break;
    }
    return b;
  }

  void ExpectQuotient(std::vector<DataT> const &a, std::vector<DataT> const &b, std::vector<DataT> const &q)
  {
    ASSERT_TRUE(Trimmed(Divide(a, b, Base2_64)) == Trimmed(q));
    ASSERT_TRUE(Trimmed(DivideAndRemainder(a, b, Base2_64, false).first) == Trimmed(q));
  }

  // All dividends q*b + r around one (q, b) pair.
  void CheckAround(std::vector<DataT> const &q, std::vector<DataT> const &b)
  {
    const std::vector<DataT> qb = Trimmed(Multiply(q, b, Base2_64));
    ExpectQuotient(qb, b, q);                                                 // r = 0
    ExpectQuotient(Trimmed(Add(qb, kOne, Base2_64)), b, q);                   // r = 1
    const std::vector<DataT> bm1 = Trimmed(Subtract(b, kOne, Base2_64));
    ExpectQuotient(Trimmed(Add(qb, bm1, Base2_64)), b, q);                    // r = b - 1
    const std::vector<DataT> below = Trimmed(Subtract(qb, kOne, Base2_64));   // q*b - 1
    if (Compare(below, b) > 0)
      ExpectQuotient(below, b, Trimmed(Subtract(q, kOne, Base2_64)));
  }
}

REGISTER_TEST(DivisionQuotientOnly, ExactAndNearMultiples)
{
  std::mt19937_64 gen(0x0707A11ULL);
  for (SizeT nb : {5, 6, 7, 12, 40, 95, 96, 100, 130, 260})
    for (SizeT ql : {(SizeT)1, (SizeT)2, (SizeT)3, nb / 2, nb - 5, nb - 4, nb - 3, nb})
    {
      if (ql == 0 || ql > nb)
        continue;
      for (int kind = 0; kind < 6; ++kind)
      {
        const std::vector<DataT> b = Divisor(nb, kind, gen);
        CheckAround(Random(ql, gen), b);
        CheckAround(std::vector<DataT>(ql, M), b); // quotient digits all B - 1
        std::vector<DataT> pow(ql, 0);             // quotient a power of the base
        pow.back() = 1;
        CheckAround(pow, b);
      }
    }
}

REGISTER_TEST(DivisionQuotientOnly, AgreesWithFullDivision)
{
  std::mt19937_64 gen(0xA6EE5ULL);
  const long before = NewtonDivision::FallbackCount();
  for (SizeT nb : {5, 9, 33, 97, 200, 700})
    for (SizeT delta : {(SizeT)0, (SizeT)1, (SizeT)2, nb / 3, nb - 6, nb - 5, nb - 4, nb - 1})
      for (int kind = 0; kind < 6; ++kind)
      {
        const std::vector<DataT> b = Divisor(nb, kind, gen);
        std::vector<DataT> a = Random(nb + delta, gen);
        if (kind == 5)
          for (auto &x : a) x = M;
        ExpectQuotient(a, b, DivideAndRemainder(a, b, Base2_64, true).first);
      }
  ASSERT_EQ(NewtonDivision::FallbackCount(), before);
}

// The top division itself goes through the dispatcher: long enough quotients
// reach Burnikel–Ziegler and Newton there.
REGISTER_TEST(DivisionQuotientOnly, LongQuotients)
{
  std::mt19937_64 gen(0x10A6ULL);
  const long before = NewtonDivision::FallbackCount();
  for (auto [nb, ql] : {std::pair<SizeT, SizeT>{700, 300}, {2500, 800}, {4000, 1800}, {8000, 3300}})
  {
    const std::vector<DataT> b = Random(nb, gen);
    CheckAround(Random(ql, gen), b);
    const std::vector<DataT> a = Random(nb + ql - 1, gen);
    ExpectQuotient(a, b, DivideAndRemainder(a, b, Base2_64, true).first);
  }
  ASSERT_EQ(NewtonDivision::FallbackCount(), before);
}

#endif
