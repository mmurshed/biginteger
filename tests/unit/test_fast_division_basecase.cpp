// FastDivision's 64-bit-limb basecase (raw-buffer Knuth D on the SubMul1
// kernel) against the reference KnuthDivision and the division identity.
//
// Shapes sweep the block remainders of the kernels, both normalization cases
// (shift == 0 and shift > 0), single- and two-limb divisors, and structured
// operands that force the qhat == B-1 special case and the add-back step.

#include "unit_test_framework.h"

#include <cstdint>
#include <random>
#include <vector>

#include "biginteger/algorithms/Addition.h"
#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/algorithms/division/FastDivision.h"
#include "biginteger/algorithms/division/KnuthDivision.h"
#include "biginteger/common/Comparator.h"

using namespace BigMath;

namespace
{
  std::vector<DataT> Trimmed(std::vector<DataT> v)
  {
    while (v.size() > 1 && v.back() == 0)
      v.pop_back();
    return v;
  }

  // q*b + r == a and r < b, plus agreement with KnuthDivision.
  void CheckDivision(std::vector<DataT> const &a, std::vector<DataT> const &b)
  {
    auto qr = FastDivision::DivideAndRemainder(a, b, Base2_64, true);
    auto q = Trimmed(qr.first);
    auto r = Trimmed(qr.second);
    ASSERT_TRUE(Compare(r, b) < 0);
    auto back = Trimmed(Add(Multiply(q, b, Base2_64), r, Base2_64));
    ASSERT_TRUE(back == Trimmed(a));

    auto ref = KnuthDivision::DivideAndRemainder(a, b, Base2_64);
    ASSERT_TRUE(q == Trimmed(ref.first));
    ASSERT_TRUE(r == Trimmed(ref.second));

    auto qOnly = FastDivision::DivideAndRemainder(a, b, Base2_64, false);
    ASSERT_TRUE(Trimmed(qOnly.first) == q);
  }

  std::vector<DataT> Random(SizeT n, std::mt19937_64 &gen)
  {
    std::vector<DataT> v(n);
    for (auto &x : v)
      x = gen();
    if (v.back() == 0)
      v.back() = 1;
    return v;
  }
}

#if BIGMATH_LIMB_64

REGISTER_TEST(FastDivisionBasecase, RandomShapes)
{
  std::mt19937_64 gen(0xFA57D1FULL);
  for (SizeT nb : {1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100, 130})
    for (SizeT extra : {0, 1, 2, 3, 5, 8, 17, 64, 191, 200})
    {
      auto a = Random(nb + extra, gen);
      auto b = Random(nb, gen);
      if (Compare(a, b) <= 0)
        a.back() |= (DataT)1 << 63, a.push_back(1);
      CheckDivision(a, b);
      // Divisor with its top bit already set: normalization shift of zero.
      b.back() |= (DataT)1 << 63;
      if (Compare(a, b) <= 0)
        a.push_back(1);
      CheckDivision(a, b);
      // Divisor with a tiny top limb: the largest shift.
      b.back() = 1;
      CheckDivision(a, b);
    }
}

REGISTER_TEST(FastDivisionBasecase, StructuredOperands)
{
  const DataT M = ~(DataT)0;
  std::mt19937_64 gen(0x57AC7ULL);
  for (SizeT nb : {2, 3, 4, 8, 16, 33, 64})
    for (SizeT na : {nb + 1, nb + 2, 2 * nb, 2 * nb + 1, 3 * nb + 3})
    {
      // All-ones dividend and divisor: maximal quotient digits.
      CheckDivision(std::vector<DataT>(na, M), std::vector<DataT>(nb, M));
      // Divisor B^(nb-1) * (B/2) exactly, and one below a power of the base.
      std::vector<DataT> b(nb, 0);
      b.back() = (DataT)1 << 63;
      CheckDivision(std::vector<DataT>(na, M), b);
      b.assign(nb, M);
      b[0] = 0;
      CheckDivision(std::vector<DataT>(na, M), b);
      // Remainder windows whose top two limbs equal the divisor's top two
      // (the qhat == B-1 special case): dividend = divisor shifted up, minus one.
      std::vector<DataT> d = Random(nb, gen);
      d.back() |= (DataT)1 << 63;
      std::vector<DataT> a(na, 0);
      for (SizeT i = 0; i < nb; ++i)
        a[na - nb + i] = d[i];
      // subtract one: borrow through the zero limbs
      for (SizeT i = 0; i < na; ++i)
      {
        if (a[i]-- != 0)
          break;
      }
      CheckDivision(a, d);
      // Sparse low half in both operands.
      std::vector<DataT> a2 = Random(na, gen), b2 = Random(nb, gen);
      for (SizeT i = 0; i < na / 2; ++i)
        a2[i] = 0;
      for (SizeT i = 0; i < nb / 2; ++i)
        b2[i] = 0;
      if (Compare(a2, b2) > 0)
        CheckDivision(a2, b2);
    }
}

#endif
