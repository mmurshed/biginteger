// Toom-3 and Toom-4 splits inside the Karatsuba recursion (64-bit limbs).
//
// Products and squares are checked against the schoolbook kernel on
// operands that exercise each branch of the evaluation: random limbs, all
// ones (largest evaluations), sparse limbs, and pieces arranged so that the
// value at -1 or -2 is exactly zero or negative. One forced Toom level is
// tested directly through MultiplyToomLevel, so Toom-4 is covered even
// though its default threshold leaves it off.

#include "unit_test_framework.h"

#include <cstdint>
#include <random>
#include <vector>

#include "biginteger/algorithms/multiplication/KaratsubaMultiplication.h"
#include "biginteger/algorithms/multiplication/KaratsubaSquare.h"
#include "biginteger/common/LimbKernels.h"

using namespace BigMath;

#if BIGMATH_LIMB_64

namespace
{
  const DataT M = ~(DataT)0;

  std::vector<DataT> Reference(std::vector<DataT> const &a, std::vector<DataT> const &b)
  {
    std::vector<DataT> r(a.size() + b.size());
    LimbKernel::MulBasecase(r.data(), a.data(), (SizeT)a.size(), b.data(), (SizeT)b.size());
    while (r.size() > 1 && r.back() == 0)
      r.pop_back();
    return r;
  }

  // kind 0 random, 1 all ones, 2 sparse, 3 pieces chosen for x(-1) = 0,
  // 4 large odd pieces / small even pieces (x(-1), x(-2) negative).
  std::vector<DataT> Operand(SizeT n, int kind, SizeT pieces, std::mt19937_64 &gen)
  {
    std::vector<DataT> v(n);
    const SizeT k = (n + pieces - 1) / pieces;
    for (SizeT i = 0; i < n; ++i)
    {
      const SizeT piece = i / k;
      switch (kind)
      {
      case 1: v[i] = M; break;
      case 2: v[i] = (gen() % 5 == 0) ? gen() : 0; break;
      case 3: v[i] = (piece == 0 || piece == 1) ? 7 : 0; break; // x0 = x1, x2 = x3 = 0 (but top nonzero below)
      case 4: v[i] = (piece % 2) ? M : 1; break;
      default: v[i] = gen(); break;
      }
    }
    if (kind == 3)
    {
      // x0 + x2 == x1 + x3 for 4 pieces: x0 = x1 = 7..., x2 = x3 = 0 except
      // the top piece must be nonzero, so put 1 in the top limb of x3 and
      // x2's matching limb.
      if (pieces == 4 && 2 * k + 0 < n)
        v[2 * k] = 1;
    }
    v.back() |= 1;
    return v;
  }
}

REGISTER_TEST(Toom, ForcedLevelMatchesSchoolbook)
{
  std::mt19937_64 gen(0x7003ULL);
  for (int way : {3, 4})
    for (SizeT la : {(SizeT)9, (SizeT)12, (SizeT)13, (SizeT)25, (SizeT)40, (SizeT)97, (SizeT)200, (SizeT)301})
      for (SizeT lb : {la, la - la / 8, la - la / 5})
        for (int kind = 0; kind < 5; ++kind)
        {
          const SizeT k = (std::max(la, lb) + way - 1) / way;
          if (la <= (SizeT)(way - 1) * k || lb <= (SizeT)(way - 1) * k)
            continue;
          auto a = Operand(la, kind, way, gen), b = Operand(lb, (kind + way) % 5, way, gen);
          ASSERT_TRUE(KaratsubaMultiplication::MultiplyToomLevel(a, b, way) == Reference(a, b));
          ASSERT_TRUE(KaratsubaMultiplication::MultiplyToomLevel(a, a, way) == Reference(a, a));
        }
}

REGISTER_TEST(Toom, RecursionAroundThresholds)
{
  std::mt19937_64 gen(0x7004ULL);
  for (SizeT la : {(SizeT)160, (SizeT)167, (SizeT)168, (SizeT)169, (SizeT)252, (SizeT)300, (SizeT)505, (SizeT)700})
    for (SizeT lb : {la, la - la / 6, la - la / 3 + 1, la / 2 + 1, la * 3})
      for (int kind = 0; kind < 5; ++kind)
      {
        auto a = Operand(la, kind, 3, gen), b = Operand(lb, (kind + 1) % 5, 3, gen);
        auto c = KaratsubaMultiplication::Multiply(a, b, Base2_64);
        ASSERT_TRUE(c == Reference(a, b));
      }
}

REGISTER_TEST(Toom, SquaresAroundThreshold)
{
  std::mt19937_64 gen(0x7005ULL);
  for (SizeT n : {(SizeT)120, (SizeT)127, (SizeT)128, (SizeT)129, (SizeT)130, (SizeT)191, (SizeT)384, (SizeT)385, (SizeT)607})
    for (int kind = 0; kind < 5; ++kind)
    {
      auto a = Operand(n, kind, 3, gen);
      ASSERT_TRUE(KaratsubaSquare::Square(a, Base2_64) == Reference(a, a));
    }
}

#endif
