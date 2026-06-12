// Seeded adversarial-limb fuzz: cross-checks the dispatch paths against
// O(n^2) references on operands biased toward 0 / 1 / 2^64-1 / 2^63 limbs.
//
// Why this exists: three silent wrong-answer bugs shipped in 2026 (MFA
// inverse twiddle #103, BZ odd-size fallback #116, FastDivision MG top-equal
// qhat) and none was catchable by uniformly random operands — each needed a
// structured limb pattern with probability ~0 under the random sweeps.
// Deterministic seeds keep CI reproducible; bump the per-test case counts
// locally for deeper runs.

#include "unit_test_framework.h"

#include <random>
#include <vector>

#include "biginteger/BigInteger.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Parser.h"
#include "biginteger/algorithms/Addition.h"
#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/algorithms/Squaring.h"
#include "biginteger/algorithms/Division.h"
#include "biginteger/algorithms/division/KnuthDivision.h"
#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"

using namespace BigMath;

namespace
{
  constexpr DataT kMax = ~(DataT)0;
  constexpr DataT kTop = (DataT)1 << 63;

  // ~25% all-ones, ~25% zero, then 1 / top-bit / uniform — the distribution
  // that surfaced the FastDivision MG bug.
  std::vector<DataT> AdversarialLimbs(SizeT n, std::mt19937_64 &rng)
  {
    std::vector<DataT> v(n);
    for (auto &x : v)
    {
      switch (rng() % 6)
      {
        case 0:
        case 1: x = kMax; break;
        case 2:
        case 3: x = 0; break;
        case 4: x = (rng() % 2) ? 1 : kTop; break;
        default: x = (DataT)rng(); break;
      }
    }
    while (v.size() > 1 && v.back() == 0)
      v.pop_back();
    if (IsZero(v))
      v[0] = 1;
    return v;
  }
}

// ── division: dispatch vs KnuthDivision reference + q*b + r == a ───────────
// Sizes span the FastDivision band and the small ends of the BZ/Newton
// bands; the Knuth O(n^2) reference keeps the runtime bounded.

REGISTER_TEST(AdversarialFuzz, DivisionDispatchVsKnuth)
{
  std::mt19937_64 rng(0xD1F1D1F1ULL);
  for (int t = 0; t < 2000; ++t)
  {
    SizeT na = 2 + (SizeT)(rng() % 40);
    SizeT nb = 1 + (SizeT)(rng() % na);
    auto a = AdversarialLimbs(na, rng);
    auto b = AdversarialLimbs(nb, rng);

    auto [q, r] = DivideAndRemainder(a, b, BigInteger::Base(), true);
    auto k = KnuthDivision::DivideAndRemainder(a, b, BigInteger::Base());
    ASSERT_EQ(Compare(q, k.first), 0);
    ASSERT_EQ(Compare(r, k.second), 0);

    auto back = Add(Multiply(q, b, BigInteger::Base()), r, BigInteger::Base());
    ASSERT_EQ(Compare(back, a), 0);
    ASSERT_EQ(Compare(r, b) < 0, true); // r < b
  }
}

REGISTER_TEST(AdversarialFuzz, DivisionLargeBandsIdentity)
{
  // Larger shapes hit BZ / Newton / QSized; verify via the identity only
  // (Knuth would dominate the runtime here).
  std::mt19937_64 rng(0xB16B00B5ULL);
  const SizeT shapes[][2] = {
      {700, 600},   // BZ near-balanced
      {1300, 640},  // BZ 2:1
      {1025, 513},  // odd sizes (BZ padding path)
      {2100, 600},  // big-and-skewed BZ band
      {1500, 160},  // FastDivision long skew
  };
  for (auto const &sh : shapes)
  {
    for (int t = 0; t < 3; ++t)
    {
      auto a = AdversarialLimbs(sh[0], rng);
      auto b = AdversarialLimbs(sh[1], rng);
      auto [q, r] = DivideAndRemainder(a, b, BigInteger::Base(), true);
      auto back = Add(Multiply(q, b, BigInteger::Base()), r, BigInteger::Base());
      ASSERT_EQ(Compare(back, a), 0);
      ASSERT_EQ(Compare(r, b) < 0, true);
    }
  }
}

// ── multiplication: dispatch vs classic schoolbook ──────────────────────────

REGISTER_TEST(AdversarialFuzz, MultiplicationDispatchVsClassic)
{
  std::mt19937_64 rng(0x5EEDBEEFULL);
  for (int t = 0; t < 1000; ++t)
  {
    SizeT na = 1 + (SizeT)(rng() % 96);
    SizeT nb = 1 + (SizeT)(rng() % 96);
    auto a = AdversarialLimbs(na, rng);
    auto b = AdversarialLimbs(nb, rng);
    auto fast = Multiply(a, b, BigInteger::Base());
    auto slow = ClassicMultiplication::Multiply(a, b, BigInteger::Base());
    ASSERT_EQ(Compare(fast, slow), 0);
  }
  // A few NTT-band cases (sum >= 1280 limbs) against classic.
  for (int t = 0; t < 3; ++t)
  {
    auto a = AdversarialLimbs(900, rng);
    auto b = AdversarialLimbs(700, rng);
    auto fast = Multiply(a, b, BigInteger::Base());
    auto slow = ClassicMultiplication::Multiply(a, b, BigInteger::Base());
    ASSERT_EQ(Compare(fast, slow), 0);
  }
}

REGISTER_TEST(AdversarialFuzz, SquareMatchesMultiply)
{
  std::mt19937_64 rng(0x5145A4E5ULL);
  for (SizeT n : {40u, 700u, 1100u}) // Karatsuba band + CRT band
  {
    for (int t = 0; t < 3; ++t)
    {
      auto a = AdversarialLimbs(n, rng);
      auto sq = Square(a, BigInteger::Base());
      auto mu = Multiply(a, a, BigInteger::Base());
      ASSERT_EQ(Compare(sq, mu), 0);
    }
  }
}

// ── decimal I/O round-trip on adversarial limbs ──────────────────────────────

REGISTER_TEST(AdversarialFuzz, ToStringParseRoundTrip)
{
  std::mt19937_64 rng(0x70537217ULL);
  for (SizeT n : {1u, 7u, 60u, 900u})
  {
    for (int t = 0; t < 3; ++t)
    {
      BigInteger v(AdversarialLimbs(n, rng), (t % 2) == 1);
      BigInteger back = BigIntegerBuilder::From(ToString(v).c_str());
      ASSERT_EQ(back.CompareTo(v), 0);
    }
  }
}
