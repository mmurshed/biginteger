// CRT NTT transform kernels vs a plain radix-2 reference.
//
// ForwardPtr/InversePtr fuse layers (radix-8 chain plus size-dependent tail
// kernels, NEON on aarch64). Every fusion is exact arithmetic mod P, so the
// result must equal the textbook radix-2 DIF/DIT chain bit for bit at every
// power-of-two length — in particular across all three log2(n) % 3 classes,
// which select different tail kernels.

#include "unit_test_framework.h"

#include <cstdint>
#include <random>
#include <vector>

#include "biginteger/algorithms/multiplication/NTTMultiplicationCrt.h"

using namespace BigMath;
using namespace BigMath::NttCrt;

namespace
{
  template <typename F>
  std::vector<UInt> RandomResidues(Int n, std::mt19937_64 &gen)
  {
    std::vector<UInt> v((size_t)n);
    for (auto &x : v)
      x = (UInt)(gen() % F::Prime);
    // Edge values: 0 and P-1 exercise every conditional-subtract boundary.
    for (Int i = 0; i < n && i < 8; ++i)
      v[(size_t)i] = (i & 1) ? F::Prime - 1 : 0;
    return v;
  }

  template <typename F>
  void ReferenceForward(std::vector<UInt> &a, const Plan<F> &plan)
  {
    Int n = (Int)a.size();
    for (Int len = n; len >= 2; len >>= 1)
      ForwardRadix2Layer<F>(a.data(), n, len, plan.forwardRoots.data());
  }

  template <typename F>
  void ReferenceInverse(std::vector<UInt> &a, const Plan<F> &plan)
  {
    Int n = (Int)a.size();
    for (Int len = 2; len <= n; len <<= 1)
      InverseRadix2Layer<F>(a.data(), n, len, plan.inverseRoots.data());
    for (auto &x : a)
      x = F::Mul(x, plan.invSize);
  }

  template <typename F, UInt G>
  void CheckField(uint64_t seed)
  {
    std::mt19937_64 gen(seed);
    for (Int logn = 1; logn <= 16; ++logn)
    {
      Int n = (Int)1 << logn;
      const Plan<F> &plan = GetPlan<F, G>(n);
      int rounds = logn <= 10 ? 6 : 2;
      for (int r = 0; r < rounds; ++r)
      {
        std::vector<UInt> input = RandomResidues<F>(n, gen);

        std::vector<UInt> ref = input;
        std::vector<UInt> got = input;
        ReferenceForward<F>(ref, plan);
        ForwardPtr<F>(got.data(), n, plan);
        ASSERT_TRUE(ref == got);

        std::vector<UInt> refInv = ref;
        std::vector<UInt> gotInv = ref;
        ReferenceInverse<F>(refInv, plan);
        InversePtr<F>(gotInv.data(), n, plan, /*scale=*/true);
        ASSERT_TRUE(refInv == gotInv);
        ASSERT_TRUE(gotInv == input); // round trip
      }
    }
  }

  // The split inverse (two half-length inverses + shared top layer, used to
  // give each prime two work units) must equal the one-piece inverse.
  template <typename F, UInt G>
  void CheckSplitInverse(uint64_t seed)
  {
    std::mt19937_64 gen(seed);
    for (Int logn = 2; logn <= 16; ++logn)
    {
      Int n = (Int)1 << logn;
      const Plan<F> &plan = GetPlan<F, G>(n);
      const Plan<F> &half = GetPlan<F, G>(n / 2);
      for (int r = 0; r < 3; ++r)
      {
        std::vector<UInt> input = RandomResidues<F>(n, gen);
        std::vector<UInt> whole = input;
        std::vector<UInt> split = input;
        InversePtr<F>(whole.data(), n, plan, /*scale=*/true);
        InverseSplitPtr<F>(split.data(), n, plan, half);
        ASSERT_TRUE(whole == split);
      }
    }
  }
}

REGISTER_TEST(NttCrtLayers, SplitInverseMatchesWholeInverse)
{
  CheckSplitInverse<F1, G1>(0x5EED0011ULL);
  CheckSplitInverse<F2, G2>(0x5EED0012ULL);
  CheckSplitInverse<F3, G3>(0x5EED0013ULL);
}

REGISTER_TEST(NttCrtLayers, ForwardInverseMatchRadix2ReferenceP1)
{
  CheckField<F1, G1>(0x5EED0001ULL);
}

REGISTER_TEST(NttCrtLayers, ForwardInverseMatchRadix2ReferenceP2)
{
  CheckField<F2, G2>(0x5EED0002ULL);
}

REGISTER_TEST(NttCrtLayers, ForwardInverseMatchRadix2ReferenceP3)
{
  CheckField<F3, G3>(0x5EED0003ULL);
}

#if BIGMATH_NTT_TRUNC && BIGMATH_NTT_MFA && BIGMATH_NTT_MFA_FUSE
REGISTER_TEST(NttCrtLayers, ProductAcrossTransformCeiling)
{
  std::mt19937_64 random(0x484486);
  // A small H exercises the same field twists and reconstruction as 2^26.
  // Cover both sides of the short-prefix cutoff and compare every output limb.
  constexpr Int H = 512;
  for (BaseT base : {Base2_32, Base2_64})
  {
    const SizeT perLimb = CoeffsPerLimb(base);
    for (int totalCoeffs : {H + 2, H + 64, H + 66, H + 128, H + 130, H + 160,
                            H + 162, H + 256, H + 258, 2 * H})
    {
      const SizeT totalLimbs = totalCoeffs / perLimb;
      // Alternate balanced and uneven splits, up to one full H-length input.
      const SizeT left = totalCoeffs % 4 ? totalLimbs / 2 : std::min<SizeT>(H / perLimb, 3 * totalLimbs / 5);
      std::vector<DataT> a(left), b(totalLimbs - left);
      const DataT mask = base == Base2_64 ? ~DataT{0} : 0xffffffffu;
      for (int pattern : {0, 1, 2})
      {
        for (auto *v : {&a, &b})
          for (auto &x : *v)
            x = pattern == 1 ? mask : (pattern == 2 && random() % 4 ? 0 : random() & mask);
        auto check = [&](const std::vector<DataT> &rhs) {
          const Int count = (Int)((a.size() + rhs.size()) * perLimb - 1);
          if (count <= H) return;
          ASSERT_TRUE(MultiplyAtCrtCeiling(a, rhs, base, H) ==
                      ClassicMultiplication::Multiply(a, rhs, base));
        };
        check(b);
        check(a);
      }
    }
  }
}

#endif

#if BIGMATH_NTT_MFA && BIGMATH_NTT_MFA_FUSE
REGISTER_TEST(NttCrtLayers, ReuseTransposePlanesAcrossPrimes)
{
  const Int n = 1 << 15;
  if (!MfaFusedWindow(n)) return;
  std::mt19937_64 random(0x2badd00d);
  std::vector<DataT> a(n / 4), b(n / 4);
  for (auto &x : a) x = random();
  for (auto &x : b) x = random();
  PackedOperandView va(a, Base2_64), vb(b, Base2_64);
  std::vector<UInt> sa(n), sb(n), f1(n), f2(n), f3(n);
  auto prime = [&]<typename F, UInt G>(std::vector<UInt> &result) {
    MfaPlanTree<F> tree;
    BuildMfaPlanTree<F, G>(n, tree);
    MfaFusedProductPrime<F>(va, vb, result.data(), sa.data(), sb.data(), n, tree);
  };
  prime.template operator()<F1, G1>(f1);
  prime.template operator()<F2, G2>(f2);
  prime.template operator()<F3, G3>(f3);
  ASSERT_TRUE(FinalizeProduct(f1, f2, f3, n - 1, Base2_64, n / 2) ==
              ClassicMultiplication::Multiply(a, b, Base2_64));
}
#endif

REGISTER_TEST(NttCrtLayers, BlockedProductsMatchSchoolbook)
{
  std::mt19937_64 random(0x320324);
  // A small transform budget forces the production block loops, carry
  // propagation, uneven tails, zero blocks, and the square diagonal.
  for (BaseT base : {Base2_32, Base2_64})
  {
    const DataT mask = base == Base2_64 ? ~DataT{0} : 0xffffffffu;
    for (size_t shortLen : {1, 7, 16, 17, 36, 79})
    {
      std::vector<DataT> a(79), b(shortLen);
      for (auto &v : a) v = random() & mask;
      for (auto &v : b) v = random() & mask;
      for (int pattern = 0; pattern < 3; ++pattern)
      {
        if (pattern == 1) { std::fill(a.begin(), a.end(), mask); std::fill(b.begin(), b.end(), mask); }
        if (pattern == 2) std::fill(a.begin() + 10, a.end() - 10, 0);
        const auto want = ClassicMultiplication::Multiply(a, b, base);
        ASSERT_TRUE(MultiplyBlocked(a, b, base, 32) == want);
        ASSERT_TRUE(MultiplyBlocked(b, a, base, 32) == want);
        ASSERT_TRUE(MultiplyBlocked(a, a, base, 32) == ClassicMultiplication::Multiply(a, a, base));
      }
    }
  }
}
