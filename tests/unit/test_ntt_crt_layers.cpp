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

#include "biginteger/algorithms/multiplication/KaratsubaMultiplication.h"
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

namespace
{
  // An operand of the given coefficient count, rounded up to whole limbs.
  // Pattern 1 is all ones, the largest convolution coefficients; pattern 2 is
  // mostly zero limbs.
  std::vector<DataT> PatternOperand(size_t coeffs, BaseT base, int pattern, std::mt19937_64 &random)
  {
    const size_t perLimb = CoeffsPerLimb(base);
    const DataT mask = base == Base2_64 ? ~DataT{0} : 0xffffffffu;
    std::vector<DataT> v(std::max<size_t>(1, (coeffs + perLimb - 1) / perLimb));
    for (auto &x : v)
      x = pattern == 1 ? mask : (pattern == 2 && random() % 4 ? 0 : random() & mask);
    if (v.back() == 0) v.back() = 1;
    return v;
  }
}

REGISTER_TEST(NttCrtLayers, ResidueProductsMatchSchoolbook)
{
  std::mt19937_64 random(0x750750);
  // Every residue count, with unsigned and biased coefficients. The sums of
  // the two coefficient counts put the top block on both sides of each
  // cyclic length and of the short-prefix cutoff; the uneven split makes one
  // operand wrap around H several times.
  constexpr Int H = 512;
  for (BaseT base : {Base2_32, Base2_64})
  {
    const size_t perLimb = CoeffsPerLimb(base);
    for (Int k = 1; k <= MaxCeilingResidues; ++k)
    {
      for (size_t over : {2, 64, 65, 66, 129, 130, 162, 258, 290, 322, 400, 512, 513})
      {
        const size_t sum = (size_t)k * H + over;
        if (sum % perLimb) continue;
        for (size_t left : {sum / 2, sum / 5})
        {
          left = std::max(perLimb, left - left % perLimb);
          for (int pattern : {0, 1, 2})
          {
            const auto a = PatternOperand(left, base, pattern, random);
            const auto b = PatternOperand(sum - left, base, pattern, random);
            const auto want = ClassicMultiplication::Multiply(a, b, base);
            for (bool biased : {false, true})
            {
              ASSERT_TRUE(MultiplyAtCrtCeiling(a, b, base, H, biased) == want);
              ASSERT_TRUE(MultiplyAtCrtCeiling(b, a, base, H, biased) == want);
            }
          }
        }
        // Squares pass the same vector twice.
        const auto a = PatternOperand(sum / 2 - (sum / 2) % perLimb, base, k % 3, random);
        const size_t count = 2 * a.size() * perLimb - 1;
        if (count <= (size_t)H) continue;
        const auto want = ClassicMultiplication::Multiply(a, a, base);
        ASSERT_TRUE(MultiplyAtCrtCeiling(a, a, base, H, false) == want);
        ASSERT_TRUE(MultiplyAtCrtCeiling(a, a, base, H, true) == want);
      }
    }
  }
}

REGISTER_TEST(NttCrtLayers, ResidueProductsUseFusedTransforms)
{
  // Residues long enough for the tiled matrix transform and the pool. The
  // shapes cover a fused cyclic tail shorter than H fed by folded operands,
  // a tail of length H read through the folding view, the short prefix
  // product, and four residues with one operand longer than 2H.
  constexpr Int H = 1 << 15;
  if (!MfaFusedWindow(H) || !MfaFusedWindow(H / 2)) return;
  std::mt19937_64 random(0x1b1b);
  const std::pair<size_t, size_t> shapes[] = {
      {38000, 39000}, {37364, 37364}, {60000, 62000}, {40000, 100000}, {17000, 49000}};
  for (const auto &[aCoeffs, bCoeffs] : shapes)
  {
    for (int pattern : {0, 1})
    {
      const auto a = PatternOperand(aCoeffs, Base2_64, pattern, random);
      const auto b = PatternOperand(bCoeffs, Base2_64, pattern, random);
      const auto want = KaratsubaMultiplication::Multiply(a, b, Base2_64);
      ASSERT_TRUE(MultiplyAtCrtCeiling(a, b, Base2_64, H, false) == want);
      ASSERT_TRUE(MultiplyAtCrtCeiling(a, b, Base2_64, H, true) == want);
    }
  }
}

REGISTER_TEST(NttCrtLayers, ProductsAboveTransformChooseResiduesOrBlocks)
{
  std::mt19937_64 random(0xab0fe);
  // Small limits walk the production choices: one set of residues with
  // unsigned or biased coefficients, blocks whose pairs are residue products
  // again, and blocks within one transform for a short second operand.
  constexpr Int H = 512;
  constexpr ULong plainTerms = 400, biasedTerms = 700;
  const std::pair<size_t, size_t> shapes[] = {
      {300, 300},   // unsigned residues
      {390, 2000},  // unsigned, long partner
      {600, 640},   // biased residues
      {700, 2800},  // biased at both limits
      {702, 800},   // shorter operand past the biased limit: two blocks of it
      {2000, 2000}, // three blocks against the whole partner; 3 x 3 as a square
      {300, 5000},  // too many residues: blocks of the longer operand
      {800, 6000},  // blocks of both
      {100, 3000},  // short operand: pairs within one transform
      {256, 4000},
      {258, 258},
  };
  for (BaseT base : {Base2_32, Base2_64})
  {
    for (const auto &[aCoeffs, bCoeffs] : shapes)
    {
      for (int pattern : {0, 1, 2})
      {
        const auto a = PatternOperand(aCoeffs, base, pattern, random);
        const auto b = PatternOperand(bCoeffs, base, pattern, random);
        const auto want = ClassicMultiplication::Multiply(a, b, base);
        ASSERT_TRUE(MultiplyAboveTransform(a, b, base, H, plainTerms, biasedTerms) == want);
        ASSERT_TRUE(MultiplyAboveTransform(b, a, base, H, plainTerms, biasedTerms) == want);
        ASSERT_TRUE(MultiplyAboveTransform(b, b, base, H, plainTerms, biasedTerms) ==
                    ClassicMultiplication::Multiply(b, b, base));
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
