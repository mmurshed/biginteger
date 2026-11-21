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
