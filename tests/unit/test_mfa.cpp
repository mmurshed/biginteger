// MFA (Bailey 6-step) transform regression tests.
//
// The MFA dispatch gate sits at 2^24 coefficients — far above every other
// correctness harness's operand sizes — so an MFA-only bug is invisible to
// the rest of the suite. That is exactly how the PR #72 inverse-twiddle
// ordering bug (fixed in PR #103) shipped silently and corrupted every
// >= 40M-digit product for 12 days: the fused inverse cross-twiddle ran
// AFTER the inverse row FFT instead of before it, and Newton division's
// fixup fallback masked the corruption as mysteriously slow runs.
//
// These tests call ForwardMFA/InverseMFA directly (bypassing the size gate)
// so the MFA code paths — including the recursive split, the fused
// cross-twiddle, and the transpose-fusion tiles — are exercised at
// unit-test-sized transforms for all three CRT primes.

#include "unit_test_framework.h"

#include <random>
#include <vector>

#include "biginteger/algorithms/multiplication/NTTMultiplicationCrt.h"

using namespace BigMath;
using namespace BigMath::NttCrt;

namespace
{
  // Forward·Inverse must be the identity. Catches ordering/addressing bugs in
  // any stage (transpose, cross-twiddle, leaf FFTs, scaling).
  template <typename F, UInt G>
  bool RoundTripOk(Int n, uint64_t seed)
  {
    std::mt19937_64 gen(seed);
    std::uniform_int_distribution<uint32_t> dis(0, F::Prime - 1);

    std::vector<UInt> orig(n), buf(n), scratch(n);
    for (Int i = 0; i < n; ++i)
      orig[i] = dis(gen);

    buf = orig;
    // parallel=false: inside the aggregated unit binary other tests have
    // already exercised the thread pool, and MFA's parallel path assumes the
    // top-level pre-warm of per-worker plan caches that NttCrt::Multiply
    // performs. Serial keeps the test deterministic and pool-state-free.
    ForwardMFA<F, G>(buf.data(), n, scratch.data(), /*parallel=*/false);
    InverseMFA<F, G>(buf.data(), n, scratch.data(), /*parallel=*/false);
    return buf == orig;
  }

  // MFA convolution must equal the plain (non-MFA) transform pipeline's:
  // forward both inputs, multiply pointwise (both sides share the same
  // MFA permutation, so index-aligned), inverse — then compare against the
  // plan-based Forward/Inverse doing the same. Catches bugs that round-trip
  // cleanly but permute the spectrum inconsistently between the two paths.
  template <typename F, UInt G>
  bool ConvolutionMatchesPlainOk(Int n, uint64_t seed)
  {
    std::mt19937_64 gen(seed);
    // Small values so plain cyclic convolution sums stay below F::Prime — the
    // comparison is transform-level, not multiply-level, so no carrying.
    std::uniform_int_distribution<uint32_t> dis(0, 1023);

    std::vector<UInt> a(n), b(n);
    for (Int i = 0; i < n; ++i)
    {
      a[i] = dis(gen);
      b[i] = dis(gen);
    }

    // MFA pipeline.
    std::vector<UInt> fa = a, fb = b, scratch(n);
    ForwardMFA<F, G>(fa.data(), n, scratch.data(), /*parallel=*/false);
    ForwardMFA<F, G>(fb.data(), n, scratch.data(), /*parallel=*/false);
    for (Int i = 0; i < n; ++i)
      fa[i] = F::Mul(fa[i], fb[i]);
    InverseMFA<F, G>(fa.data(), n, scratch.data(), /*parallel=*/false);

    // Plain pipeline.
    const Plan<F> &plan = GetPlan<F, G>(n);
    std::vector<UInt> ga = a, gb = b;
    ForwardPtr<F>(ga.data(), n, plan);
    ForwardPtr<F>(gb.data(), n, plan);
    for (Int i = 0; i < n; ++i)
      ga[i] = F::Mul(ga[i], gb[i]);
    InversePtr<F>(ga.data(), n, plan, /*scale=*/true);

    return fa == ga;
  }
}

// n must exceed BIGMATH_NTT_MFA_LEAF (2^13) for the MFA decomposition (and
// its cross-twiddle) to actually run; 2^14 splits once, 2^15 covers the
// uneven n1 != n2 factorization, 2^17 adds a deeper recursion level.

REGISTER_TEST(MfaTransform, RoundTripP1)
{
  ASSERT_TRUE((RoundTripOk<F1, G1>(1 << 14, 0xAA1)));
  ASSERT_TRUE((RoundTripOk<F1, G1>(1 << 15, 0xAA2)));
  ASSERT_TRUE((RoundTripOk<F1, G1>(1 << 17, 0xAA3)));
}

REGISTER_TEST(MfaTransform, RoundTripP2)
{
  ASSERT_TRUE((RoundTripOk<F2, G2>(1 << 14, 0xBB1)));
  ASSERT_TRUE((RoundTripOk<F2, G2>(1 << 15, 0xBB2)));
}

REGISTER_TEST(MfaTransform, RoundTripP3)
{
  ASSERT_TRUE((RoundTripOk<F3, G3>(1 << 14, 0xCC1)));
  ASSERT_TRUE((RoundTripOk<F3, G3>(1 << 15, 0xCC2)));
}

REGISTER_TEST(MfaTransform, ConvolutionMatchesPlainP1)
{
  ASSERT_TRUE((ConvolutionMatchesPlainOk<F1, G1>(1 << 14, 0xDD1)));
  ASSERT_TRUE((ConvolutionMatchesPlainOk<F1, G1>(1 << 15, 0xDD2)));
}

REGISTER_TEST(MfaTransform, ConvolutionMatchesPlainP2)
{
  ASSERT_TRUE((ConvolutionMatchesPlainOk<F2, G2>(1 << 14, 0xEE1)));
}
