// CRT NTT truncated transform length (the H + L scheme).
//
// Products whose coefficient count c lands just above a power of two H are
// computed mod (x^H + 1) and mod (x^L − 1) and recombined, instead of with a
// cyclic transform of length 2H. These tests pin the product against
// Karatsuba across the band where the scheme engages, at the L boundaries,
// for skewed operands whose own coefficients spill past H, for squaring, and
// just outside the band (where the full transform must still be used).

#include "unit_test_framework.h"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

#include "biginteger/algorithms/multiplication/KaratsubaMultiplication.h"
#include "biginteger/algorithms/multiplication/NTTMultiplicationCrt.h"
#include "biginteger/common/Parallel.h"

using namespace BigMath;

namespace
{
  std::vector<DataT> Limbs(SizeT n, int kind, std::mt19937_64 &gen)
  {
    std::vector<DataT> v(n);
    for (auto &x : v)
    {
      switch (kind)
      {
        case 0: x = gen(); break;
        case 1: x = ~(DataT)0; break;                       // every coefficient 2^32 − 1
        default: x = (gen() % 7 == 0) ? gen() : (DataT)0; break; // sparse
      }
    }
    if (v.back() == 0)
      v.back() = 1;
    return v;
  }

  ULong CoeffCount(SizeT la, SizeT lb) { return (ULong)2 * la + (ULong)2 * lb - 1; }

  void CheckProduct(SizeT la, SizeT lb, int kind, std::mt19937_64 &gen)
  {
    auto a = Limbs(la, kind, gen);
    auto b = Limbs(lb, kind, gen);
    auto expect = KaratsubaMultiplication::Multiply(a, b, Base2_64);
    ASSERT_TRUE(NttCrt::Multiply(a, b, Base2_64) == expect);
    ASSERT_TRUE(NttCrt::Multiply(b, a, Base2_64) == expect);
  }
}

#if BIGMATH_NTT_TRUNC && BIGMATH_LIMB_64

REGISTER_TEST(NttCrtTruncated, BalancedAcrossTheBand)
{
  std::mt19937_64 gen(0x7256AC7EDULL);
  for (Int H : {2048, 4096, 8192, 16384})
  {
    if (H < (Int)BIGMATH_NTT_TRUNC_MIN || H > (Int)BIGMATH_NTT_TRUNC_MAX)
      continue;
    // Total limbs la + lb with 2(la + lb) − 1 coefficients: from the first
    // size above H through the top of the band.
    const SizeT first = (SizeT)(H / 2 + 1);
    const SizeT last = (SizeT)((H + (H >> BIGMATH_NTT_TRUNC_BAND_SHIFT) + 1) / 2);
    for (SizeT total : {first, first + 1, first + 31, first + 32, first + 33, first + 64, first + 65,
                        (SizeT)(H / 2 + H / 32), (SizeT)(H / 2 + H / 16), (SizeT)(H / 2 + H / 16 + 1),
                        (SizeT)(H / 2 + H / 8), (SizeT)(H / 2 + H / 8 + 1), last - 1, last})
    {
      SizeT la = total / 2, lb = total - la;
      ASSERT_TRUE(NttCrt::TruncCyclicLength(CoeffCount(la, lb), H) != 0);
      for (int kind = 0; kind < 3; ++kind)
        CheckProduct(la, lb, kind, gen);
    }
  }
}

REGISTER_TEST(NttCrtTruncated, SkewedOperandSpillsPastH)
{
  // The long operand alone has more than H coefficients, so its upper part
  // wraps (negated) into the weighted half and folds into the cyclic half.
  std::mt19937_64 gen(0x5B111ULL);
  for (Int H : {2048, 4096, 8192})
  {
    if (H < (Int)BIGMATH_NTT_TRUNC_MIN || H > (Int)BIGMATH_NTT_TRUNC_MAX)
      continue;
    for (SizeT lb : {2, 7, 30, 100})
    {
      for (SizeT extra : {1, 9, 40})
      {
        SizeT la = (SizeT)(H / 2) + extra; // 2·la > H coefficients
        if (NttCrt::TruncCyclicLength(CoeffCount(la, lb), H) == 0)
          continue;
        for (int kind = 0; kind < 3; ++kind)
          CheckProduct(la, lb, kind, gen);
      }
    }
    // And make sure at least one such shape really is in the band.
    ASSERT_TRUE(NttCrt::TruncCyclicLength(CoeffCount((SizeT)(H / 2) + 9, 30), H) != 0);
  }
}

REGISTER_TEST(NttCrtTruncated, Squaring)
{
  std::mt19937_64 gen(0x50AAEULL);
  for (Int H : {2048, 4096, 8192})
  {
    if (H < (Int)BIGMATH_NTT_TRUNC_MIN || H > (Int)BIGMATH_NTT_TRUNC_MAX)
      continue;
    for (SizeT la : {(SizeT)(H / 4 + 1), (SizeT)(H / 4 + 20), (SizeT)(H / 4 + H / 32)})
    {
      ASSERT_TRUE(NttCrt::TruncCyclicLength(CoeffCount(la, la), H) != 0);
      for (int kind = 0; kind < 3; ++kind)
      {
        auto a = Limbs(la, kind, gen);
        auto expect = KaratsubaMultiplication::Multiply(a, a, Base2_64);
        ASSERT_TRUE(NttCrt::Multiply(a, a, Base2_64) == expect); // &a == &b: self-multiply path
      }
    }
  }
}

REGISTER_TEST(NttCrtTruncated, BandEdges)
{
  // c == H exactly and c just past the top of the band use the full
  // transform; both sides of each edge must agree with Karatsuba.
  std::mt19937_64 gen(0xED6E5ULL);
  const Int H = 4096;
  const Int band = H >> BIGMATH_NTT_TRUNC_BAND_SHIFT;
  ASSERT_EQ(NttCrt::TruncCyclicLength((ULong)H, H), 0);
  ASSERT_EQ(NttCrt::TruncCyclicLength((ULong)(H + band + 1), H), 0);
  ASSERT_TRUE(NttCrt::TruncCyclicLength((ULong)(H + 1), H) == 64 || H < (Int)BIGMATH_NTT_TRUNC_MIN || H > (Int)BIGMATH_NTT_TRUNC_MAX);
  for (SizeT total : {(SizeT)(H / 2), (SizeT)(H / 2 + 1), (SizeT)((H + band) / 2), (SizeT)((H + band) / 2 + 1), (SizeT)((H + band) / 2 + 2)})
  {
    SizeT la = total / 2, lb = total - la;
    for (int kind = 0; kind < 3; ++kind)
      CheckProduct(la, lb, kind, gen);
  }
}

// ─── Single-thread cache-blocked variant (CrtTruncSerial) ────────────────────
//
// NttCrt::Multiply takes this path only above BIGMATH_NTT_TRUNC_MAX and only
// when the process has one thread, so these tests call MultiplyTruncSerial
// directly: the same pipeline, in any build, with explicit H and L and an MFA
// cutover low enough to exercise the fused stages at unit-test sizes.

#if BIGMATH_NTT_MFA

namespace
{
  // v mod m, limbs little-endian.
  uint64_t Residue(const std::vector<DataT> &v, uint64_t m)
  {
    unsigned __int128 r = 0;
    for (SizeT i = v.size(); i-- > 0;)
      r = ((r << 64) | v[i]) % m;
    return (uint64_t)r;
  }

  // a · b ≡ c modulo two 61/62-bit moduli. Linear time, so it can vouch for
  // products far too large to recompute with Karatsuba.
  bool ResiduesAgree(const std::vector<DataT> &a, const std::vector<DataT> &b, const std::vector<DataT> &c)
  {
    for (uint64_t m : {(uint64_t)0x1FFFFFFFFFFFFFFFULL, (uint64_t)0x3FFFFFFFFFFFFFC7ULL}) // 2^61 − 1, 2^62 − 57
    {
      const unsigned __int128 lhs = (unsigned __int128)Residue(a, m) * Residue(b, m) % m;
      if ((uint64_t)lhs != Residue(c, m))
        return false;
    }
    return true;
  }

  void CheckSerial(SizeT la, SizeT lb, Int H, Int L, Int mfaMin, int kind, std::mt19937_64 &gen)
  {
    auto a = Limbs(la, kind, gen);
    auto b = Limbs(lb, kind, gen);
    auto expect = KaratsubaMultiplication::Multiply(a, b, Base2_64);
    ASSERT_TRUE(NttCrt::MultiplyTruncSerial(a, b, Base2_64, H, L, mfaMin) == expect);
    ASSERT_TRUE(NttCrt::MultiplyTruncSerial(b, a, Base2_64, H, L, mfaMin) == expect);
  }
}

REGISTER_TEST(NttCrtTruncatedSerial, WholeArrayChain)
{
  // mfaMin above every length here: both transforms on the radix-8 chain.
  std::mt19937_64 gen(0x5E21A1ULL);
  const Int mfaMin = 1 << 30;
  for (Int H : {2048, 8192})
  {
    for (Int L : {(Int)64, H / 8, H / 2, H})
    {
      // c = 2·total − 1 coefficients: the first size above H, and the last
      // that fits H + L.
      for (SizeT total : {(SizeT)(H / 2 + 1), (SizeT)((H + L) / 2)})
      {
        SizeT la = total / 2, lb = total - la;
        for (int kind = 0; kind < 3; ++kind)
          CheckSerial(la, lb, H, L, mfaMin, kind, gen);
      }
    }
  }
}

REGISTER_TEST(NttCrtTruncatedSerial, FusedMfaStages)
{
  // mfaMin = 2^14: the length-H transform always runs the fused MFA stages
  // (2^15 = 2^8 · 2^7, 2^16 = 2^8 · 2^8); L covers the chain (64, and 2^13 =
  // the MFA leaf), MFA below H, and L = H.
  std::mt19937_64 gen(0x5E21A2ULL);
  const Int mfaMin = 1 << 14;
  for (Int H : {1 << 15, 1 << 16})
  {
    for (Int L : {(Int)64, (Int)(1 << 13), (Int)(1 << 14), H / 2, H})
    {
      for (SizeT total : {(SizeT)(H / 2 + 1), (SizeT)((H + L) / 2)})
      {
        SizeT la = total / 2, lb = total - la;
        CheckSerial(la, lb, H, L, mfaMin, 0, gen);
      }
    }
    // Carry-heavy and sparse inputs at the widest shape.
    for (int kind = 1; kind < 3; ++kind)
      CheckSerial((SizeT)(H / 2), (SizeT)(H / 2), H, H, mfaMin, kind, gen);
  }
}

REGISTER_TEST(NttCrtTruncatedSerial, SkewedOperandSpillsPastH)
{
  // The long operand alone has more than H coefficients.
  std::mt19937_64 gen(0x5E21A3ULL);
  const Int H = 1 << 15, mfaMin = 1 << 14;
  for (SizeT lb : {2, 30, 1000})
  {
    for (SizeT extra : {1, 40, 3000})
    {
      const SizeT la = (SizeT)(H / 2) + extra; // 2·la > H coefficients
      const ULong c = CoeffCount(la, lb);
      const Int L = std::max<Int>((Int)std::bit_ceil(c - (ULong)H), 64);
      for (int kind = 0; kind < 3; ++kind)
        CheckSerial(la, lb, H, L, mfaMin, kind, gen);
    }
  }
}

REGISTER_TEST(NttCrtTruncatedSerial, Squaring)
{
  // The same vector twice takes the one-operand flow (one forward transform,
  // spectrum squared in place).
  std::mt19937_64 gen(0x5E21A4ULL);
  for (Int mfaMin : {1 << 14, 1 << 30})
  {
    const Int H = 1 << 15;
    for (Int L : {(Int)64, H / 4, H})
    {
      const SizeT la = (SizeT)((H + L) / 4);
      for (int kind = 0; kind < 3; ++kind)
      {
        auto a = Limbs(la, kind, gen);
        auto expect = KaratsubaMultiplication::Multiply(a, a, Base2_64);
        ASSERT_TRUE(NttCrt::MultiplyTruncSerial(a, a, Base2_64, H, L, mfaMin) == expect);
      }
    }
  }
}

REGISTER_TEST(NttCrtTruncatedSerial, RejectsLengthsThatDoNotFit)
{
  std::mt19937_64 gen(0x5E21A5ULL);
  auto a = Limbs(600, 0, gen); // 2399 coefficients with itself
  auto expectThrow = [&](Int H, Int L) {
    bool threw = false;
    try { (void)NttCrt::MultiplyTruncSerial(a, a, Base2_64, H, L); }
    catch (const std::invalid_argument &) { threw = true; }
    return threw;
  };
  ASSERT_TRUE(expectThrow(4096, 64));   // product fits H: nothing to truncate
  ASSERT_TRUE(expectThrow(2048, 256));  // 2399 > H + L
  ASSERT_TRUE(expectThrow(2048, 4096)); // L > H
  ASSERT_TRUE(expectThrow(2048, 384));  // L not a power of two
  ASSERT_FALSE(expectThrow(2048, 512));
}

REGISTER_TEST(NttCrtTruncatedSerial, DispatchGate)
{
  // At or below BIGMATH_NTT_TRUNC_MAX the whole-array scheme owns the band.
  const Int low = (Int)BIGMATH_NTT_TRUNC_MAX;
  ASSERT_EQ(NttCrt::TruncSerialCyclicLength((ULong)low + 1, low), 0);
  // Above it: single-threaded only.
  const Int H = low * 2;
  if (H <= (Int)BIGMATH_NTT_TRUNC_SERIAL_MAX)
  {
    const bool serial = ParallelNumThreads() <= 1;
    ASSERT_EQ(NttCrt::TruncSerialCyclicLength((ULong)H + 1, H), serial ? 64 : 0);
    ASSERT_EQ(NttCrt::TruncSerialCyclicLength((ULong)H + (ULong)(H >> 2), H), serial ? (H >> 2) : 0);
    ASSERT_EQ(NttCrt::TruncSerialCyclicLength((ULong)H, H), 0); // fits H: full transform of length H
    // The ordinary scheme stays out of this range in every build.
    ASSERT_EQ(NttCrt::TruncCyclicLength((ULong)H + 1, H), 0);
  }
  ASSERT_EQ(NttCrt::TruncSerialCyclicLength(((ULong)BIGMATH_NTT_TRUNC_SERIAL_MAX << 1) + 1,
                                            (Int)BIGMATH_NTT_TRUNC_SERIAL_MAX << 1), 0);
}

REGISTER_TEST(NttCrtTruncatedSerial, DispatchSizes)
{
  // The lengths the dispatcher really uses: H = 2^21 with the production MFA
  // cutover. Too large for Karatsuba, so check residues, and compare with
  // NttCrt::Multiply — in a threaded build that is the full fused-MFA
  // transform, an independent pipeline.
  std::mt19937_64 gen(0x5E21A6ULL);
  const Int H = 1 << 21;
  if (H > (Int)BIGMATH_NTT_TRUNC_MAX && H <= (Int)BIGMATH_NTT_TRUNC_SERIAL_MAX)
  {
    // L = 2^19 stays on the chain (below the cutover); L = H is MFA.
    for (Int L : {1 << 19, H})
    {
      const SizeT la = (SizeT)((H + L) / 4), lb = (SizeT)((H + L) / 4) - 3;
      auto a = Limbs(la, 0, gen);
      auto b = Limbs(lb, 0, gen);
      auto c = NttCrt::MultiplyTruncSerial(a, b, Base2_64, H, L);
      ASSERT_TRUE(ResiduesAgree(a, b, c));
      ASSERT_TRUE(c == NttCrt::Multiply(a, b, Base2_64));
    }
    auto a = Limbs((SizeT)(H / 4 + H / 16), 0, gen);
    auto sq = NttCrt::MultiplyTruncSerial(a, a, Base2_64, H, H / 4);
    ASSERT_TRUE(ResiduesAgree(a, a, sq));
    ASSERT_TRUE(sq == NttCrt::Multiply(a, a, Base2_64));
  }
}

#endif // BIGMATH_NTT_MFA

#endif

#if BIGMATH_NTT_TRUNC && BIGMATH_NTT_MFA && BIGMATH_NTT_MFA_FUSE && BIGMATH_USE_THREADS
REGISTER_TEST(NttCrtTruncatedThreaded, PackedProducts)
{
  const Int H = 1 << 16;
  if (!NttCrt::MfaFusedWindow(H))
    return;
  std::mt19937_64 gen(0x7A1EAD);
  NttCrt::ProductWorkspace work;
  for (BaseT base : {Base2_64, Base2_32})
  {
    const SizeT cpl = NttCrt::CoeffsPerLimb(base);
    for (Int L : {64, 512, 1024, H / 4})
    {
      const SizeT total = (H + L) / cpl;
      auto a = Limbs(total / 2, L == 512 ? 1 : 0, gen);
      auto b = Limbs(total - a.size(), L == H / 4 ? 2 : 0, gen);
      if (base == Base2_32)
      {
        for (auto &x : a) x &= 0xffffffffu;
        for (auto &x : b) x &= 0xffffffffu;
      }
      NttCrt::CrtTruncThreaded(NttCrt::PackedOperandView(a, base), NttCrt::PackedOperandView(b, base),
                               H, L, work);
      const auto actual = NttCrt::FinalizeProduct(work.fa1, work.fa2, work.fa3, cpl * total - 1,
                                                  base, total + 2);
      ASSERT_TRUE(actual == KaratsubaMultiplication::Multiply(a, b, base));
    }
  }
}

REGISTER_TEST(NttCrtTruncatedThreaded, DispatchBand)
{
  const Int H = 1 << 21;
  if (!NttCrt::MfaFusedWindow(H) || H <= BIGMATH_NTT_TRUNC_MAX || ParallelAvailableThreads(true) <= 1)
    return;
  const SizeT a = H / 4;
  ASSERT_EQ(NttCrt::TruncThreadedCyclicLength(H - 1, H, a, a), 0);
  ASSERT_EQ(NttCrt::TruncThreadedCyclicLength(H + 1, H, a, a + 1), 64);
  ASSERT_EQ(NttCrt::TruncThreadedCyclicLength(H + 65, H, a + 16, a + 17), 128);
  ASSERT_EQ(NttCrt::TruncThreadedCyclicLength(H + H / 4 + 1, H, a + H / 16, a + H / 16 + 1), 0);
  ASSERT_EQ(NttCrt::TruncThreadedCyclicLength(H + 1, H, a + 2, a - 1), 0);

  Int nested[2] = {-1, -1};
  ParallelDo(2, [&](Int first, Int end) {
    for (Int i = first; i < end; ++i)
      nested[i] = NttCrt::TruncThreadedCyclicLength(H + 1, H, a, a + 1);
  });
  ASSERT_EQ(nested[0], 0);
  ASSERT_EQ(nested[1], 0);
}
#endif

#if BIGMATH_NTT_MFA
REGISTER_TEST(NttCrtTruncatedSerial, FactoredWeightsCrossBlockBoundaries)
{
  std::mt19937_64 gen(0xC01DF00D);
  // An odd coefficient count exercises the scalar tail after complete blocks.
  for (BaseT base : {Base2_32, Base2_64})
  {
    const Int h = 2048;
    const SizeT count = h / (2 * NttCrt::CoeffsPerLimb(base)) + 1;
    auto a = Limbs(count, 0, gen);
    auto b = Limbs(count, 0, gen);
    if (base == Base2_32)
    {
      for (auto &x : a) x &= 0xffffffffu;
      for (auto &x : b) x &= 0xffffffffu;
    }
    ASSERT_TRUE(NttCrt::MultiplyTruncSerial(a, b, base, h, 64) ==
                ClassicMultiplication::Multiply(a, b, base));
  }
}
#endif
