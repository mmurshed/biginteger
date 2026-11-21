// CRT NTT truncated transform length (the H + L scheme).
//
// Products whose coefficient count c lands just above a power of two H are
// computed mod (x^H + 1) and mod (x^L − 1) and recombined, instead of with a
// cyclic transform of length 2H. These tests pin the product against
// Karatsuba across the band where the scheme engages, at the L boundaries,
// for skewed operands whose own coefficients spill past H, for squaring, and
// just outside the band (where the full transform must still be used).

#include "unit_test_framework.h"

#include <cstdint>
#include <random>
#include <vector>

#include "biginteger/algorithms/multiplication/KaratsubaMultiplication.h"
#include "biginteger/algorithms/multiplication/NTTMultiplicationCrt.h"

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

#endif
