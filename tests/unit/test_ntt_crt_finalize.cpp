// CRT NTT finalize (Garner reconstruction + carry) vs a serial scalar
// reference.
//
// FinalizeProduct computes the mixed-radix digits four coefficients at a time
// (NEON on aarch64) and splits the carry chain into per-thread limb ranges
// that are stitched afterwards. Both must reproduce the plain one-coefficient-
// at-a-time loop exactly, including when a carry ripples across every range
// boundary and out of the top limb.

#include "unit_test_framework.h"

#include <cstdint>
#include <random>
#include <vector>

#include "biginteger/algorithms/multiplication/NTTMultiplicationCrt.h"

using namespace BigMath;
using namespace BigMath::NttCrt;

namespace
{
  struct Planes
  {
    std::vector<UInt> f1, f2, f3;
  };

  // Residue planes for the given per-coefficient values (each < p1·p2·p3).
  Planes FromValues(const std::vector<ULong128> &x)
  {
    Planes p;
    p.f1.resize(x.size());
    p.f2.resize(x.size());
    p.f3.resize(x.size());
    for (size_t i = 0; i < x.size(); ++i)
    {
      p.f1[i] = (UInt)(x[i] % P1);
      p.f2[i] = (UInt)(x[i] % P2);
      p.f3[i] = (UInt)(x[i] % P3);
    }
    return p;
  }

  // The original serial finalize: one scalar Garner per coefficient, 32-bit
  // digits packed two per 64-bit limb.
  std::vector<DataT> ReferenceFinalize64(const Planes &p, size_t coeffCount)
  {
    const InvTable &inv = GarnerInverses();
    std::vector<DataT> result;
    ULong128 carry = 0;
    ULong acc = 0;
    int slot = 0;
    auto emit = [&](ULong digit) {
      acc |= digit << (32 * slot);
      if (++slot == 2)
      {
        result.push_back((DataT)acc);
        acc = 0;
        slot = 0;
      }
    };
    for (size_t i = 0; i < coeffCount; ++i)
    {
      ULong128 total = Garner(p.f1[i], p.f2[i], p.f3[i], inv) + carry;
      emit((ULong)(total & 0xFFFFFFFFULL));
      carry = total >> 32;
    }
    while (carry > 0)
    {
      emit((ULong)(carry & 0xFFFFFFFFULL));
      carry >>= 32;
    }
    if (slot != 0)
      result.push_back((DataT)acc);
    while (!result.empty() && result.back() == 0)
      result.pop_back();
    return result;
  }

  void Check(const Planes &p, size_t coeffCount)
  {
    std::vector<DataT> expect = ReferenceFinalize64(p, coeffCount);
    std::vector<DataT> got = FinalizeProduct(p.f1, p.f2, p.f3, (ULong)coeffCount, Base2_64, (SizeT)(coeffCount / 2 + 2));
    ASSERT_TRUE(expect == got);
  }

  const size_t kCounts[] = {1, 2, 3, 4, 5, 7, 8, 9, 63, 64, 65, 1023, 4095, 8191, 8192, 8193, 20001, 70001, 262143};
}

REGISTER_TEST(NttCrtFinalize, RandomResiduesMatchScalarReference)
{
  std::mt19937_64 gen(0xF17A11E5ULL);
  for (size_t count : kCounts)
  {
    Planes p;
    p.f1.resize(count);
    p.f2.resize(count);
    p.f3.resize(count);
    for (size_t i = 0; i < count; ++i)
    {
      p.f1[i] = (UInt)(gen() % P1);
      p.f2[i] = (UInt)(gen() % P2);
      p.f3[i] = (UInt)(gen() % P3);
    }
    Check(p, count);
  }
}

REGISTER_TEST(NttCrtFinalize, ExtremeResiduesMatchScalarReference)
{
  for (size_t count : kCounts)
  {
    // Every plane at its maximum, at zero, and mixed: the conditional
    // subtracts inside the digit computation sit on these boundaries.
    for (int pattern = 0; pattern < 4; ++pattern)
    {
      Planes p;
      p.f1.assign(count, (pattern & 1) ? P1 - 1 : 0);
      p.f2.assign(count, (pattern & 2) ? P2 - 1 : 0);
      p.f3.assign(count, (pattern == 1 || pattern == 2) ? P3 - 1 : 0);
      Check(p, count);
    }
  }
}

REGISTER_TEST(NttCrtFinalize, CarryRipplesAcrossEveryRange)
{
  // All-ones limbs everywhere, plus one extra 2^64 injected low down: the
  // carry it creates has to ripple through every later limb — across every
  // parallel range boundary — and out of the top.
  for (size_t count : kCounts)
  {
    if (count < 8)
      continue;
    for (size_t inject : {(size_t)0, (size_t)3, count / 2, count - 3})
    {
      std::vector<ULong128> x(count, (ULong128)0xFFFFFFFFULL);
      x[inject] += (ULong128)1 << 64;
      Check(FromValues(x), count);
    }
  }
}

REGISTER_TEST(NttCrtFinalize, LargestCoefficientValues)
{
  // x = p1·p2·p3 − 1 in every slot: the largest value Garner can return, so
  // the largest carries the chain ever sees.
  const ULong128 maxX = (ULong128)P1 * P2 * P3 - 1;
  for (size_t count : kCounts)
  {
    std::vector<ULong128> x(count, maxX);
    Check(FromValues(x), count);
  }
}
