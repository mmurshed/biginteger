// Limb kernels (common/LimbKernels.h) vs one-limb-at-a-time references.
//
// The kernels process four limbs per straight-line block with flag-chained
// carries; the references below are the plain 128-bit-accumulator loops they
// replaced. Sizes sweep every block remainder (n % 4) and the patterns hit
// the long-carry cases (all-ones limbs, maximal multiplier).

#include "unit_test_framework.h"

#include <cstdint>
#include <random>
#include <vector>

#include "biginteger/common/LimbKernels.h"

using namespace BigMath;

namespace
{
  std::vector<DataT> Pattern(SizeT n, int kind, std::mt19937_64 &gen)
  {
    std::vector<DataT> v(n);
    for (auto &x : v)
    {
      switch (kind)
      {
        case 0: x = gen(); break;
        case 1: x = ~(DataT)0; break;
        case 2: x = 0; break;
        default: x = (gen() & 1) ? ~(DataT)0 : (DataT)(gen() & 3); break; // runs of all-ones
      }
    }
    return v;
  }

  DataT RefAdd(std::vector<DataT> &r, const std::vector<DataT> &a, const std::vector<DataT> &b, DataT carry)
  {
    for (size_t i = 0; i < a.size(); ++i)
    {
      ULong128 s = (ULong128)a[i] + b[i] + carry;
      r[i] = (DataT)s;
      carry = (DataT)(s >> 64);
    }
    return carry;
  }

  DataT RefSub(std::vector<DataT> &r, const std::vector<DataT> &a, const std::vector<DataT> &b, DataT borrow)
  {
    for (size_t i = 0; i < a.size(); ++i)
    {
      ULong128 d = (ULong128)a[i] - b[i] - borrow;
      r[i] = (DataT)d;
      borrow = (DataT)((d >> 64) & 1);
    }
    return borrow;
  }

  // r += a * b, one limb at a time; returns the carry limb.
  DataT RefAddMul1(std::vector<DataT> &r, const std::vector<DataT> &a, DataT b)
  {
    DataT carry = 0;
    for (size_t j = 0; j < a.size(); ++j)
    {
      ULong128 p = (ULong128)a[j] * b + r[j] + carry;
      r[j] = (DataT)p;
      carry = (DataT)(p >> 64);
    }
    return carry;
  }

  std::vector<DataT> RefMul(const std::vector<DataT> &a, const std::vector<DataT> &b)
  {
    std::vector<DataT> r(a.size() + b.size(), 0);
    for (size_t i = 0; i < b.size(); ++i)
    {
      DataT carry = 0;
      for (size_t j = 0; j < a.size(); ++j)
      {
        ULong128 p = (ULong128)a[j] * b[i] + r[i + j] + carry;
        r[i + j] = (DataT)p;
        carry = (DataT)(p >> 64);
      }
      r[i + a.size()] = carry;
    }
    return r;
  }
}

REGISTER_TEST(LimbKernels, AddNSubNMatchReference)
{
  std::mt19937_64 gen(0xADD5B8ULL);
  for (SizeT n = 0; n <= 41; ++n)
    for (int ka = 0; ka < 4; ++ka)
      for (int kb = 0; kb < 4; ++kb)
        for (DataT cin = 0; cin <= 1; ++cin)
        {
          auto a = Pattern(n, ka, gen), b = Pattern(n, kb, gen);
          std::vector<DataT> ref(n), got(n);

          DataT rc = RefAdd(ref, a, b, cin);
          DataT gc = LimbKernel::AddN(got.data(), a.data(), b.data(), n, cin);
          ASSERT_TRUE(ref == got);
          ASSERT_EQ(rc, gc);

          DataT rb = RefSub(ref, a, b, cin);
          DataT gb = LimbKernel::SubN(got.data(), a.data(), b.data(), n, cin);
          ASSERT_TRUE(ref == got);
          ASSERT_EQ(rb, gb);

          // Exact aliasing: r == a.
          RefAdd(ref, a, b, cin);
          got = a;
          LimbKernel::AddN(got.data(), got.data(), b.data(), n, cin);
          ASSERT_TRUE(ref == got);
          RefSub(ref, a, b, cin);
          got = a;
          LimbKernel::SubN(got.data(), got.data(), b.data(), n, cin);
          ASSERT_TRUE(ref == got);
        }
}

REGISTER_TEST(LimbKernels, ExtendedAndInPlaceForms)
{
  std::mt19937_64 gen(0x1A9E57ULL);
  for (SizeT la = 0; la <= 21; ++la)
    for (SizeT lb = 0; lb <= 21; ++lb)
      for (int kind = 0; kind < 4; ++kind)
      {
        auto a = Pattern(la, kind, gen), b = Pattern(lb, (kind + 1) & 3, gen);
        SizeT mx = std::max(la, lb);

        // AddExt: r has two limbs of headroom; the first takes the carry.
        std::vector<DataT> pa = a, pb = b;
        pa.resize(mx + 2, 0);
        pb.resize(mx + 2, 0);
        std::vector<DataT> ref(mx + 2), got(mx + 2, 0xDEADBEEFULL);
        RefAdd(ref, pa, pb, 0);
        LimbKernel::AddExt(got.data(), mx + 2, a.data(), la, b.data(), lb);
        ASSERT_TRUE(ref == got);

        // AddInPlace / SubInPlace: dest is the longer operand plus headroom.
        if (la >= lb)
        {
          std::vector<DataT> dest = pa;
          DataT rc = RefAdd(ref, pa, pb, 0);
          DataT gc = LimbKernel::AddInPlace(dest.data(), mx + 2, b.data(), lb);
          ASSERT_TRUE(ref == dest);
          ASSERT_EQ(rc, gc);

          dest = pa;
          DataT rb = RefSub(ref, pa, pb, 0);
          DataT gb = LimbKernel::SubInPlace(dest.data(), mx + 2, b.data(), lb);
          ASSERT_TRUE(ref == dest);
          ASSERT_EQ(rb, gb);
        }
      }
}

REGISTER_TEST(LimbKernels, Mul1AddMul1MatchReference)
{
  std::mt19937_64 gen(0x4D554C31ULL);
  const DataT multipliers[] = {0, 1, 2, ~(DataT)0, (DataT)1 << 63, 0x9E3779B97F4A7C15ULL};
  for (SizeT n = 0; n <= 41; ++n)
    for (int ka = 0; ka < 4; ++ka)
      for (int kr = 0; kr < 4; ++kr)
        for (DataT m : multipliers)
        {
          auto a = Pattern(n, ka, gen), r0 = Pattern(n, kr, gen);

          std::vector<DataT> ref = r0, got = r0;
          DataT rc = RefAddMul1(ref, a, m);
          DataT gc = LimbKernel::AddMul1(got.data(), a.data(), n, m);
          ASSERT_TRUE(ref == got);
          ASSERT_EQ(rc, gc);

          // mul_1 == addmul_1 into zeros; also in place.
          std::vector<DataT> zero(n, 0);
          ref = zero;
          rc = RefAddMul1(ref, a, m);
          got.assign(n, 0x5555555555555555ULL);
          gc = LimbKernel::Mul1(got.data(), a.data(), n, m);
          ASSERT_TRUE(ref == got);
          ASSERT_EQ(rc, gc);
          got = a;
          gc = LimbKernel::Mul1(got.data(), got.data(), n, m);
          ASSERT_TRUE(ref == got);
          ASSERT_EQ(rc, gc);
        }
}

REGISTER_TEST(LimbKernels, BasecaseMultiplyAndSquareMatchReference)
{
  std::mt19937_64 gen(0xBA5ECA5EULL);
  for (SizeT la = 1; la <= 37; ++la)
    for (SizeT lb = 1; lb <= 37; lb += (la % 3) + 1)
      for (int kind = 0; kind < 4; ++kind)
      {
        auto a = Pattern(la, kind, gen), b = Pattern(lb, (kind * 3 + 1) & 3, gen);
        std::vector<DataT> ref = RefMul(a, b);
        std::vector<DataT> got(la + lb, 0xABABABABABABABABULL); // r need not be initialized
        LimbKernel::MulBasecase(got.data(), a.data(), la, b.data(), lb);
        ASSERT_TRUE(ref == got);
      }

  for (SizeT n = 1; n <= 70; ++n)
    for (int kind = 0; kind < 4; ++kind)
    {
      auto a = Pattern(n, kind, gen);
      std::vector<DataT> ref = RefMul(a, a);
      std::vector<DataT> got(2 * n, 0xABABABABABABABABULL);
      LimbKernel::SqrBasecase(got.data(), a.data(), n);
      ASSERT_TRUE(ref == got);
    }
}
