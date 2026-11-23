/**
 * BigMath: magnitude helpers for the BigInteger operators, with a stack-buffer
 * fast path for operands short enough to be stored inline.
 *
 * Operands of at most MaxLimbs (16) limbs are combined in a local array with
 * the limb kernels and the result is built with BigInteger::FromLimbs: no
 * allocation at all when the result fits the inline storage
 * (BigInteger::InlineLimbs), otherwise exactly one, for the result. At these
 * sizes the dispatchers would run the same schoolbook multiply and Knuth
 * division; the stack path skips their vectors and bookkeeping. Everything
 * else goes through the vector-based dispatchers. The fast paths need true 64-bit
 * limbs (the kernels' carry chains); the legacy 32-in-64 layout always takes
 * the dispatcher.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_SMALL_OPS
#define BIGINTEGER_SMALL_OPS

#include <algorithm>
#include <stdexcept>
#include <utility>

#include "../BigInteger.h"
#include "../algorithms/Addition.h"
#include "../algorithms/ParallelAddSub.h"
#include "../algorithms/Subtraction.h"
#include "../common/BitShifts.h"
#include "../common/LimbKernels.h"
#include "../algorithms/division/ClassicDivision.h"
#include "../algorithms/division/FastDivision.h"

namespace BigMath
{
  namespace SmallOps
  {
#if BIGMATH_LIMB_64
    inline constexpr bool Enabled = true;
#else
    inline constexpr bool Enabled = false;
#endif
#ifndef BIGMATH_SMALL_OPS_MAX_LIMBS
#define BIGMATH_SMALL_OPS_MAX_LIMBS 16
#endif
    inline constexpr SizeT MaxLimbs = BIGMATH_SMALL_OPS_MAX_LIMBS;
    static_assert(MaxLimbs >= BigInteger::InlineLimbs, "stack path must cover inline values");
    inline constexpr DataT ZeroLimb = 0;
    inline constexpr DataT OneLimb = 1;

    inline bool Fits(BigInteger const &a)
    {
      return Enabled && a.size() <= MaxLimbs;
    }

    inline bool Fits(BigInteger const &a, BigInteger const &b)
    {
      return Enabled && a.size() <= MaxLimbs && b.size() <= MaxLimbs;
    }

    // |a| + |b|, signed `negative`.
    inline BigInteger AddMagnitudes(BigInteger const &a, BigInteger const &b, bool negative)
    {
      if (Fits(a, b))
      {
        auto x = a.Limbs();
        auto y = b.Limbs();
        if (x.size() < y.size())
          std::swap(x, y);
        const SizeT la = (SizeT)x.size(), lb = (SizeT)y.size();
        DataT t[MaxLimbs + 1];
        DataT carry = LimbKernel::AddN(t, x.data(), y.data(), lb);
        if (la > lb)
          carry = LimbKernel::AddCarry(t + lb, x.data() + lb, la - lb, carry);
        t[la] = carry;
        return BigInteger::FromLimbs(std::span<const DataT>(t, la + 1), negative);
      }
      const SizeT blocks = AddSubDetail::DispatchBlocks<false>(a.size(), b.size());
      if (blocks > 1)
        return BigInteger(AddSubDetail::ParallelResult<false>(a.Limbs(), b.Limbs(), blocks), negative);
      return BigInteger(detail::AddLimbs(a.Limbs(), b.Limbs(), BigInteger::Base()), negative);
    }

    // |a| - |b| for |a| >= |b|, signed `negative`.
    inline BigInteger SubMagnitudes(BigInteger const &a, BigInteger const &b, bool negative)
    {
      if (Fits(a, b))
      {
        auto x = a.Limbs();
        auto y = b.Limbs();
        const SizeT la = (SizeT)x.size(), lb = (SizeT)y.size();
        DataT t[MaxLimbs];
        DataT borrow = LimbKernel::SubN(t, x.data(), y.data(), lb);
        if (la > lb)
          LimbKernel::SubBorrow(t + lb, x.data() + lb, la - lb, borrow);
        return BigInteger::FromLimbs(std::span<const DataT>(t, la), negative);
      }
      const SizeT blocks = AddSubDetail::DispatchBlocks<true>(a.size(), b.size());
      if (blocks > 1)
        return BigInteger(AddSubDetail::ParallelResult<true>(a.Limbs(), b.Limbs(), blocks), negative);
      return BigInteger(detail::SubtractLimbs(a.Limbs(), b.Limbs(), BigInteger::Base()), negative);
    }

    // |a| * |b| into the stack; caller checks Fits(a, b).
    inline BigInteger MulMagnitudes(BigInteger const &a, BigInteger const &b, bool negative)
    {
      auto x = a.Limbs();
      auto y = b.Limbs();
      DataT t[2 * MaxLimbs];
      if (&a == &b) // x * x: the squaring kernel does about half the products
        LimbKernel::SqrBasecase(t, x.data(), (SizeT)x.size());
      else
        LimbKernel::MulBasecase(t, x.data(), (SizeT)x.size(), y.data(), (SizeT)y.size());
      return BigInteger::FromLimbs(std::span<const DataT>(t, x.size() + y.size()), negative);
    }

    // |a| * b into the stack; caller checks Fits(a).
    inline BigInteger MulScalar(BigInteger const &a, DataT b, bool negative)
    {
      auto x = a.Limbs();
      DataT t[MaxLimbs + 1];
      t[x.size()] = LimbKernel::Mul1(t, x.data(), (SizeT)x.size(), b);
      return BigInteger::FromLimbs(std::span<const DataT>(t, x.size() + 1), negative);
    }

    // Truncated division of small operands: q = a / b rounded toward zero
    // (sign a.neg != b.neg), r = a - q*b (sign of a). Either output may be
    // null. Caller checks Fits(a, b). Throws on b == 0.
    inline void DivRem(BigInteger const &a, BigInteger const &b, BigInteger *q, BigInteger *r)
    {
      if (b.Zero())
        throw std::invalid_argument("Division by zero");
      const bool qNeg = a.IsNegative() != b.IsNegative();
      const bool rNeg = a.IsNegative();
      auto x = a.Limbs();
      auto y = b.Limbs();

      const Int cmp = Compare(x, y);
      if (cmp <= 0)
      {
        if (q)
          *q = BigInteger::FromLimbs(std::span<const DataT>(cmp == 0 ? &OneLimb : &ZeroLimb, 1), qNeg);
        if (r)
          *r = cmp == 0 ? BigInteger() : BigInteger::FromLimbs(x, rNeg);
        return;
      }

      const SizeT na = (SizeT)x.size(), nb = (SizeT)y.size();
      DataT qt[MaxLimbs];
      if (nb == 1)
      {
        const DataT rem = ClassicDivision::DivRem1(qt, x.data(), na, y[0]);
        if (q)
          *q = BigInteger::FromLimbs(std::span<const DataT>(qt, na), qNeg);
        if (r)
          *r = BigInteger::FromLimbs(std::span<const DataT>(&rem, 1), rNeg);
        return;
      }

      // Knuth D on normalized copies: shift so the divisor's top bit is set,
      // with the bits shifted out of the dividend in one extra top limb.
      const auto s = static_cast<UInt>(NormalizationShiftBits(y[nb - 1], 64));
      DataT vn[MaxLimbs];
      DataT un[MaxLimbs + 1];
      if (s > 0)
      {
        LimbKernel::Lshift(vn, y.data(), nb, s);
        un[na] = LimbKernel::Lshift(un, x.data(), na, s);
      }
      else
      {
        std::copy(y.begin(), y.end(), vn);
        std::copy(x.begin(), x.end(), un);
        un[na] = 0;
      }
      const ULong vinv = FastDivision::Reciprocal3by2((ULong)vn[nb - 1], (ULong)vn[nb - 2]);
      FastDivision::DivideNormalized64(qt, un, na, vn, nb, vinv);
      if (q)
        *q = BigInteger::FromLimbs(std::span<const DataT>(qt, na - nb + 1), qNeg);
      if (r)
      {
        if (s > 0)
          LimbKernel::Rshift(un, un, nb, s);
        *r = BigInteger::FromLimbs(std::span<const DataT>(un, nb), rNeg);
      }
    }

  }
}

#endif
