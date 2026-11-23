/**
 * BigInteger Class
 * Schoolbook squaring: half the partial products of Multiply(a,a).
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef CLASSIC_SQUARE
#define CLASSIC_SQUARE

#include <vector>
#include <span>
#include <cstring>

#include "../../common/LimbKernels.h"
#include "../../common/Util.h"

namespace BigMath
{
  class ClassicSquare
  {
  public:
    // Pointer variant: r[0..2n-1] = a[0..n-1]^2. r is zeroed here; callers
    // may hand in dirty buffers.
    static void SquarePtr(const DataT *a, SizeT n, DataT *r, BaseT base)
    {
      if (base == Base2_64 && n != 0)
      {
        // Off-diagonal rows on mul_1/addmul_1, then one doubling + diagonal
        // pass; writes every limb of r. See common/LimbKernels.h.
        LimbKernel::SqrBasecase(r, a, n);
        return;
      }

      std::memset(r, 0, 2 * n * sizeof(DataT));

      if (base == Base2_32)
      {
        // 1. Upper triangle: r[i+j] += a[i]*a[j]  for j > i.
        for (SizeT i = 0; i < n; ++i)
        {
          if (a[i] == 0) continue;
          ULong ai = a[i];
          ULong carry = 0;
          for (SizeT j = i + 1; j < n; ++j)
          {
            ULong prod = ai * (ULong)a[j] + (ULong)r[i + j] + carry;
            r[i + j] = (DataT)(prod & 0xFFFFFFFFULL);
            carry = prod >> 32;
          }
          r[i + n] = (DataT)carry;
        }

        // 2. Multiply by 2 (left shift by 1 bit).
        ULong c = 0;
        for (SizeT i = 0; i < 2 * n; ++i)
        {
          ULong v = ((ULong)r[i] << 1) | c;
          r[i] = (DataT)(v & 0xFFFFFFFFULL);
          c = v >> 32;
        }

        // 3. Add diagonal a[i]^2 at position 2i.
        ULong carry = 0;
        for (SizeT i = 0; i < n; ++i)
        {
          ULong sq = (ULong)a[i] * (ULong)a[i];
          ULong lo = (ULong)r[2 * i] + (sq & 0xFFFFFFFFULL) + carry;
          r[2 * i] = (DataT)(lo & 0xFFFFFFFFULL);
          ULong hi = (ULong)r[2 * i + 1] + (sq >> 32) + (lo >> 32);
          r[2 * i + 1] = (DataT)(hi & 0xFFFFFFFFULL);
          carry = hi >> 32;
        }
      }
      else if (base == Base2_64)
      {
        // Only reached with n == 0; nothing to do.
      }
      else
      {
        for (SizeT i = 0; i < n; ++i)
        {
          if (a[i] == 0) continue;
          ULong carry = 0;
          for (SizeT j = 0; j < n; ++j)
          {
            ULong prod = (ULong)a[i] * (ULong)a[j] + (ULong)r[i + j] + carry;
            r[i + j] = (DataT)(prod % base);
            carry = prod / base;
          }
          r[i + n] = (DataT)carry;
        }
      }
    }

    static std::vector<DataT> Square(std::vector<DataT> const &a, BaseT base)
    {
      return SquareLimbs(a, base);
    }

    static std::vector<DataT> SquareLimbs(std::span<const DataT> a, BaseT base)
    {
      if (IsZero(a))
        return std::vector<DataT>{0};
      SizeT n = (SizeT)a.size();
      std::vector<DataT> r(2 * n, 0);
      SquarePtr(a.data(), n, r.data(), base);
      TrimZeros(r);
      return r;
    }
  };
}

#endif
