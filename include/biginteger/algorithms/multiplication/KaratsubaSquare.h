/**
 * BigInteger Class
 * Karatsuba squaring: a² = al² + ((al+ah)² − al² − ah²) · B^m + ah² · B^(2m).
 * Pointer-based with shared workspace, mirroring KaratsubaMultiplication.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef KARATSUBA_SQUARE
#define KARATSUBA_SQUARE

#include <vector>
#include <algorithm>
#include <cstring>
#include <memory>

#include "../../common/LimbKernels.h"
#include "../../common/Util.h"
#include "ClassicSquare.h"
#include "ToomKernels.h"

namespace BigMath
{
  class KaratsubaSquare
  {
  private:
#ifndef BIGMATH_KARATSUBA_SQUARE_THRESHOLD
#define BIGMATH_KARATSUBA_SQUARE_THRESHOLD 48
#endif
    static const SizeT THRESHOLD = BIGMATH_KARATSUBA_SQUARE_THRESHOLD;
// Toom-3 split for squares (64-bit limbs): operands of at least this many
// limbs are split in three (five squares of about n/3 limbs).
#ifndef BIGMATH_SQR_TOOM3_THRESHOLD
#define BIGMATH_SQR_TOOM3_THRESHOLD 160 // canonical value: build/DispatchThresholds.h
#endif
    static const SizeT TOOM3_THRESHOLD = BIGMATH_SQR_TOOM3_THRESHOLD;

    // r[0..lenR-1] = a[0..lenA-1] + b[0..lenB-1].
    static void AddPtr(
        const DataT *a, SizeT lenA,
        const DataT *b, SizeT lenB,
        DataT *r, SizeT lenR,
        BaseT base)
    {
      if (base == Base2_64)
      {
        // Flag-chain kernels — see common/LimbKernels.h. Every call site
        // passes lenR >= max(lenA, lenB).
        LimbKernel::AddExt(r, lenR, a, lenA, b, lenB);
        return;
      }
      ULong carry = 0;
      for (SizeT i = 0; i < lenR; ++i)
      {
        ULong sum = carry;
        if (i < lenA) sum += a[i];
        if (i < lenB) sum += b[i];
        if (base == Base2_32)
        {
          r[i] = (DataT)(sum & 0xFFFFFFFFULL);
          carry = sum >> 32;
        }
        else
        {
          r[i] = (DataT)(sum % base);
          carry = sum / base;
        }
      }
    }

    // dest[0..destLen-1] += src[0..srcLen-1]
    static void AddToPtr(
        DataT *dest, SizeT destLen,
        const DataT *src, SizeT srcLen,
        BaseT base)
    {
      if (base == Base2_64)
      {
        LimbKernel::AddInPlace(dest, destLen, src, srcLen);
        return;
      }
      ULong carry = 0;
      for (SizeT i = 0; i < destLen; ++i)
      {
        ULong sum = (ULong)dest[i] + carry;
        if (i < srcLen) sum += src[i];
        if (base == Base2_32)
        {
          dest[i] = (DataT)(sum & 0xFFFFFFFFULL);
          carry = sum >> 32;
        }
        else
        {
          dest[i] = (DataT)(sum % base);
          carry = sum / base;
        }
        if (carry == 0 && i >= srcLen) break;
      }
    }

    // dest[0..destLen-1] -= src[0..srcLen-1]  (assumes dest >= src)
    static void SubFromPtr(
        DataT *dest, SizeT destLen,
        const DataT *src, SizeT srcLen,
        BaseT base)
    {
      if (base == Base2_64)
      {
        LimbKernel::SubInPlace(dest, destLen, src, srcLen);
        return;
      }
      Long borrow = 0;
      for (SizeT i = 0; i < destLen; ++i)
      {
        Long diff = (Long)dest[i] - borrow;
        if (i < srcLen) diff -= (Long)src[i];
        if (diff < 0) { diff += base; borrow = 1; }
        else borrow = 0;
        dest[i] = (DataT)diff;
        if (borrow == 0 && i >= srcLen) break;
      }
    }

    // c[0..2n-1] = a[0..n-1]^2.  w is workspace (≥ ~4n limbs, see public Square).
    static void SquareRec(
        const DataT *a, SizeT n,
        DataT *c,
        DataT *w,
        BaseT base)
    {
      if (n <= THRESHOLD)
      {
        ClassicSquare::SquarePtr(a, n, c, base);
        return;
      }

      if (base == Base2_64 && n >= TOOM3_THRESHOLD && n >= 6)
      {
        // Toom-3: evaluate once, square at 0, 1, -1, 2, inf (the value at
        // -1 is squared, so its sign is +). Workspace 9k + 9 here.
        const SizeT k = (n + 2) / 3, r = n - 2 * k;
        const SizeT e = k + 1, L = 2 * e;
        DataT *s1 = w, *sm1 = s1 + e, *s2 = sm1 + e;
        DataT *v1 = s2 + e, *vm1 = v1 + L, *v2 = vm1 + L, *next = v2 + L;
        ToomKernel::Toom3Evaluate(a, k, r, s1, sm1, s2);
        SquareRec(a, k, c, next, base);
        SquareRec(a + 2 * k, r, c + 4 * k, next, base);
        SquareRec(s1, e, v1, next, base);
        SquareRec(sm1, e, vm1, next, base);
        SquareRec(s2, e, v2, next, base);
        ToomKernel::Toom3Interpolate(c, k, 2 * r, 2 * n, v1, vm1, v2, L, 1);
        return;
      }

      SizeT m = (n + 1) / 2;
      SizeT lenAl = m;
      SizeT lenAh = n - m;
      SizeT lenSum = m + 1;
      SizeT lenPs = 2 * lenSum;

      // Layout c (size 2n):
      //   c[0..2m-1]  = al² (p0)
      //   c[2m..2n-1] = ah² (p2)
      // Then add (ps − p0 − p2) shifted by m into c[m..].
      std::memset(c, 0, 2 * n * sizeof(DataT));

      // 1. p0 = al² into c[0..2m-1]
      SquareRec(a, lenAl, c, w, base);

      // 2. p2 = ah² into c[2m..2n-1]
      SquareRec(a + m, lenAh, c + 2 * m, w, base);

      // Workspace layout: sum (lenSum) | psBuf (lenPs) | nextW...
      DataT *sum   = w;
      DataT *psBuf = w + lenSum;
      DataT *nextW = psBuf + lenPs;

      // 3. sum = al + ah  (lenSum limbs incl. carry)
      AddPtr(a, lenAl, a + m, lenAh, sum, lenSum, base);

      // 4. psBuf = sum²  (size lenPs)
      SquareRec(sum, lenSum, psBuf, nextW, base);

      // 5. psBuf -= p0 ; psBuf -= p2     (result = 2·al·ah ≥ 0)
      SubFromPtr(psBuf, lenPs, c, 2 * m, base);
      SubFromPtr(psBuf, lenPs, c + 2 * m, 2 * lenAh, base);

      // 6. c[m..2n-1] += psBuf
      AddToPtr(c + m, 2 * n - m, psBuf, lenPs, base);
    }

  public:
    static std::vector<DataT> Square(std::vector<DataT> const &a, BaseT base)
    {
      if (IsZero(a))
        return std::vector<DataT>{0};

      SizeT n = (SizeT)a.size();
      if (n <= THRESHOLD)
        return ClassicSquare::Square(a, base);

      std::vector<DataT> c(2 * n, 0);

      // Workspace bound: per level uses 3m+3 ≈ 1.5n; recursion sum ≈ 3n. Use 8n for safety
      // (matches KaratsubaMultiplication).
      std::unique_ptr<DataT[]> w(new DataT[8 * n]);

      SquareRec(a.data(), n, c.data(), w.get(), base);
      TrimZeros(c);
      return c;
    }
  };
}

#endif
