/**
 * BigMath: Toom-Cook evaluation and interpolation kernels on raw limbs.
 *
 * Shared by KaratsubaMultiplication (products) and KaratsubaSquare
 * (squares), which split a recursion level three or four ways at the
 * BIGMATH_MUL_TOOM33/44 and BIGMATH_SQR_TOOM3 thresholds. 64-bit limbs only.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGMATH_TOOM_KERNELS
#define BIGMATH_TOOM_KERNELS

#include <algorithm>
#include <cstring>

#include "../../common/LimbKernels.h"

namespace BigMath
{
  namespace ToomKernel
  {
    // r[0..n) = a[0..n) + 2 * b[0..n); returns the carry limb (0..2).
    inline DataT AddLsh1N(DataT *r, const DataT *a, const DataT *b, SizeT n)
    {
        DataT shifted = 0, carry = 0;
        for (SizeT i = 0; i < n; ++i)
        {
            const DataT bi = b[i];
            const DataT twice = (bi << 1) | shifted;
            shifted = bi >> 63;
            r[i] = LimbKernel::Adc(a[i], twice, carry, carry);
        }
        return shifted + carry;
    }

    // r[0..n) /= 3, exact (r must be a multiple of 3). Multiplies by the
    // inverse of 3 mod 2^64, carrying the borrow of each limb's
    // quotient into the next (GMP's divexact_by3).
    inline void DivExact3(DataT *r, SizeT n)
    {
        const DataT inv3 = 0xAAAAAAAAAAAAAAABULL;
        DataT c = 0;
        for (SizeT i = 0; i < n; ++i)
        {
            const DataT s = r[i];
            const DataT l = s - c;
            const DataT q = l * inv3;
            r[i] = q;
            c = (DataT)(l > s) + (DataT)(q > 0x5555555555555555ULL) + (DataT)(q > 0xAAAAAAAAAAAAAAAAULL);
        }
    }

    // Evaluate x = x0 + x1*X + x2*X^2 (x0, x1: k limbs; x2: r limbs,
    // 1 <= r <= k) at X = 1, -1 and 2, each into k + 1 limbs. Returns
    // the sign of x(-1) (+1 or -1); sm1 holds its magnitude.
    inline int Toom3Evaluate(const DataT *x, SizeT k, SizeT r, DataT *s1, DataT *sm1, DataT *s2)
    {
        const DataT *x0 = x, *x1 = x + k, *x2 = x + 2 * k;
        // t = x0 + x2 (k + 1 limbs), parked in s2 until s1 and sm1 are done.
        DataT *t = s2;
        DataT c = LimbKernel::AddN(t, x0, x2, r);
        if (k > r)
            c = LimbKernel::AddCarry(t + r, x0 + r, k - r, c);
        t[k] = c;
        // s1 = t + x1
        s1[k] = t[k] + LimbKernel::AddN(s1, t, x1, k);
        // sm1 = |t - x1|
        int sign = 1;
        bool less = t[k] == 0;
        if (less)
        {
            SizeT i = k;
            while (i > 0 && t[i - 1] == x1[i - 1])
                --i;
            less = i > 0 && t[i - 1] < x1[i - 1];
        }
        if (less)
        {
            LimbKernel::SubN(sm1, x1, t, k);
            sm1[k] = 0;
            sign = -1;
        }
        else
        {
            const DataT borrow = LimbKernel::SubN(sm1, t, x1, k);
            sm1[k] = t[k] - borrow;
        }
        // s2 = x0 + 2*(x1 + 2*x2): u = x1 + 2*x2 into s2, then s2 = x0 + 2*u.
        DataT *u = s2;
        DataT cu = AddLsh1N(u, x1, x2, r);
        if (k > r)
            cu = LimbKernel::AddCarry(u + r, x1 + r, k - r, cu);
        u[k] = cu;
        const DataT top = u[k];
        s2[k] = (top << 1) + AddLsh1N(s2, x0, u, k);
        return sign;
    }

    // r[0..n) = x + (y << s) with x zero-extended past lx and y past ly
    // (lx, ly <= n, 0 <= s < 64). The value must fit in n limbs.
    inline void AddShifted(DataT *r, SizeT n, const DataT *x, SizeT lx, const DataT *y, SizeT ly, unsigned s)
    {
        DataT prev = 0, carry = 0;
        for (SizeT i = 0; i < n; ++i)
        {
            const DataT xi = i < lx ? x[i] : 0;
            const DataT yi = i < ly ? y[i] : 0;
            const DataT sh = s ? (yi << s) | (prev >> (64 - s)) : yi;
            prev = yi;
            r[i] = LimbKernel::Adc(xi, sh, carry, carry);
        }
    }

    // r[0..n) = |x - y| for n-limb x, y; returns +1 if x >= y, else -1.
    inline int AbsDiff(DataT *r, const DataT *x, const DataT *y, SizeT n)
    {
        SizeT i = n;
        while (i > 0 && x[i - 1] == y[i - 1])
            --i;
        if (i > 0 && x[i - 1] < y[i - 1])
        {
            LimbKernel::SubN(r, y, x, n);
            return -1;
        }
        LimbKernel::SubN(r, x, y, n);
        return 1;
    }

    // r[0..n) /= 5, exact.
    inline void DivExact5(DataT *r, SizeT n)
    {
        const DataT inv5 = 0xCCCCCCCCCCCCCCCDULL; // 5 * inv5 == 1 mod 2^64
        DataT c = 0;
        for (SizeT i = 0; i < n; ++i)
        {
            const DataT s = r[i];
            const DataT l = s - c;
            const DataT q = l * inv5;
            r[i] = q;
            c = (DataT)(l > s) + (DataT)(((ULong128)q * 5) >> 64);
        }
    }

    // dst[0..L) -= (src[0..ls) << s), src zero-extended (value fits).
    inline void SubShifted(DataT *dst, SizeT L, const DataT *src, SizeT ls, unsigned s, DataT *tmp)
    {
        AddShifted(tmp, L, nullptr, 0, src, ls, s);
        LimbKernel::SubN(dst, dst, tmp, L);
    }

    // Evaluate x = x0 + x1 X + x2 X^2 + x3 X^3 (x0..x2: k limbs, x3: r
    // limbs) at 1, -1, 2, -2 and 1/2 (as 8 x(1/2)), each into e = k + 1
    // limbs. sm1/sm2 hold |x(-1)|, |x(-2)|; the signs are returned via
    // sign1/sign2. t1, t2 are e-limb scratch.
    inline void Toom4Evaluate(const DataT *x, SizeT k, SizeT r,
                              DataT *p1, DataT *pm1, DataT *p2, DataT *pm2, DataT *ph,
                              int &sign1, int &sign2, DataT *t1, DataT *t2)
    {
        const SizeT e = k + 1;
        const DataT *x0 = x, *x1 = x + k, *x2 = x + 2 * k, *x3 = x + 3 * k;
        // +-1: t1 = x0 + x2, t2 = x1 + x3
        AddShifted(t1, e, x0, k, x2, k, 0);
        AddShifted(t2, e, x1, k, x3, r, 0);
        LimbKernel::AddN(p1, t1, t2, e);
        sign1 = AbsDiff(pm1, t1, t2, e);
        // +-2: t1 = x0 + 4 x2, t2 = 2 x1 + 8 x3 = 2 (x1 + 4 x3)
        AddShifted(t1, e, x0, k, x2, k, 2);
        AddShifted(ph, e, x1, k, x3, r, 2);       // ph as scratch
        AddShifted(t2, e, nullptr, 0, ph, e, 1);
        LimbKernel::AddN(p2, t1, t2, e);
        sign2 = AbsDiff(pm2, t1, t2, e);
        // 8 x(1/2) = ((2 x0 + x1) 2 + x2) 2 + x3
        AddShifted(t1, e, x1, k, x0, k, 1);
        AddShifted(t2, e, x2, k, t1, e, 1);
        AddShifted(ph, e, x3, r, t2, e, 1);
    }

    // Toom-3 interpolation (Bodrato), points 0, 1, -1, 2, inf. On entry
    // c[0, 2k) = v0 and c[4k, 4k + linf) = vinf, and v1, vm1, v2 hold the
    // L-limb products at 1, -1 (magnitude; `sign` is its sign) and 2.
    // Leaves the product in c[0, total). With the product's coefficients
    // r0..r4:
    //   v2  <- (v2 - vm1) / 3   = r1 + r2 + 3 r3 + 5 r4
    //   vm1 <- (v1 - vm1) / 2   = r1 + r3
    //   v1  <- v1 - v0          = r1 + r2 + r3 + r4
    //   v2  <- (v2 - v1) / 2    = r3 + 2 r4
    //   v1  <- v1 - vm1 - vinf  = r2
    //   v2  <- v2 - 2 vinf      = r3
    //   vm1 <- vm1 - v2         = r1
    // Every intermediate is a non-negative combination of the r_i, so the
    // sequence runs on unsigned limbs.
    inline void Toom3Interpolate(DataT *c, SizeT k, SizeT linf, SizeT total,
                                 DataT *v1, DataT *vm1, DataT *v2, SizeT L, int sign)
    {
      DataT *v0 = c, *vinf = c + 4 * k;
      const SizeT lv0 = 2 * k;
      // 1. v2 = (v2 - vm1) / 3
      if (sign > 0)
          LimbKernel::SubN(v2, v2, vm1, L);
      else
          LimbKernel::AddN(v2, v2, vm1, L);
      DivExact3(v2, L);
      // 2. vm1 = (v1 - vm1) / 2
      if (sign > 0)
          LimbKernel::SubN(vm1, v1, vm1, L);
      else
          LimbKernel::AddN(vm1, v1, vm1, L);
      LimbKernel::Rshift(vm1, vm1, L, 1);
      // 3. v1 -= v0
      LimbKernel::SubInPlace(v1, L, v0, lv0);
      // 4. v2 = (v2 - v1) / 2
      LimbKernel::SubN(v2, v2, v1, L);
      LimbKernel::Rshift(v2, v2, L, 1);
      // 5. v1 -= vm1 + vinf
      LimbKernel::SubN(v1, v1, vm1, L);
      LimbKernel::SubInPlace(v1, L, vinf, linf);
      // 6. v2 -= 2 vinf
      LimbKernel::SubInPlace(v2, L, vinf, linf);
      LimbKernel::SubInPlace(v2, L, vinf, linf);
      // 7. vm1 -= v2
      LimbKernel::SubN(vm1, vm1, v2, L);

      // c = v0 + r1 B^k + r2 B^2k + r3 B^3k + vinf B^4k. c[2k, 4k) is
      // still unwritten. Each r_i fits below the top of c, so limbs of
      // the L-limb buffers that would land past it are zero.
      std::memset(c + 2 * k, 0, 2 * k * sizeof(DataT));
      LimbKernel::AddInPlace(c + k, total - k, vm1, std::min<SizeT>(L, total - k));
      LimbKernel::AddInPlace(c + 2 * k, total - 2 * k, v1, std::min<SizeT>(L, total - 2 * k));
      LimbKernel::AddInPlace(c + 3 * k, total - 3 * k, v2, std::min<SizeT>(L, total - 3 * k));
    }

    // Toom-4 interpolation, points 0, +-1, +-2, 1/2, inf. On entry
    // c[0, 2k) = r0 and c[6k, 6k + l6) = r6; P1, Pm1, P2, Pm2, Q hold the
    // L-limb products at 1, -1, 2, -2 (magnitudes; s1, s2 their signs) and
    // 64 P(1/2). T1..T3 are L-limb scratch. With coefficients r0..r6:
    //   E1 = (P(1) + P(-1)) / 2 - r0 - r6        = r2 + r4
    //   O1 = (P(1) - P(-1)) / 2                  = r1 + r3 + r5
    //   E2 = (P(2) + P(-2)) / 2 - r0 - 64 r6     = 4 r2 + 16 r4
    //   O2 = (P(2) - P(-2)) / 4                  = r1 + 4 r3 + 16 r5
    //   r4 = (E2 - 4 E1) / 12,  r2 = E1 - r4
    //   Q  = (64 P(1/2) - 64 r0 - 16 r2 - 4 r4 - r6) / 2 = 16 r1 + 4 r3 + r5
    //   X  = (O2 - O1) / 3        = r3 + 5 r5
    //   Y  = (16 O1 - Q) / 3      = 4 r3 + 5 r5
    //   r3 = (Y - X) / 3,  r5 = (X - r3) / 5,  r1 = O1 - r3 - r5
    // The signs of P(-1), P(-2) matter only in the first sums; every later
    // value is a non-negative combination of the r_i.
    inline void Toom4Interpolate(DataT *c, SizeT k, SizeT l6, SizeT total,
                                 DataT *P1, DataT *Pm1, DataT *P2, DataT *Pm2, DataT *Q,
                                 DataT *T1, DataT *T2, DataT *T3, SizeT L, int s1, int s2)
    {
      DataT *r0 = c, *r6 = c + 6 * k;
      const SizeT l0 = 2 * k;
      // E1 -> T1, O1 -> P1
      if (s1 > 0) { LimbKernel::AddN(T1, P1, Pm1, L); LimbKernel::SubN(P1, P1, Pm1, L); }
      else        { LimbKernel::SubN(T1, P1, Pm1, L); LimbKernel::AddN(P1, P1, Pm1, L); }
      LimbKernel::Rshift(T1, T1, L, 1);
      LimbKernel::Rshift(P1, P1, L, 1);
      LimbKernel::SubInPlace(T1, L, r0, l0);
      LimbKernel::SubInPlace(T1, L, r6, l6);
      // E2 -> T2, O2 -> P2
      if (s2 > 0) { LimbKernel::AddN(T2, P2, Pm2, L); LimbKernel::SubN(P2, P2, Pm2, L); }
      else        { LimbKernel::SubN(T2, P2, Pm2, L); LimbKernel::AddN(P2, P2, Pm2, L); }
      LimbKernel::Rshift(T2, T2, L, 1);
      LimbKernel::Rshift(P2, P2, L, 2);
      LimbKernel::SubInPlace(T2, L, r0, l0);
      SubShifted(T2, L, r6, l6, 6, T3);
      // r4 = (E2 - 4 E1) / 12 -> Pm2;  r2 = E1 - r4 -> T1
      SubShifted(T2, L, T1, L, 2, T3);
      LimbKernel::Rshift(T2, T2, L, 2);
      DivExact3(T2, L);
      DataT *R4 = T2, *R2 = T1;
      LimbKernel::SubN(R2, R2, R4, L);
      // Q = (Q - 64 r0 - 16 r2 - 4 r4 - r6) / 2
      SubShifted(Q, L, r0, l0, 6, T3);
      SubShifted(Q, L, R2, L, 4, T3);
      SubShifted(Q, L, R4, L, 2, T3);
      LimbKernel::SubInPlace(Q, L, r6, l6);
      LimbKernel::Rshift(Q, Q, L, 1);
      // X = (O2 - O1) / 3 -> P2;  Y = (16 O1 - Q) / 3 -> Pm1
      LimbKernel::SubN(P2, P2, P1, L);
      DivExact3(P2, L);
      AddShifted(Pm1, L, nullptr, 0, P1, L, 4);
      LimbKernel::SubN(Pm1, Pm1, Q, L);
      DivExact3(Pm1, L);
      // r3 = (Y - X) / 3 -> Pm1;  r5 = (X - r3) / 5 -> P2;  r1 = O1 - r3 - r5 -> P1
      LimbKernel::SubN(Pm1, Pm1, P2, L);
      DivExact3(Pm1, L);
      LimbKernel::SubN(P2, P2, Pm1, L);
      DivExact5(P2, L);
      LimbKernel::SubN(P1, P1, Pm1, L);
      LimbKernel::SubN(P1, P1, P2, L);
      DataT *R1 = P1, *R3 = Pm1, *R5 = P2;

      std::memset(c + 2 * k, 0, 4 * k * sizeof(DataT));
      const DataT *mid[5] = {R1, R2, R3, R4, R5};
      for (SizeT i = 1; i <= 5; ++i)
          LimbKernel::AddInPlace(c + i * k, total - i * k, mid[i - 1], std::min<SizeT>(L, total - i * k));
    }
  }
}

#endif
