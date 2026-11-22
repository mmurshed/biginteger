/**
 * BigMath: limb-level arithmetic kernels for full 64-bit limbs (Base2_64).
 *
 * The basecase loops of multiplication — schoolbook rows, Karatsuba's
 * add/subtract passes — all reduce to a handful of mpn-style primitives over
 * little-endian `DataT` arrays: add_n, sub_n, mul_1, addmul_1, submul_1,
 * shifts. This header is the one place they live. Schoolbook division
 * (FastDivision) and Burnikel–Ziegler use the same kernels.
 *
 * Why they are written this way. The obvious formulation
 *
 *     unsigned __int128 p = (unsigned __int128)a[j] * b + r[j] + carry;
 *     r[j] = (uint64_t)p;  carry = (uint64_t)(p >> 64);
 *
 * threads a *limb-valued* carry through every iteration: the next limb cannot
 * start its additions until the previous limb's high half is known, and each
 * step materializes a carry flag back into a register. On an M1 P-core that
 * costs 2.5-3.3 cycles per limb for addmul_1 and 2.2-3.3 for a plain add.
 *
 * The kernels below instead
 *   - work in straight-line blocks of four limbs, so the carry stays in the
 *     flags register from one limb to the next (an `adcs` chain) and is
 *     materialized once per block, and
 *   - for the multiply kernels, send the low and high halves of the products
 *     down two separate carry chains (lo[j] + hi[j-1], then += r[j]) instead
 *     of one chain through the multiplier.
 *
 * Measured on the same core: addmul_1 1.64-1.76 cycles/limb, add_n 1.1-1.7,
 * schoolbook n x n 1.5-1.77x faster — within 1.14-1.28x of GMP's hand-written
 * assembly basecase. No assembly here: the carry builtins compile to
 * adds/adcs on aarch64 and add/adc on x86-64.
 *
 * All kernels require true 64-bit limbs. The legacy 32-in-64 layout
 * (BIGMATH_LIMB_64=0) keeps its own loops.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGMATH_LIMB_KERNELS
#define BIGMATH_LIMB_KERNELS

#include <algorithm>

#include "Constants.h"

// Carry primitive selection. Clang and GCC >= 14 provide the multiprecision
// builtins; older GCC on x86-64 gets the ADX-style intrinsics, which it turns
// into the same adc chains; anything else falls back to two overflow checks.
#if defined(__clang__) || (defined(__GNUC__) && __GNUC__ >= 14)
#define BIGMATH_CARRY_BUILTINS 1
#elif defined(__x86_64__) && defined(__GNUC__)
#define BIGMATH_CARRY_X86 1
#include <immintrin.h>
#endif

namespace BigMath
{
  namespace LimbKernel
  {
    // a + b + carryIn (carryIn is 0 or 1); carryOut receives 0 or 1.
    static inline DataT Adc(DataT a, DataT b, DataT carryIn, DataT &carryOut)
    {
#if defined(BIGMATH_CARRY_BUILTINS)
      unsigned long long co;
      unsigned long long r = __builtin_addcll(a, b, carryIn, &co);
      carryOut = (DataT)co;
      return (DataT)r;
#elif defined(BIGMATH_CARRY_X86)
      unsigned long long r;
      carryOut = (DataT)_addcarry_u64((unsigned char)carryIn, a, b, &r);
      return (DataT)r;
#else
      DataT s, t;
      bool c1 = __builtin_add_overflow(a, b, &s);
      bool c2 = __builtin_add_overflow(s, carryIn, &t);
      carryOut = (DataT)(c1 | c2);
      return t;
#endif
    }

    // a - b - borrowIn (borrowIn is 0 or 1); borrowOut receives 0 or 1.
    static inline DataT Sbb(DataT a, DataT b, DataT borrowIn, DataT &borrowOut)
    {
#if defined(BIGMATH_CARRY_BUILTINS)
      unsigned long long bo;
      unsigned long long r = __builtin_subcll(a, b, borrowIn, &bo);
      borrowOut = (DataT)bo;
      return (DataT)r;
#elif defined(BIGMATH_CARRY_X86)
      unsigned long long r;
      borrowOut = (DataT)_subborrow_u64((unsigned char)borrowIn, a, b, &r);
      return (DataT)r;
#else
      DataT s, t;
      bool b1 = __builtin_sub_overflow(a, b, &s);
      bool b2 = __builtin_sub_overflow(s, borrowIn, &t);
      borrowOut = (DataT)(b1 | b2);
      return t;
#endif
    }

    // r[0..n) = a[0..n) + b[0..n) + carry. Returns the carry out. r may alias
    // a and/or b exactly (same start).
    static inline DataT AddN(DataT *r, const DataT *a, const DataT *b, SizeT n, DataT carry = 0)
    {
      SizeT i = 0;
      for (; i + 4 <= n; i += 4)
      {
        r[i] = Adc(a[i], b[i], carry, carry);
        r[i + 1] = Adc(a[i + 1], b[i + 1], carry, carry);
        r[i + 2] = Adc(a[i + 2], b[i + 2], carry, carry);
        r[i + 3] = Adc(a[i + 3], b[i + 3], carry, carry);
      }
      for (; i < n; ++i)
        r[i] = Adc(a[i], b[i], carry, carry);
      return carry;
    }

    // r[0..n) = a[0..n) - b[0..n) - borrow. Returns the borrow out. r may
    // alias a and/or b exactly.
    static inline DataT SubN(DataT *r, const DataT *a, const DataT *b, SizeT n, DataT borrow = 0)
    {
      SizeT i = 0;
      for (; i + 4 <= n; i += 4)
      {
        r[i] = Sbb(a[i], b[i], borrow, borrow);
        r[i + 1] = Sbb(a[i + 1], b[i + 1], borrow, borrow);
        r[i + 2] = Sbb(a[i + 2], b[i + 2], borrow, borrow);
        r[i + 3] = Sbb(a[i + 3], b[i + 3], borrow, borrow);
      }
      for (; i < n; ++i)
        r[i] = Sbb(a[i], b[i], borrow, borrow);
      return borrow;
    }

    // r[0..n) = a[0..n) + carry. Returns the carry out. Copies the tail once
    // the carry dies (unless r == a, where there is nothing left to do).
    static inline DataT AddCarry(DataT *r, const DataT *a, SizeT n, DataT carry)
    {
      SizeT i = 0;
      for (; i < n && carry; ++i)
      {
        DataT s = a[i] + carry;
        carry = s < carry ? 1 : 0;
        r[i] = s;
      }
      if (r != a)
        for (; i < n; ++i)
          r[i] = a[i];
      return carry;
    }

    // r[0..n) = a[0..n) - borrow. Returns the borrow out.
    static inline DataT SubBorrow(DataT *r, const DataT *a, SizeT n, DataT borrow)
    {
      SizeT i = 0;
      for (; i < n && borrow; ++i)
      {
        DataT ai = a[i];
        r[i] = ai - borrow;
        borrow = ai < borrow ? 1 : 0;
      }
      if (r != a)
        for (; i < n; ++i)
          r[i] = a[i];
      return borrow;
    }

    // r[0..lenR) = a[0..lenA) + b[0..lenB). Requires lenR >= max(lenA, lenB);
    // limbs above the longer operand receive the carry, then zeros. r must
    // not partially overlap the inputs.
    static inline void AddExt(DataT *r, SizeT lenR, const DataT *a, SizeT lenA, const DataT *b, SizeT lenB)
    {
      if (lenA < lenB)
      {
        std::swap(a, b);
        std::swap(lenA, lenB);
      }
      DataT carry = AddN(r, a, b, lenB);
      carry = AddCarry(r + lenB, a + lenB, lenA - lenB, carry);
      for (SizeT i = lenA; i < lenR; ++i)
      {
        r[i] = carry;
        carry = 0;
      }
    }

    // dest[0..destLen) += src[0..srcLen). Returns the carry out of the top
    // limb. srcLen beyond destLen is ignored.
    static inline DataT AddInPlace(DataT *dest, SizeT destLen, const DataT *src, SizeT srcLen)
    {
      SizeT shared = std::min(srcLen, destLen);
      DataT carry = AddN(dest, dest, src, shared);
      return AddCarry(dest + shared, dest + shared, destLen - shared, carry);
    }

    // dest[0..destLen) -= src[0..srcLen). Returns the borrow out of the top
    // limb (0 whenever dest >= src). srcLen beyond destLen is ignored.
    static inline DataT SubInPlace(DataT *dest, SizeT destLen, const DataT *src, SizeT srcLen)
    {
      SizeT shared = std::min(srcLen, destLen);
      DataT borrow = SubN(dest, dest, src, shared);
      return SubBorrow(dest + shared, dest + shared, destLen - shared, borrow);
    }

    // r[0..n) = a[0..n) * b. Returns the high limb. r may alias a exactly.
    static inline DataT Mul1(DataT *r, const DataT *a, SizeT n, DataT b)
    {
      DataT high = 0; // high half of the previous product, plus its carry
      SizeT j = 0;
      for (; j + 4 <= n; j += 4)
      {
        ULong128 p0 = (ULong128)a[j] * b;
        ULong128 p1 = (ULong128)a[j + 1] * b;
        ULong128 p2 = (ULong128)a[j + 2] * b;
        ULong128 p3 = (ULong128)a[j + 3] * b;
        DataT c;
        r[j] = Adc((DataT)p0, high, 0, c);
        r[j + 1] = Adc((DataT)p1, (DataT)(p0 >> 64), c, c);
        r[j + 2] = Adc((DataT)p2, (DataT)(p1 >> 64), c, c);
        r[j + 3] = Adc((DataT)p3, (DataT)(p2 >> 64), c, c);
        high = (DataT)(p3 >> 64) + c; // cannot overflow: a high half is <= 2^64 - 2
      }
      for (; j < n; ++j)
      {
        ULong128 p = (ULong128)a[j] * b + high;
        r[j] = (DataT)p;
        high = (DataT)(p >> 64);
      }
      return high;
    }

    // r[0..n) += a[0..n) * b. Returns the carry limb. r must not overlap a.
    static inline DataT AddMul1(DataT *r, const DataT *a, SizeT n, DataT b)
    {
      DataT high = 0;   // chain A: high half of the previous product (+ its carry)
      DataT carryR = 0; // chain B: carry of the accumulation into r
      SizeT j = 0;
      for (; j + 4 <= n; j += 4)
      {
        ULong128 p0 = (ULong128)a[j] * b;
        ULong128 p1 = (ULong128)a[j + 1] * b;
        ULong128 p2 = (ULong128)a[j + 2] * b;
        ULong128 p3 = (ULong128)a[j + 3] * b;
        // Chain A: t_k = lo(p_k) + hi(p_{k-1}).
        DataT c;
        DataT t0 = Adc((DataT)p0, high, 0, c);
        DataT t1 = Adc((DataT)p1, (DataT)(p0 >> 64), c, c);
        DataT t2 = Adc((DataT)p2, (DataT)(p1 >> 64), c, c);
        DataT t3 = Adc((DataT)p3, (DataT)(p2 >> 64), c, c);
        high = (DataT)(p3 >> 64) + c;
        // Chain B: r_k += t_k.
        r[j] = Adc(r[j], t0, carryR, carryR);
        r[j + 1] = Adc(r[j + 1], t1, carryR, carryR);
        r[j + 2] = Adc(r[j + 2], t2, carryR, carryR);
        r[j + 3] = Adc(r[j + 3], t3, carryR, carryR);
      }
      for (; j < n; ++j)
      {
        ULong128 p = (ULong128)a[j] * b;
        DataT c;
        DataT t = Adc((DataT)p, high, 0, c);
        high = (DataT)(p >> 64) + c;
        r[j] = Adc(r[j], t, carryR, carryR);
      }
      // r + a*b < B^(n+1), so the carry limb high + carryR cannot overflow.
      return high + carryR;
    }

    // r[0..n) -= a[0..n) * b. Returns the borrow limb (what must still be
    // subtracted from r[n]). r must not overlap a. Same two-chain shape as
    // AddMul1: the products are formed along an add chain, then subtracted
    // from r along a borrow chain — no data-dependent branch per limb.
    static inline DataT SubMul1(DataT *r, const DataT *a, SizeT n, DataT b)
    {
      DataT high = 0;    // chain A: high half of the previous product (+ its carry)
      DataT borrowR = 0; // chain B: borrow of the subtraction from r
      SizeT j = 0;
      for (; j + 4 <= n; j += 4)
      {
        ULong128 p0 = (ULong128)a[j] * b;
        ULong128 p1 = (ULong128)a[j + 1] * b;
        ULong128 p2 = (ULong128)a[j + 2] * b;
        ULong128 p3 = (ULong128)a[j + 3] * b;
        DataT c;
        DataT t0 = Adc((DataT)p0, high, 0, c);
        DataT t1 = Adc((DataT)p1, (DataT)(p0 >> 64), c, c);
        DataT t2 = Adc((DataT)p2, (DataT)(p1 >> 64), c, c);
        DataT t3 = Adc((DataT)p3, (DataT)(p2 >> 64), c, c);
        high = (DataT)(p3 >> 64) + c;
        r[j] = Sbb(r[j], t0, borrowR, borrowR);
        r[j + 1] = Sbb(r[j + 1], t1, borrowR, borrowR);
        r[j + 2] = Sbb(r[j + 2], t2, borrowR, borrowR);
        r[j + 3] = Sbb(r[j + 3], t3, borrowR, borrowR);
      }
      for (; j < n; ++j)
      {
        ULong128 p = (ULong128)a[j] * b;
        DataT c;
        DataT t = Adc((DataT)p, high, 0, c);
        high = (DataT)(p >> 64) + c;
        r[j] = Sbb(r[j], t, borrowR, borrowR);
      }
      // a*b < B^(n+1) and r < B^n, so the borrow limb high + borrowR <= B - 1.
      return high + borrowR;
    }

    // r[0..n) = a[0..n) << count, 0 < count < 64. Returns the bits shifted
    // out of the top limb. n >= 1. r may alias a exactly, or sit above it.
    static inline DataT Lshift(DataT *r, const DataT *a, SizeT n, unsigned count)
    {
      const unsigned back = 64 - count;
      DataT high = a[n - 1];
      const DataT out = high >> back;
      for (SizeT i = n - 1; i > 0; --i)
      {
        DataT low = a[i - 1];
        r[i] = (high << count) | (low >> back);
        high = low;
      }
      r[0] = high << count;
      return out;
    }

    // r[0..n) = a[0..n) >> count, 0 < count < 64. Returns the bits shifted
    // out of the bottom limb, left-aligned. n >= 1. r may alias a exactly,
    // or sit below it.
    static inline DataT Rshift(DataT *r, const DataT *a, SizeT n, unsigned count)
    {
      const unsigned back = 64 - count;
      DataT low = a[0];
      const DataT out = low << back;
      for (SizeT i = 0; i + 1 < n; ++i)
      {
        DataT high = a[i + 1];
        r[i] = (low >> count) | (high << back);
        low = high;
      }
      r[n - 1] = low >> count;
      return out;
    }

    // A fixed row width lets the compiler schedule the carry chains across
    // the whole schoolbook leaf and eliminate the loop and remainder tests.
    template <SizeT Width>
    static inline void MulFixedRows(DataT *r, const DataT *a, const DataT *b, SizeT lb)
    {
      r[Width] = Mul1(r, a, Width, b[0]);
      for (SizeT i = 1; i < lb; ++i)
        r[Width + i] = AddMul1(r + i, a, Width, b[i]);
    }

    // r[0..la+lb) = a[0..la) * b[0..lb), schoolbook. la, lb >= 1. r must not
    // overlap a or b and need not be initialized. The inner loop runs over
    // the longer operand.
    static inline void MulBasecase(DataT *r, const DataT *a, SizeT la, const DataT *b, SizeT lb)
    {
      if (la < lb)
      {
        std::swap(a, b);
        std::swap(la, lb);
      }
#if defined(__APPLE__) && defined(__aarch64__) && defined(__clang__)
      switch (la)
      {
        case 12: return MulFixedRows<12>(r, a, b, lb);
        case 13: return MulFixedRows<13>(r, a, b, lb);
        case 14: return MulFixedRows<14>(r, a, b, lb);
        case 15: return MulFixedRows<15>(r, a, b, lb);
        case 16: return MulFixedRows<16>(r, a, b, lb);
        case 17: return MulFixedRows<17>(r, a, b, lb);
        case 18: return MulFixedRows<18>(r, a, b, lb);
        case 19: return MulFixedRows<19>(r, a, b, lb);
        case 20: return MulFixedRows<20>(r, a, b, lb);
        case 21: return MulFixedRows<21>(r, a, b, lb);
        case 22: return MulFixedRows<22>(r, a, b, lb);
        case 23: return MulFixedRows<23>(r, a, b, lb);
        case 24: return MulFixedRows<24>(r, a, b, lb);
      }
#endif
      r[la] = Mul1(r, a, la, b[0]);
      for (SizeT i = 1; i < lb; ++i)
        r[la + i] = AddMul1(r + i, a, la, b[i]);
    }

    // r[0..2n) = a[0..n)^2, schoolbook: off-diagonal products once, doubled,
    // plus the diagonal squares. n >= 1. r must not overlap a.
    static inline void SqrBasecase(DataT *r, const DataT *a, SizeT n)
    {
      // Off-diagonal triangle into r[1 .. 2n-1): row i adds a[i] * a[i+1..n).
      r[0] = 0;
      if (n > 1)
      {
        r[n] = Mul1(r + 1, a + 1, n - 1, a[0]);
        for (SizeT i = 1; i + 1 < n; ++i)
          r[n + i] = AddMul1(r + 2 * i + 1, a + i + 1, n - i - 1, a[i]);
      }
      r[2 * n - 1] = 0;

      // r = 2*r + sum a[i]^2 * B^(2i), one pass: the doubling carry and the
      // addition carry are separate chains.
      DataT shiftIn = 0;
      DataT carry = 0;
      for (SizeT i = 0; i < n; ++i)
      {
        ULong128 sq = (ULong128)a[i] * a[i];
        DataT lo = r[2 * i];
        DataT hi = r[2 * i + 1];
        DataT dlo = (lo << 1) | shiftIn;
        DataT dhi = (hi << 1) | (lo >> 63);
        shiftIn = hi >> 63;
        r[2 * i] = Adc(dlo, (DataT)sq, carry, carry);
        r[2 * i + 1] = Adc(dhi, (DataT)(sq >> 64), carry, carry);
      }
    }
  }
}

#endif
