/**
 * BigInteger Class
 * Version 12.0
 * Toom-3 multiplication: 5 evaluation points {0, 1, -1, 2, ∞}.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef TOOMCOOK_MULTIPLICATION
#define TOOMCOOK_MULTIPLICATION

#include <vector>
#include <algorithm>
using namespace std;

#include "../../common/Util.h"
#include "../../common/Comparator.h"
#include "../Addition.h"
#include "../Subtraction.h"
#include "../Shift.h"
#include "../multiplication/ClassicMultiplication.h"
#include "../multiplication/KaratsubaMultiplication.h"

namespace BigMath
{
  class ToomCookMultiplication
  {
  private:
#ifndef BIGMATH_TOOM3_THRESHOLD
#define BIGMATH_TOOM3_THRESHOLD 256
#endif
    static const SizeT TOOM3_THRESHOLD = BIGMATH_TOOM3_THRESHOLD;

    struct Signed
    {
      vector<DataT> mag;
      int sign;
    };

    static Signed MakePositive(vector<DataT> v)
    {
      TrimZeros(v);
      return {std::move(v), 1};
    }

    static Signed SubSigned(vector<DataT> const &a, vector<DataT> const &b, BaseT base)
    {
      Int c = Compare(a, b);
      if (c == 0) return {vector<DataT>{0}, 1};
      if (c > 0) return {Subtract(a, b, base), 1};
      return {Subtract(b, a, base), -1};
    }

    static Signed AddSigned(Signed const &a, Signed const &b, BaseT base)
    {
      if (IsZero(a.mag)) return b;
      if (IsZero(b.mag)) return a;
      if (a.sign == b.sign)
        return {Add(a.mag, b.mag, base), a.sign};
      Int c = Compare(a.mag, b.mag);
      if (c == 0) return {vector<DataT>{0}, 1};
      if (c > 0) return {Subtract(a.mag, b.mag, base), a.sign};
      return {Subtract(b.mag, a.mag, base), b.sign};
    }

    static Signed SubSigned(Signed const &a, Signed const &b, BaseT base)
    {
      Signed nb = b;
      if (!IsZero(nb.mag)) nb.sign = -nb.sign;
      return AddSigned(a, nb, base);
    }

    // In-place divisions used during interpolation.
    static void HalveInPlace(vector<DataT> &a, BaseT base)
    {
      if (a.empty()) return;
      if (base == Base2_32)
      {
        ULong carry = 0;
        for (Int i = (Int)a.size() - 1; i >= 0; --i)
        {
          ULong v = ((ULong)a[i]) | (carry << 32);
          a[i] = (DataT)(v >> 1);
          carry = v & 1ULL;
        }
      }
      else if (base == Base2_64)
      {
        // carry is 0 or 1 (v % 2). v = (carry << 64) | a[i] fits ULong128.
        ULong carry = 0;
        for (Int i = (Int)a.size() - 1; i >= 0; --i)
        {
          ULong128 v = (ULong128)a[i] | ((ULong128)carry << 64);
          a[i] = (DataT)(v >> 1);
          carry = (ULong)(v & 1);
        }
      }
      else
      {
        ULong carry = 0;
        for (Int i = (Int)a.size() - 1; i >= 0; --i)
        {
          ULong v = a[i] + carry * base;
          a[i] = (DataT)(v / 2);
          carry = v % 2;
        }
      }
      TrimZeros(a);
    }

    static void DivBy3InPlace(vector<DataT> &a, BaseT base)
    {
      if (a.empty()) return;
      if (base == Base2_32)
      {
        ULong carry = 0;
        for (Int i = (Int)a.size() - 1; i >= 0; --i)
        {
          ULong v = ((ULong)a[i]) | (carry << 32);
          a[i] = (DataT)(v / 3);
          carry = v % 3;
        }
      }
      else if (base == Base2_64)
      {
        // carry is 0, 1, or 2 (v % 3). v = (carry << 64) | a[i] fits ULong128.
        ULong carry = 0;
        for (Int i = (Int)a.size() - 1; i >= 0; --i)
        {
          ULong128 v = (ULong128)a[i] | ((ULong128)carry << 64);
          a[i] = (DataT)(v / 3);
          carry = (ULong)(v % 3);
        }
      }
      else
      {
        ULong carry = 0;
        for (Int i = (Int)a.size() - 1; i >= 0; --i)
        {
          ULong v = a[i] + carry * base;
          a[i] = (DataT)(v / 3);
          carry = v % 3;
        }
      }
      TrimZeros(a);
    }

    // Split x into x0, x1, x2 of size k limbs each (last may be shorter / empty).
    static void Split3(vector<DataT> const &x, SizeT k,
                      vector<DataT> &x0, vector<DataT> &x1, vector<DataT> &x2)
    {
      SizeT s = x.size();
      x0.assign(x.begin(), x.begin() + std::min(k, s));
      if (s > k)
        x1.assign(x.begin() + k, x.begin() + std::min((SizeT)(2 * k), s));
      else
        x1.assign(1, 0);
      if (s > 2 * k)
        x2.assign(x.begin() + 2 * k, x.end());
      else
        x2.assign(1, 0);
      TrimZeros(x0);
      TrimZeros(x1);
      TrimZeros(x2);
    }

  public:
    static vector<DataT> Multiply(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base)
    {
      if (IsZero(a) || IsZero(b))
        return vector<DataT>{0};
      if (a.size() == 1)
        return ClassicMultiplication::Multiply(b, a[0], base);
      if (b.size() == 1)
        return ClassicMultiplication::Multiply(a, b[0], base);

      SizeT n = (SizeT)std::max(a.size(), b.size());
      if (n < TOOM3_THRESHOLD)
        return KaratsubaMultiplication::Multiply(a, b, base);

      // Split each into 3 parts of k = ceil(n/3) limbs.
      SizeT k = (n + 2) / 3;

      vector<DataT> a0, a1, a2, b0, b1, b2;
      Split3(a, k, a0, a1, a2);
      Split3(b, k, b0, b1, b2);

      // --- Evaluate A and B at {0, 1, -1, 2, ∞} ---
      // p0 = a0 ;  pinf = a2
      // p1 = a0 + a1 + a2
      // pm1 = (a0 + a2) - a1   (signed)
      // p2 = a0 + 2 a1 + 4 a2
      vector<DataT> a0_plus_a2 = Add(a0, a2, base);
      vector<DataT> A1 = Add(a0_plus_a2, a1, base);
      Signed Am1 = SubSigned(a0_plus_a2, a1, base);
      vector<DataT> two_a1 = ClassicMultiplication::Multiply(a1, 2, base);
      vector<DataT> four_a2 = ClassicMultiplication::Multiply(a2, 4, base);
      vector<DataT> A2 = Add(a0, Add(two_a1, four_a2, base), base);

      vector<DataT> b0_plus_b2 = Add(b0, b2, base);
      vector<DataT> B1 = Add(b0_plus_b2, b1, base);
      Signed Bm1 = SubSigned(b0_plus_b2, b1, base);
      vector<DataT> two_b1 = ClassicMultiplication::Multiply(b1, 2, base);
      vector<DataT> four_b2 = ClassicMultiplication::Multiply(b2, 4, base);
      vector<DataT> B2 = Add(b0, Add(two_b1, four_b2, base), base);

      // --- 5 pointwise products (recursive). ---
      vector<DataT> r0   = Multiply(a0, b0, base);
      vector<DataT> r1   = Multiply(A1, B1, base);
      vector<DataT> rm1_mag = Multiply(Am1.mag, Bm1.mag, base);
      int rm1_sign = (IsZero(rm1_mag)) ? 1 : (Am1.sign * Bm1.sign);
      Signed rm1 = {std::move(rm1_mag), rm1_sign};
      vector<DataT> r2   = Multiply(A2, B2, base);
      vector<DataT> rinf = Multiply(a2, b2, base);

      // --- Interpolate. ---
      // Final polynomial coefficients c0..c4 satisfy P(x)=sum ci x^i.
      // c0 = r0, c4 = rinf already.
      // c2 = (r1 + rm1)/2 - c0 - c4
      // s  = (r1 - rm1)/2          → c1 + c3
      // t  = (r2 - c0 - 4 c2 - 16 c4)/2  → c1 + 4 c3
      // c3 = (t - s)/3
      // c1 = s - c3
      Signed R1 = MakePositive(r1);                  // r1
      Signed R0 = MakePositive(r0);                  // r0 (= c0)
      Signed Rinf = MakePositive(rinf);              // rinf (= c4)

      // c2 = (r1 + rm1)/2 - r0 - rinf
      Signed c2 = AddSigned(R1, rm1, base);
      HalveInPlace(c2.mag, base);
      c2 = SubSigned(c2, R0, base);
      c2 = SubSigned(c2, Rinf, base);

      // s = (r1 - rm1)/2
      Signed s = SubSigned(R1, rm1, base);
      HalveInPlace(s.mag, base);

      // t = (r2 - r0 - 4 c2 - 16 rinf)/2
      Signed R2 = MakePositive(r2);
      Signed four_c2 = c2;
      four_c2.mag = ClassicMultiplication::Multiply(four_c2.mag, 4, base);
      Signed sixteen_rinf = Rinf;
      sixteen_rinf.mag = ClassicMultiplication::Multiply(sixteen_rinf.mag, 16, base);
      Signed t = SubSigned(R2, R0, base);
      t = SubSigned(t, four_c2, base);
      t = SubSigned(t, sixteen_rinf, base);
      HalveInPlace(t.mag, base);

      // c3 = (t - s)/3
      Signed c3 = SubSigned(t, s, base);
      DivBy3InPlace(c3.mag, base);

      // c1 = s - c3
      Signed c1 = SubSigned(s, c3, base);

      // --- Combine: result = c0 + c1·B^k + c2·B^2k + c3·B^3k + c4·B^4k ---
      // All ci must be non-negative for a valid product.
      Signed result = R0;
      Signed term;

      term = {ShiftLeft(c1.mag, k), c1.sign};
      result = AddSigned(result, term, base);

      term = {ShiftLeft(c2.mag, 2 * k), c2.sign};
      result = AddSigned(result, term, base);

      term = {ShiftLeft(c3.mag, 3 * k), c3.sign};
      result = AddSigned(result, term, base);

      term = {ShiftLeft(Rinf.mag, 4 * k), Rinf.sign};
      result = AddSigned(result, term, base);

      TrimZeros(result.mag);
      return result.mag;
    }
  };
}

#endif
