/**
 * BigInteger Class
 * Toom-5 multiplication: 5-way split, evaluation at {0, ±1, ±2, ±3, 4, ∞}.
 *
 * This implementation is intentionally kept out of dispatch until benchmarks
 * show a stable win over Karatsuba before the NTT crossover.
 */

#ifndef TOOM5_MULTIPLICATION
#define TOOM5_MULTIPLICATION

#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>
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
  class Toom5Multiplication
  {
  private:
#ifndef BIGMATH_TOOM5_THRESHOLD
#define BIGMATH_TOOM5_THRESHOLD 512
#endif
    static const SizeT TOOM5_THRESHOLD = BIGMATH_TOOM5_THRESHOLD;

    struct Signed
    {
      vector<DataT> mag;
      int sign;
    };

    struct Fraction
    {
      long long num;
      long long den;
    };

    static long long GcdAbs(long long a, long long b)
    {
      a = a < 0 ? -a : a;
      b = b < 0 ? -b : b;
      return std::gcd(a, b);
    }

    static Fraction Normalize(Fraction f)
    {
      if (f.den == 0)
        throw std::runtime_error("invalid Toom-5 interpolation denominator");
      if (f.den < 0)
      {
        f.num = -f.num;
        f.den = -f.den;
      }
      long long g = GcdAbs(f.num, f.den);
      if (g != 0)
      {
        f.num /= g;
        f.den /= g;
      }
      return f;
    }

    static Signed MakePositive(vector<DataT> v)
    {
      TrimZeros(v);
      return {std::move(v), 1};
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

    static Signed NegSigned(Signed v)
    {
      if (!IsZero(v.mag))
        v.sign = -v.sign;
      return v;
    }

    static Signed SubSigned(Signed const &a, Signed const &b, BaseT base)
    {
      return AddSigned(a, NegSigned(b), base);
    }

    static Signed MulSmall(Signed v, ULong m, BaseT base)
    {
      if (m == 0 || IsZero(v.mag))
        return {vector<DataT>{0}, 1};
      v.mag = ClassicMultiplication::Multiply(v.mag, m, base);
      return v;
    }

    static bool DivSmallInPlace(vector<DataT> &a, ULong d, BaseT base)
    {
      if (d == 1 || a.empty())
        return true;
      ULong carry = 0;
      if (base == Base2_32)
      {
        for (Int i = (Int)a.size() - 1; i >= 0; --i)
        {
          ULong v = ((ULong)a[i]) | (carry << 32);
          a[i] = (DataT)(v / d);
          carry = v % d;
        }
      }
      else
      {
        for (Int i = (Int)a.size() - 1; i >= 0; --i)
        {
          ULong v = a[i] + carry * base;
          a[i] = (DataT)(v / d);
          carry = v % d;
        }
      }
      TrimZeros(a);
      return carry == 0;
    }

    static Signed DivSmall(Signed v, ULong d, BaseT base)
    {
      if (!DivSmallInPlace(v.mag, d, base))
        throw std::runtime_error("non-exact Toom-5 interpolation division");
      if (IsZero(v.mag))
        v.sign = 1;
      return v;
    }

    static void Split5(vector<DataT> const &x, SizeT k, array<vector<DataT>, 5> &parts)
    {
      SizeT s = x.size();
      for (SizeT p = 0; p < 5; ++p)
      {
        SizeT start = p * k;
        if (start < s)
          parts[p].assign(x.begin() + start, x.begin() + std::min(s, start + k));
        else
          parts[p].assign(1, 0);
        TrimZeros(parts[p]);
      }
    }

    static Signed Evaluate(array<vector<DataT>, 5> const &parts, int x, BaseT base)
    {
      Signed acc = MakePositive(parts[4]);
      for (Int i = 3; i >= 0; --i)
      {
        if (x < 0)
          acc = NegSigned(MulSmall(acc, (ULong)(-x), base));
        else
          acc = MulSmall(acc, (ULong)x, base);
        acc = AddSigned(acc, MakePositive(parts[(SizeT)i]), base);
      }
      return acc;
    }

    static const array<array<Fraction, 7>, 7> &InverseInterpolation()
    {
      static const array<array<Fraction, 7>, 7> inv = {{
        {{{1, 1}, {-3, 5}, {-3, 10}, {1, 10}, {1, 15}, {-1, 105}, {-1, 140}}},
        {{{3, 4}, {3, 4}, {-3, 40}, {-3, 40}, {1, 180}, {1, 180}, {0, 1}}},
        {{{-11, 18}, {1, 15}, {89, 240}, {-71, 720}, {-4, 45}, {1, 90}, {7, 720}}},
        {{{-13, 48}, {-13, 48}, {1, 12}, {1, 12}, {-1, 144}, {-1, 144}, {0, 1}}},
        {{{17, 144}, {3, 80}, {-3, 40}, {-1, 360}, {17, 720}, {-1, 720}, {-1, 360}}},
        {{{1, 48}, {1, 48}, {-1, 120}, {-1, 120}, {1, 720}, {1, 720}, {0, 1}}},
        {{{-1, 144}, {-1, 240}, {1, 240}, {1, 720}, {-1, 720}, {-1, 5040}, {1, 5040}}}
      }};
      return inv;
    }

    static Signed LinearCombination(
        array<Signed, 7> const &ys,
        array<Fraction, 7> const &weights,
        BaseT base)
    {
      ULong common = 1;
      for (Fraction const &w : weights)
        common = std::lcm(common, (ULong)w.den);

      Signed total = {vector<DataT>{0}, 1};
      for (SizeT i = 0; i < 7; ++i)
      {
        if (weights[i].num == 0)
          continue;
        ULong scale = common / (ULong)weights[i].den * (ULong)(weights[i].num < 0 ? -weights[i].num : weights[i].num);
        Signed term = MulSmall(ys[i], scale, base);
        if (weights[i].num < 0)
          term = NegSigned(term);
        total = AddSigned(total, term, base);
      }

      if (!DivSmallInPlace(total.mag, common, base))
        throw std::runtime_error("non-exact Toom-5 interpolation division");
      if (IsZero(total.mag))
        total.sign = 1;
      return total;
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
      if (n < TOOM5_THRESHOLD)
        return KaratsubaMultiplication::Multiply(a, b, base);

      SizeT k = (n + 4) / 5;
      array<vector<DataT>, 5> ap, bp;
      Split5(a, k, ap);
      Split5(b, k, bp);

      vector<DataT> c0 = Multiply(ap[0], bp[0], base);
      vector<DataT> c8 = Multiply(ap[4], bp[4], base);

      constexpr int X[7] = {1, -1, 2, -2, 3, -3, 4};
      array<Signed, 7> y{};
      Signed C0 = MakePositive(c0);
      Signed C8 = MakePositive(c8);
      for (SizeT i = 0; i < 7; ++i)
      {
        Signed ea = Evaluate(ap, X[i], base);
        Signed eb = Evaluate(bp, X[i], base);
        vector<DataT> prod = Multiply(ea.mag, eb.mag, base);
        int prodSign = IsZero(prod) ? 1 : ea.sign * eb.sign;
        Signed value = {std::move(prod), prodSign};

        Signed x8c8 = MulSmall(C8, (ULong)std::llabs((long long)X[i] * X[i] * X[i] * X[i] * X[i] * X[i] * X[i] * X[i]), base);
        y[i] = SubSigned(SubSigned(value, C0, base), x8c8, base);
      }

      array<Signed, 9> coeff{};
      coeff[0] = C0;
      coeff[8] = C8;
      const array<array<Fraction, 7>, 7> &inv = InverseInterpolation();
      for (SizeT i = 0; i < 7; ++i)
      {
        try
        {
          coeff[i + 1] = LinearCombination(y, inv[i], base);
        }
        catch (std::runtime_error const &)
        {
          throw std::runtime_error("non-exact Toom-5 interpolation");
        }
      }

      Signed result = coeff[0];
      for (SizeT i = 1; i < 9; ++i)
      {
        Signed term = {ShiftLeft(coeff[i].mag, i * k), coeff[i].sign};
        result = AddSigned(result, term, base);
      }

      TrimZeros(result.mag);
      return result.mag;
    }
  };
}

#endif
