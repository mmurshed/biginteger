#ifndef BURNIKELZIEGLER_DIVISION
#define BURNIKELZIEGLER_DIVISION

#include <cstring>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>
using namespace std;

#include "../../common/Comparator.h"
#include "../../common/Util.h"
#include "../Addition.h"
#include "../Multiplication.h"
#include "../Shift.h"
#include "FastDivision.h"

namespace BigMath
{
#ifndef BIGMATH_BZ_RECURSION_THRESHOLD
#define BIGMATH_BZ_RECURSION_THRESHOLD 512
#endif

  class BurnikelZieglerDivision
  {
  private:
    static const SizeT BZ_THRESHOLD = BIGMATH_BZ_RECURSION_THRESHOLD;

    static vector<DataT> NormalizeZero(vector<DataT> v)
    {
      TrimZerosToOne(v);
      return v;
    }

    static vector<DataT> Slice(vector<DataT> const &v, SizeT start, SizeT count)
    {
      if (start >= v.size() || count == 0)
        return vector<DataT>{0};
      SizeT end = min((SizeT)v.size(), (SizeT)(start + count));
      return NormalizeZero(vector<DataT>(v.begin() + start, v.begin() + end));
    }

    static vector<DataT> FromSpan(span<const DataT> v)
    {
      return NormalizeZero(vector<DataT>(v.begin(), v.end()));
    }

    static vector<DataT> MaxBlock(SizeT limbs, BaseT base)
    {
      if (limbs == 0)
        return vector<DataT>{0};
      // base - 1 underflows when base == Base2_64 sentinel (0); use LimbMask.
      const DataT limbMax = (base == Base2_64) ? (DataT)0xFFFFFFFFFFFFFFFFULL
                                                : (DataT)(base - 1);
      return vector<DataT>(limbs, limbMax);
    }

    static void PadTo(vector<DataT> &v, SizeT limbs)
    {
      if (IsZero(v))
        v.clear();
      while (v.size() < limbs)
        v.push_back(0);
    }

    // Bit-shift left by [1..LimbBits] bits. Caller ensures shift ∈ [0,LimbBits].
    static vector<DataT> ShiftLeftBits(vector<DataT> const &v, Int shift)
    {
      if (shift == 0 || IsZero(v))
        return v;
      vector<DataT> out(v.size() + 1, 0);
      for (SizeT i = 0; i < v.size(); ++i)
      {
#if BIGMATH_LIMB_64
        ULong128 cur = (ULong128)v[i] << shift;
        out[i] |= (DataT)(cur & 0xFFFFFFFFFFFFFFFFULL);
        out[i + 1] |= (DataT)(cur >> 64);
#else
        ULong cur = (ULong)v[i] << shift;
        out[i] |= (DataT)(cur & 0xFFFFFFFFULL);
        out[i + 1] |= (DataT)(cur >> 32);
#endif
      }
      TrimZerosToOne(out);
      return out;
    }

    // Bit-shift right by [0..LimbBits] bits.
    static vector<DataT> ShiftRightBits(vector<DataT> const &v, Int shift)
    {
      if (shift == 0 || IsZero(v))
        return v;
      vector<DataT> out(v.size(), 0);
      for (Int i = 0; i < (Int)v.size(); ++i)
      {
#if BIGMATH_LIMB_64
        ULong128 cur = (ULong128)v[i];
        if (i + 1 < (Int)v.size())
          cur |= ((ULong128)v[i + 1]) << 64;
        out[i] = (DataT)((cur >> shift) & 0xFFFFFFFFFFFFFFFFULL);
#else
        ULong cur = v[i];
        if (i + 1 < (Int)v.size())
          cur |= ((ULong)v[i + 1]) << 32;
        out[i] = (DataT)((cur >> shift) & 0xFFFFFFFFULL);
#endif
      }
      TrimZerosToOne(out);
      return out;
    }

    // Number of bits to shift left so that b's limb count becomes b.size()+1.
    // Pushes b's MSB into a new high limb. Returns shift ∈ [1,LimbBits].
    static Int BitsToBumpLimbCount(vector<DataT> const &b)
    {
      DataT top = b.back();
#if BIGMATH_LIMB_64
      Int bits_top = 64 - __builtin_clzll((unsigned long long)top);
      return 65 - bits_top; // 1..64
#else
      Int bits_top = 32 - __builtin_clz((unsigned int)top);
      return 33 - bits_top; // 1..32
#endif
    }

    static void Decrement(vector<DataT> &v, BaseT base)
    {
      // limbMax is base-1 for power-of-two bases. The BaseT sentinel for
      // Base2_64 is 0, so (base - 1) would underflow — use LimbMask.
      const DataT limbMax = (base == Base2_64) ? (DataT)0xFFFFFFFFFFFFFFFFULL
                                                : (DataT)(base - 1);
      SizeT i = 0;
      while (i < v.size())
      {
        if (v[i] > 0)
        {
          --v[i];
          break;
        }
        v[i] = limbMax;
        ++i;
      }
      TrimZerosToOne(v);
    }

    // out = (high << shift_limbs) + low. One allocation.
    static vector<DataT> AddShifted(
        vector<DataT> const &high, SizeT shift,
        const DataT *lowData, SizeT lowSize,
        BaseT base)
    {
      SizeT outSize = max(high.size() + shift, (size_t)lowSize) + 1;
      vector<DataT> out(outSize, 0);

      std::memcpy(out.data(), lowData, lowSize * sizeof(DataT));

      if (base == Base2_64)
      {
        ULong128 carry = 0;
        SizeT i = 0;
        for (; i < high.size(); ++i)
        {
          ULong128 sum = (ULong128)out[shift + i] + high[i] + carry;
          out[shift + i] = (DataT)(sum & 0xFFFFFFFFFFFFFFFFULL);
          carry = sum >> 64;
        }
        SizeT idx = shift + i;
        while (carry)
        {
          if (idx >= out.size())
            out.push_back(0);
          ULong128 sum = (ULong128)out[idx] + carry;
          out[idx] = (DataT)(sum & 0xFFFFFFFFFFFFFFFFULL);
          carry = sum >> 64;
          ++idx;
        }
        return out;
      }

      ULong carry = 0;
      SizeT i = 0;
      for (; i < high.size(); ++i)
      {
        ULong sum = (ULong)out[shift + i] + high[i] + carry;
        if (base == Base2_32)
        {
          out[shift + i] = (DataT)(sum & 0xFFFFFFFFULL);
          carry = sum >> 32;
        }
        else
        {
          out[shift + i] = (DataT)(sum % base);
          carry = sum / base;
        }
      }
      SizeT idx = shift + i;
      while (carry)
      {
        if (idx >= out.size())
          out.push_back(0);
        ULong sum = (ULong)out[idx] + carry;
        if (base == Base2_32)
        {
          out[idx] = (DataT)(sum & 0xFFFFFFFFULL);
          carry = sum >> 32;
        }
        else
        {
          out[idx] = (DataT)(sum % base);
          carry = sum / base;
        }
        ++idx;
      }
      return out;
    }

    static vector<DataT> CombineShifted(
        vector<DataT> const &high,
        SizeT shift,
        vector<DataT> const &low,
        BaseT base)
    {
      vector<DataT> out = AddShifted(high, shift, low.data(), (SizeT)low.size(), base);
      return NormalizeZero(std::move(out));
    }

    // Full Burnikel-Ziegler 2n-by-n division. For n = 2m, split the dividend
    // into four m-limb blocks and the divisor into two m-limb blocks, then use
    // two 3n-by-2n divisions.
    static pair<vector<DataT>, vector<DataT>> Divide2nByN(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base,
        bool computeRemainder)
    {
      if (IsZero(b))
        throw invalid_argument("Division by zero");

      if (Compare(a, b) < 0)
        return {vector<DataT>{0}, computeRemainder ? a : vector<DataT>()};
      if (Compare(a, b) == 0)
        return {vector<DataT>{1}, computeRemainder ? vector<DataT>{0} : vector<DataT>()};

      SizeT n = (SizeT)b.size();
      if (n <= BZ_THRESHOLD || n % 2 != 0 || a.size() <= BZ_THRESHOLD)
        return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

      if (a.size() > 2 * n)
        return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

      SizeT m = n / 2;
      vector<DataT> a0 = Slice(a, 0, m);
      vector<DataT> a1 = Slice(a, m, m);
      vector<DataT> a2 = Slice(a, 2 * m, m);
      vector<DataT> a3 = Slice(a, 3 * m, m);
      vector<DataT> b0 = Slice(b, 0, m);
      vector<DataT> b1 = Slice(b, m, m);

      // Divide A = x2*B^(2m) + x1*B^m + x0 by
      // D = b1*B^m + b0, returning Q < B^m and R < D.
      auto divide3nBy2n = [&](vector<DataT> const &x2,
                              vector<DataT> const &x1,
                              vector<DataT> const &x0) {
        vector<DataT> top = CombineShifted(x2, m, x1, base);
        vector<DataT> q;
        vector<DataT> r;

        if (Compare(x2, b1) < 0)
        {
          auto qr = Divide2nByN(top, b1, base, true);
          q = std::move(qr.first);
          r = std::move(qr.second);
        }
        else
        {
          q = MaxBlock(m, base);
          vector<DataT> product = Multiply(q, b1, base);
          r = Subtract(top, product, base);
        }

        vector<DataT> d = Multiply(q, b0, base);
        r = CombineShifted(r, m, x0, base);
        vector<DataT> divisor = CombineShifted(b1, m, b0, base);

        while (Compare(r, d) < 0)
        {
          Decrement(q, base);
          r = Add(r, divisor, base);
        }

        r = Subtract(r, d, base);
        return pair<vector<DataT>, vector<DataT>>{
            NormalizeZero(std::move(q)),
            NormalizeZero(std::move(r))};
      };

      auto high = divide3nBy2n(a3, a2, a1);
      vector<DataT> r1Low = Slice(high.second, 0, m);
      vector<DataT> r1High = Slice(high.second, m, m);
      auto low = divide3nBy2n(r1High, r1Low, a0);

      PadTo(low.first, m);
      vector<DataT> q = CombineShifted(high.first, m, low.first, base);
      TrimZerosToOne(q);

      return {q, computeRemainder ? low.second : vector<DataT>()};
    }

    // Arbitrary-size wrapper: process the dividend in base B^n blocks and use
    // the balanced 2n-by-n primitive for each combined remainder/block step.
    static pair<vector<DataT>, vector<DataT>> DivideRecursive(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base,
        bool computeRemainder)
    {
      if (a.size() <= BZ_THRESHOLD || b.size() <= BZ_THRESHOLD || b.size() % 2 != 0)
        return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

      SizeT n = (SizeT)b.size();
      SizeT blocks = (a.size() + n - 1) / n;
      vector<DataT> rem{0};
      vector<DataT> quotient(blocks * n, 0);

      for (SizeT block = blocks; block > 0; --block)
      {
        SizeT start = (block - 1) * n;
        vector<DataT> low = Slice(a, start, n);
        vector<DataT> combined = CombineShifted(rem, n, low, base);

        auto qr = Divide2nByN(combined, b, base, true);
        vector<DataT> qBlock = std::move(qr.first);
        rem = std::move(qr.second);

        if (qBlock.size() > n)
          return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

        std::memcpy(quotient.data() + start, qBlock.data(), qBlock.size() * sizeof(DataT));
      }

      TrimZerosToOne(quotient);

      return {quotient, computeRemainder ? NormalizeZero(std::move(rem)) : vector<DataT>()};
    }

  public:
    static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
        span<const DataT> a,
        span<const DataT> b,
        BaseT base,
        bool computeRemainder = true)
    {
      if (IsZero(b))
        throw invalid_argument("Division by zero");

      if (IsZero(a))
        return {vector<DataT>{0}, computeRemainder ? vector<DataT>{0} : vector<DataT>()};

      Int cmp = Compare(a, b);
      if (cmp < 0)
        return {vector<DataT>{0}, computeRemainder ? vector<DataT>(a.begin(), a.end()) : vector<DataT>()};
      if (cmp == 0)
        return {vector<DataT>{1}, computeRemainder ? vector<DataT>{0} : vector<DataT>()};

      if (a.size() <= BZ_THRESHOLD || b.size() <= BZ_THRESHOLD)
        return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

      // Knuth-normalize the divisor (top bit of b.back() set). Required for BZ's
      // 3n/2n correction-loop bound of 2 iterations; without it the loop runs O(B)
      // times when the recursive q-estimate is wildly off (PR #30 LIMB_64 sparse-input
      // hang was here — divisor.back() = 1 forced ~2^64 correction iters).
      DataT b_top = b[b.size() - 1];
      Int shift = 0;
#if BIGMATH_LIMB_64
      const DataT topBitMask = 0x8000000000000000ULL;
#else
      const DataT topBitMask = 0x80000000U;
#endif
      while ((b_top & topBitMask) == 0)
      {
        b_top <<= 1;
        ++shift;
      }

      vector<DataT> aNorm = (shift > 0)
                                ? ShiftLeftBits(vector<DataT>(a.begin(), a.end()), shift)
                                : vector<DataT>(a.begin(), a.end());
      vector<DataT> bNorm = (shift > 0)
                                ? ShiftLeftBits(vector<DataT>(b.begin(), b.end()), shift)
                                : vector<DataT>(b.begin(), b.end());

      // BZ recursion needs even divisor size for the 2n-by-n split. If post-normalize
      // size is odd, fall to FastDivision — which has its own internal normalization
      // and isn't sensitive to BZ's correction-loop bound.
      if (bNorm.size() % 2 != 0)
      {
        auto qr = FastDivision::DivideAndRemainder(aNorm, bNorm, base, computeRemainder);
        if (computeRemainder && shift > 0)
          qr.second = ShiftRightBits(qr.second, shift);
        return qr;
      }

      auto qr = DivideRecursive(aNorm, bNorm, base, computeRemainder);
      if (computeRemainder && shift > 0)
        qr.second = ShiftRightBits(qr.second, shift);
      return qr;
    }

    static vector<DataT> Divide(span<const DataT> a, span<const DataT> b, BaseT base)
    {
      return DivideAndRemainder(a, b, base, false).first;
    }

    // Vector overloads — backward compat.
    static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base,
        bool computeRemainder = true)
    {
      return DivideAndRemainder(span<const DataT>(a), span<const DataT>(b), base, computeRemainder);
    }

    static vector<DataT> Divide(vector<DataT> const &a, vector<DataT> const &b, BaseT base)
    {
      return DivideAndRemainder(a, b, base, false).first;
    }
  };
}

#endif
