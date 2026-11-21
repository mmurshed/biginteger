#ifndef BURNIKELZIEGLER_DIVISION
#define BURNIKELZIEGLER_DIVISION

#include <cstring>
#include <memory>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../../common/BitShifts.h"
#include "../../common/Comparator.h"
#include "../../common/LimbKernels.h"
#include "../../common/Util.h"
#include "../Addition.h"
#include "../Multiplication.h"
#include "../Shift.h"
#include "FastDivision.h"

namespace BigMath
{
// Basecase floor for the legacy (vector-based, 32-in-64 layout) 2n-by-n
// recursion below. Swept 2026-06-12: 128 beats the old 512 by 1.2-1.5×.
#ifndef BIGMATH_BZ_RECURSION_THRESHOLD
#define BIGMATH_BZ_RECURSION_THRESHOLD 128
#endif

// Basecase floor for the in-place 64-bit-limb recursion: a half of at most
// this many limbs is divided by the schoolbook kernel
// (FastDivision::DivideNormalized64) instead of being split again. Swept
// 2026-11-21 over 24-96 on divisors of 64-1536 limbs at ratios 1.3/2/4: flat
// to within 1-3%, 64 best or tied most often.
#ifndef BIGMATH_DC_DIV_THRESHOLD
#define BIGMATH_DC_DIV_THRESHOLD 64
#endif

  class BurnikelZieglerDivision
  {
  private:
    static const SizeT BZ_THRESHOLD = BIGMATH_BZ_RECURSION_THRESHOLD;

    static std::vector<DataT> NormalizeZero(std::vector<DataT> v)
    {
      TrimZerosToOne(v);
      return v;
    }

    static std::vector<DataT> Slice(std::vector<DataT> const &v, SizeT start, SizeT count)
    {
      if (start >= v.size() || count == 0)
        return std::vector<DataT>{0};
      SizeT end = std::min((SizeT)v.size(), (SizeT)(start + count));
      return NormalizeZero(std::vector<DataT>(v.begin() + start, v.begin() + end));
    }

    static std::vector<DataT> FromSpan(std::span<const DataT> v)
    {
      return NormalizeZero(std::vector<DataT>(v.begin(), v.end()));
    }

    static std::vector<DataT> MaxBlock(SizeT limbs, BaseT base)
    {
      if (limbs == 0)
        return std::vector<DataT>{0};
      // base - 1 underflows when base == Base2_64 sentinel (0); use LimbMask.
      const DataT limbMax = (base == Base2_64) ? (DataT)0xFFFFFFFFFFFFFFFFULL
                                                : (DataT)(base - 1);
      return std::vector<DataT>(limbs, limbMax);
    }

    static void PadTo(std::vector<DataT> &v, SizeT limbs)
    {
      if (IsZero(v))
        v.clear();
      while (v.size() < limbs)
        v.push_back(0);
    }

    static void Decrement(std::vector<DataT> &v, BaseT base)
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
    static std::vector<DataT> AddShifted(
        std::vector<DataT> const &high, SizeT shift,
        const DataT *lowData, SizeT lowSize,
        BaseT base)
    {
      SizeT outSize = std::max(high.size() + shift, (size_t)lowSize) + 1;
      std::vector<DataT> out(outSize, 0);

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

    static std::vector<DataT> CombineShifted(
        std::vector<DataT> const &high,
        SizeT shift,
        std::vector<DataT> const &low,
        BaseT base)
    {
      std::vector<DataT> out = AddShifted(high, shift, low.data(), (SizeT)low.size(), base);
      return NormalizeZero(std::move(out));
    }

    // Full Burnikel-Ziegler 2n-by-n division. For n = 2m, split the dividend
    // into four m-limb blocks and the divisor into two m-limb blocks, then use
    // two 3n-by-2n divisions.
    static std::pair<std::vector<DataT>, std::vector<DataT>> Divide2nByN(
        std::vector<DataT> const &a,
        std::vector<DataT> const &b,
        BaseT base,
        bool computeRemainder)
    {
      if (IsZero(b))
        throw std::invalid_argument("Division by zero");

      if (Compare(a, b) < 0)
        return {std::vector<DataT>{0}, computeRemainder ? a : std::vector<DataT>()};
      if (Compare(a, b) == 0)
        return {std::vector<DataT>{1}, computeRemainder ? std::vector<DataT>{0} : std::vector<DataT>()};

      SizeT n = (SizeT)b.size();
      if (n <= BZ_THRESHOLD || n % 2 != 0 || a.size() <= BZ_THRESHOLD)
        return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

      if (a.size() > 2 * n)
        return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

      SizeT m = n / 2;
      std::vector<DataT> a0 = Slice(a, 0, m);
      std::vector<DataT> a1 = Slice(a, m, m);
      std::vector<DataT> a2 = Slice(a, 2 * m, m);
      std::vector<DataT> a3 = Slice(a, 3 * m, m);
      std::vector<DataT> b0 = Slice(b, 0, m);
      std::vector<DataT> b1 = Slice(b, m, m);

      // Divide A = x2*B^(2m) + x1*B^m + x0 by
      // D = b1*B^m + b0, returning Q < B^m and R < D.
      auto divide3nBy2n = [&](std::vector<DataT> const &x2,
                              std::vector<DataT> const &x1,
                              std::vector<DataT> const &x0) {
        std::vector<DataT> top = CombineShifted(x2, m, x1, base);
        std::vector<DataT> q;
        std::vector<DataT> r;

        if (Compare(x2, b1) < 0)
        {
          auto qr = Divide2nByN(top, b1, base, true);
          q = std::move(qr.first);
          r = std::move(qr.second);
        }
        else
        {
          q = MaxBlock(m, base);
          std::vector<DataT> product = Multiply(q, b1, base);
          r = Subtract(top, product, base);
        }

        std::vector<DataT> d = Multiply(q, b0, base);
        r = CombineShifted(r, m, x0, base);
        std::vector<DataT> divisor = CombineShifted(b1, m, b0, base);

        while (Compare(r, d) < 0)
        {
          Decrement(q, base);
          r = Add(r, divisor, base);
        }

        r = Subtract(r, d, base);
        return std::pair<std::vector<DataT>, std::vector<DataT>>{
            NormalizeZero(std::move(q)),
            NormalizeZero(std::move(r))};
      };

      auto high = divide3nBy2n(a3, a2, a1);
      std::vector<DataT> r1Low = Slice(high.second, 0, m);
      std::vector<DataT> r1High = Slice(high.second, m, m);
      auto low = divide3nBy2n(r1High, r1Low, a0);

      PadTo(low.first, m);
      std::vector<DataT> q = CombineShifted(high.first, m, low.first, base);
      TrimZerosToOne(q);

      return {q, computeRemainder ? low.second : std::vector<DataT>()};
    }

    // Arbitrary-size wrapper: process the dividend in base B^n blocks and use
    // the balanced 2n-by-n primitive for each combined remainder/block step.
    static std::pair<std::vector<DataT>, std::vector<DataT>> DivideRecursive(
        std::vector<DataT> const &a,
        std::vector<DataT> const &b,
        BaseT base,
        bool computeRemainder)
    {
      if (a.size() <= BZ_THRESHOLD || b.size() <= BZ_THRESHOLD || b.size() % 2 != 0)
        return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

      SizeT n = (SizeT)b.size();
      SizeT blocks = (a.size() + n - 1) / n;
      std::vector<DataT> rem{0};
      std::vector<DataT> quotient(blocks * n, 0);

      for (SizeT block = blocks; block > 0; --block)
      {
        SizeT start = (block - 1) * n;
        std::vector<DataT> low = Slice(a, start, n);
        std::vector<DataT> combined = CombineShifted(rem, n, low, base);

        auto qr = Divide2nByN(combined, b, base, true);
        std::vector<DataT> qBlock = std::move(qr.first);
        rem = std::move(qr.second);

        if (qBlock.size() > n)
          return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

        std::memcpy(quotient.data() + start, qBlock.data(), qBlock.size() * sizeof(DataT));
      }

      TrimZerosToOne(quotient);

      return {quotient, computeRemainder ? NormalizeZero(std::move(rem)) : std::vector<DataT>()};
    }

    // ─── In-place divide-and-conquer division, 64-bit limbs ─────────────────
    //
    // Same recursion as Burnikel–Ziegler's 2n/n ↔ 3n/2n pair, written the way
    // GMP's mpn_dcpi1_div_qr does it: everything happens inside the (shifted)
    // dividend buffer, the divisor is split into low and high parts by
    // pointer arithmetic, and an odd length just splits into floor/ceil
    // halves — no slices, no per-level vectors, no padding to a power of two.
    //
    // To divide 2n limbs by n (quotient n limbs), with hi = ceil(n/2),
    // lo = floor(n/2):
    //   1. divide the top 2·hi limbs by the divisor's top hi limbs → q_hi;
    //   2. subtract q_hi · (the divisor's low lo limbs) from the running
    //      remainder; while that went negative, decrement q_hi and add the
    //      divisor back (at most a couple of times, since the divisor is
    //      normalized);
    //   3. the same for the next 2·lo limbs against the top lo limbs → q_lo,
    //      corrected with the divisor's low hi limbs.
    // Cost: D(n) = 2·D(n/2) + 2·M(n/2) → O(M(n) log n), bottoming out in the
    // schoolbook kernel at BIGMATH_DC_DIV_THRESHOLD.

    static const SizeT DC_THRESHOLD = BIGMATH_DC_DIV_THRESHOLD;

    // tp[0 .. la+lb) = a[0..la) · b[0..lb). Small products go straight to the
    // schoolbook kernel. Karatsuba writes directly into tp below the Toom/NTT
    // window. Larger shapes use the normal multiplication dispatcher.
    static void MulInto(DataT *tp, const DataT *a, SizeT la, const DataT *b, SizeT lb)
    {
      if (la == 0 || lb == 0)
        return;
      const SizeT shorter = std::min(la, lb);
      if (shorter <= 24 || detail::ClassicProductShape(shorter, std::max(la, lb)))
      {
        LimbKernel::MulBasecase(tp, a, la, b, lb);
        return;
      }
      // Below either transform/Toom entry no pool-dependent decision is needed.
      // Respect overridden classic and skew thresholds on tuned builds.
      if ((ULong)la + lb < std::min(NTT_MULTIPLICATION_THRESHOLD, TOOM3_MULTIPLICATION_THRESHOLD))
      {
        static thread_local std::vector<DataT> scratch;
        const size_t count = KaratsubaMultiplication::ScratchLimbs(la, lb);
        if (scratch.size() < count) scratch.resize(count);
        KaratsubaMultiplication::MultiplyInto(a, la, b, lb, tp, scratch.data(), Base2_64);
        return;
      }
      static thread_local std::vector<DataT> av, bv;
      av.assign(a, a + la);
      bv.assign(b, b + lb);
      std::vector<DataT> product = Multiply(av, bv, Base2_64);
      const SizeT total = la + lb;
      const SizeT have = std::min((SizeT)product.size(), total);
      std::memcpy(tp, product.data(), (size_t)have * sizeof(DataT));
      if (have < total)
        std::memset(tp + have, 0, (size_t)(total - have) * sizeof(DataT));
    }

    // Schoolbook 2k-by-k division: np[0..2k) / dp[0..k), k >= 2, dp normalized.
    // Quotient low k limbs in qp, remainder in np[0..k); returns the quotient's
    // high bit (set when the top k limbs of np were >= dp).
    static DataT SbDivQr(DataT *qp, DataT *np, const DataT *dp, SizeT k, ULong vinv)
    {
      DataT qh = 0;
      if (LimbCompare(np + k, dp, k) >= 0)
      {
        LimbKernel::SubN(np + k, np + k, dp, k);
        qh = 1;
      }
      FastDivision::DivideNormalized64(qp, np, 2 * k - 1, dp, k, vinv);
      return qh;
    }

    static int LimbCompare(const DataT *a, const DataT *b, SizeT n)
    {
      for (SizeT i = n; i-- > 0;)
        if (a[i] != b[i])
          return a[i] > b[i] ? 1 : -1;
      return 0;
    }

    // np[0..2n) / dp[0..n): quotient low n limbs in qp, remainder in np[0..n),
    // returns the quotient's high bit. tp is n limbs of scratch. n >= 4.
    static DataT DivQrN(DataT *qp, DataT *np, const DataT *dp, SizeT n, ULong vinv, DataT *tp)
    {
      const SizeT lo = n >> 1;  // floor(n/2)
      const SizeT hi = n - lo;  // ceil(n/2)

      DataT qh = (hi <= DC_THRESHOLD || hi < 4)
                     ? SbDivQr(qp + lo, np + 2 * lo, dp + lo, hi, vinv)
                     : DivQrN(qp + lo, np + 2 * lo, dp + lo, hi, vinv, tp);

      MulInto(tp, qp + lo, hi, dp, lo);
      DataT cy = LimbKernel::SubN(np + lo, np + lo, tp, n);
      if (qh != 0)
        cy += LimbKernel::SubN(np + n, np + n, dp, lo);
      while (cy != 0)
      {
        qh -= LimbKernel::SubBorrow(qp + lo, qp + lo, hi, 1);
        cy -= LimbKernel::AddN(np + lo, np + lo, dp, n);
      }

      DataT ql = (lo <= DC_THRESHOLD || lo < 4)
                     ? SbDivQr(qp, np + hi, dp + hi, lo, vinv)
                     : DivQrN(qp, np + hi, dp + hi, lo, vinv, tp);

      MulInto(tp, dp, hi, qp, lo);
      cy = LimbKernel::SubN(np, np, tp, n);
      if (ql != 0)
        cy += LimbKernel::SubN(np + lo, np + lo, dp, hi);
      while (cy != 0)
      {
        LimbKernel::SubBorrow(qp, qp, lo, 1);
        cy -= LimbKernel::AddN(np, np, dp, n);
      }
      return qh;
    }

    // One block: wp[0 .. dn+k) / dp[0..dn) with 1 <= k <= dn and the top dn
    // limbs of wp already < dp, so the quotient is exactly k limbs (qp[0..k))
    // and the remainder ends up in wp[0..dn). tp is dn limbs of scratch.
    static void DivBlock(DataT *qp, DataT *wp, const DataT *dp, SizeT dn, SizeT k, ULong vinv, DataT *tp)
    {
      if (k <= DC_THRESHOLD || k < 4 || dn <= DC_THRESHOLD)
      {
        // Short quotient or short divisor: schoolbook against the whole divisor.
        FastDivision::DivideNormalized64(qp, wp, dn + k - 1, dp, dn, vinv);
        return;
      }
      if (k == dn)
      {
        DivQrN(qp, wp, dp, dn, vinv, tp); // high bit is 0 by the precondition
        return;
      }
      // k < dn: divide the top 2k limbs by the divisor's top k limbs, then
      // correct with the remaining dn - k low limbs of the divisor.
      const SizeT low = dn - k;
      DataT qh = DivQrN(qp, wp + low, dp + low, k, vinv, tp);
      MulInto(tp, qp, k, dp, low);
      DataT cy = LimbKernel::SubN(wp, wp, tp, dn);
      if (qh != 0)
        cy += LimbKernel::SubN(wp + k, wp + k, dp, low);
      while (cy != 0)
      {
        qh -= LimbKernel::SubBorrow(qp, qp, k, 1);
        cy -= LimbKernel::AddN(wp, wp, dp, dn);
      }
    }

    // a > b > 0, b has at least 2 limbs after trimming.
    static std::pair<std::vector<DataT>, std::vector<DataT>> DivideAndRemainder64(
        std::span<const DataT> a,
        std::span<const DataT> b,
        bool computeRemainder)
    {
      SizeT na = (SizeT)a.size();
      while (na > 1 && a[na - 1] == 0) --na;
      SizeT dn = (SizeT)b.size();
      while (dn > 1 && b[dn - 1] == 0) --dn;

      // Workspace: shifted dividend with one extra top limb, shifted divisor,
      // and dn limbs of product scratch.
      const size_t need = (size_t)na + 1 + 2 * (size_t)dn;
      std::unique_ptr<DataT[]> work(new DataT[need]);
      DataT *np = work.get();
      DataT *dp = np + na + 1;
      DataT *tp = dp + dn;

      const unsigned shift = (unsigned)__builtin_clzll((unsigned long long)b[dn - 1]);
      if (shift != 0)
      {
        np[na] = LimbKernel::Lshift(np, a.data(), na, shift);
        LimbKernel::Lshift(dp, b.data(), dn, shift);
      }
      else
      {
        std::memcpy(np, a.data(), (size_t)na * sizeof(DataT));
        np[na] = 0;
        std::memcpy(dp, b.data(), (size_t)dn * sizeof(DataT));
      }
      const ULong vinv = FastDivision::Reciprocal3by2((ULong)dp[dn - 1], (ULong)dp[dn - 2]);

      // np now has nn = na + 1 limbs whose top dn limbs are < dp (the extra
      // limb holds fewer than 64 bits), so the quotient is exactly
      // qn = nn - dn limbs. Peel it off from the top: a first block of
      // qn mod dn limbs (or dn), then full dn-limb blocks, each leaving its
      // remainder as the top of the next window.
      const SizeT nn = na + 1;
      SizeT remaining = nn - dn;
      std::vector<DataT> q(remaining);
      SizeT k = remaining % dn;
      if (k == 0)
        k = dn;
      while (remaining > 0)
      {
        const SizeT pos = remaining - k;
        DivBlock(q.data() + pos, np + pos, dp, dn, k, vinv, tp);
        remaining = pos;
        k = dn;
      }
      TrimZerosToOne(q);

      std::vector<DataT> r;
      if (computeRemainder)
      {
        r.resize(dn);
        if (shift != 0)
          LimbKernel::Rshift(r.data(), np, dn, shift);
        else
          std::memcpy(r.data(), np, (size_t)dn * sizeof(DataT));
        TrimZerosToOne(r);
      }
      return {std::move(q), std::move(r)};
    }

  public:
    static std::pair<std::vector<DataT>, std::vector<DataT>> DivideAndRemainder(
        std::span<const DataT> a,
        std::span<const DataT> b,
        BaseT base,
        bool computeRemainder = true)
    {
      if (IsZero(b))
        throw std::invalid_argument("Division by zero");

      if (IsZero(a))
        return {std::vector<DataT>{0}, computeRemainder ? std::vector<DataT>{0} : std::vector<DataT>()};

      Int cmp = Compare(a, b);
      if (cmp < 0)
        return {std::vector<DataT>{0}, computeRemainder ? std::vector<DataT>(a.begin(), a.end()) : std::vector<DataT>()};
      if (cmp == 0)
        return {std::vector<DataT>{1}, computeRemainder ? std::vector<DataT>{0} : std::vector<DataT>()};

      if (base == Base2_64)
      {
        SizeT dn = (SizeT)b.size();
        while (dn > 1 && b[dn - 1] == 0) --dn;
        if (dn < 2)
          return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);
        return DivideAndRemainder64(a, b, computeRemainder);
      }

      if (a.size() <= BZ_THRESHOLD || b.size() <= BZ_THRESHOLD)
        return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

      // Knuth-normalize the divisor (top bit of b.back() set). Required for BZ's
      // 3n/2n correction-loop bound of 2 iterations; without it the loop runs O(B)
      // times when the recursive q-estimate is wildly off (PR #30 LIMB_64 sparse-input
      // hang was here — divisor.back() = 1 forced ~2^64 correction iters).
      const int limbBits = LimbBitsFor(base);
      Int shift = NormalizationShiftBits(b[b.size() - 1], limbBits);

      std::vector<DataT> aNorm = ShiftLeftBits(a, shift, limbBits);
      std::vector<DataT> bNorm = ShiftLeftBits(b, shift, limbBits);

      // BZ recursion needs the divisor size divisible by 2 at every level of
      // the 2n-by-n split until the basecase; any odd size on the way down
      // used to fall back to FastDivision wholesale — O(n·Δ) at full size, the
      // root cause of the "2^k+1 family" pathology (97.7× at 2^18+1) and a
      // 1.3-10× loss on ~half of real divisor sizes (odd limb counts).
      // Instead pad both operands by k zero limbs at the bottom so the divisor
      // size becomes ceil(n/2^d)·2^d with d = halvings to reach the basecase:
      // a·B^k = q·(b·B^k) + r·B^k, so the quotient is unchanged and the
      // remainder's bottom k limbs are exactly zero. k < 2^d ≈ n/BZ_THRESHOLD,
      // so the padding overhead is bounded by ~1/BZ_THRESHOLD (≈0.2%).
      SizeT n = (SizeT)bNorm.size();
      SizeT levels = 0;
      SizeT basecase = n;
      while (basecase > BZ_THRESHOLD)
      {
        basecase = (basecase + 1) / 2;
        ++levels;
      }
      SizeT padded = basecase << levels;
      SizeT k = padded - n;
      if (k > 0)
      {
        aNorm.insert(aNorm.begin(), k, 0);
        bNorm.insert(bNorm.begin(), k, 0);
      }

      auto qr = DivideRecursive(aNorm, bNorm, base, computeRemainder);
      if (computeRemainder && k > 0)
      {
        if (qr.second.size() > k)
          qr.second.erase(qr.second.begin(), qr.second.begin() + k);
        else
          qr.second.assign(1, 0);
      }
      if (computeRemainder && shift > 0)
        qr.second = ShiftRightBits(qr.second, shift, limbBits);
      return qr;
    }

    static std::vector<DataT> Divide(std::span<const DataT> a, std::span<const DataT> b, BaseT base)
    {
      return DivideAndRemainder(a, b, base, false).first;
    }

    // Vector overloads — backward compat.
    static std::pair<std::vector<DataT>, std::vector<DataT>> DivideAndRemainder(
        std::vector<DataT> const &a,
        std::vector<DataT> const &b,
        BaseT base,
        bool computeRemainder = true)
    {
      return DivideAndRemainder(std::span<const DataT>(a), std::span<const DataT>(b), base, computeRemainder);
    }

    static std::vector<DataT> Divide(std::vector<DataT> const &a, std::vector<DataT> const &b, BaseT base)
    {
      return DivideAndRemainder(a, b, base, false).first;
    }
  };
}

#endif
