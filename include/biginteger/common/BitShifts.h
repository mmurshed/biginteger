/**
 * BigMath: Limb-vector bit-shift helpers for power-of-two bases.
 *
 * Shared by FastDivision, BurnikelZieglerDivision and NewtonDivision, which
 * previously carried three near-identical private copies (two of them
 * compile-time-width, silently wrong for the non-default base). The limb
 * width is a runtime parameter derived from the base, so the same code
 * serves Base2_32 and Base2_64 builds and cross-base calls.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGMATH_BIT_SHIFTS
#define BIGMATH_BIT_SHIFTS

#include <span>
#include <vector>

#include "Constants.h"
#include "Util.h"

namespace BigMath
{
  // Limb width in bits for a power-of-two base. Callers guarantee
  // base ∈ {Base2_32, Base2_64}.
  inline int LimbBitsFor(BaseT base)
  {
    return base == Base2_64 ? 64 : 32;
  }

  // Knuth-style normalization count: bits to shift left so `top` (nonzero)
  // gets its high bit set within a limbBits-wide limb. Replaces the
  // topBitMask while-loops that were copy-pasted across the dividers.
  inline Int NormalizationShiftBits(DataT top, int limbBits)
  {
    return limbBits == 64 ? (Int)__builtin_clzll((unsigned long long)top)
                          : (Int)__builtin_clz((unsigned)top);
  }

  // v << bits at the bit level, bits ∈ [0, limbBits]. Result is trimmed to
  // the canonical (≥ 1 limb) form. 128-bit intermediates make bits == limbBits
  // well-defined for both widths.
  inline std::vector<DataT> ShiftLeftBits(std::span<const DataT> v, Int bits, int limbBits)
  {
    if (bits == 0 || IsZero(v))
      return std::vector<DataT>(v.begin(), v.end());

    const DataT mask = (limbBits == 64) ? (DataT)0xFFFFFFFFFFFFFFFFULL
                                        : (DataT)0xFFFFFFFFULL;
    std::vector<DataT> out(v.size() + 1, 0);
    for (SizeT i = 0; i < (SizeT)v.size(); ++i)
    {
      ULong128 cur = (ULong128)v[i] << bits;
      out[i] |= (DataT)(cur & mask);
      out[i + 1] |= (DataT)((cur >> limbBits) & mask);
    }
    TrimZerosToOne(out);
    return out;
  }

  // v >> bits at the bit level, bits ∈ [0, limbBits). Result is trimmed to
  // the canonical (≥ 1 limb) form.
  inline std::vector<DataT> ShiftRightBits(std::span<const DataT> v, Int bits, int limbBits)
  {
    if (bits == 0 || IsZero(v))
      return std::vector<DataT>(v.begin(), v.end());

    const DataT mask = (limbBits == 64) ? (DataT)0xFFFFFFFFFFFFFFFFULL
                                        : (DataT)0xFFFFFFFFULL;
    const SizeT n = (SizeT)v.size();
    std::vector<DataT> out(n, 0);
    for (SizeT i = 0; i < n; ++i)
    {
      ULong128 cur = (ULong128)v[i];
      if (i + 1 < n)
        cur |= (ULong128)v[i + 1] << limbBits;
      out[i] = (DataT)((cur >> bits) & mask);
    }
    TrimZerosToOne(out);
    return out;
  }

  // r >>= bits in place, bits ∈ [0, limbBits). Does not trim — remainder
  // denormalization wants the limb count preserved.
  inline void ShiftRightBitsInPlace(std::vector<DataT> &r, Int bits, int limbBits)
  {
    if (bits == 0)
      return;
    const DataT mask = (limbBits == 64) ? (DataT)0xFFFFFFFFFFFFFFFFULL
                                        : (DataT)0xFFFFFFFFULL;
    const SizeT n = (SizeT)r.size();
    for (SizeT i = 0; i < n; ++i)
    {
      DataT hi = (i + 1 < n) ? r[i + 1] : 0;
      ULong128 cur = (ULong128)r[i] | ((ULong128)hi << limbBits);
      r[i] = (DataT)((cur >> bits) & mask);
    }
  }
}

#endif
