/**
 * BigMath: Division dispatcher implementation.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/algorithms/Division.h"
#include "biginteger/algorithms/division/QuotientSizedDivision.h"

#include <stdexcept>

namespace BigMath
{
  const SizeT NEWTON_MEDIUM_B = BIGMATH_NEWTON_MEDIUM_B;
  const SizeT NEWTON_RATIO2_B = BIGMATH_NEWTON_RATIO2_B;
  const SizeT NEWTON_RATIO2_NUMERATOR = BIGMATH_NEWTON_RATIO2_NUMERATOR;
  const SizeT NEWTON_RATIO2_DENOMINATOR = BIGMATH_NEWTON_RATIO2_DENOMINATOR;
  const SizeT BZ_DIVISOR_THRESHOLD = BIGMATH_BZ_DIVISOR_THRESHOLD;
  const SizeT NEWTON_SKEW_NUMERATOR = BIGMATH_NEWTON_SKEW_NUMERATOR;
  const SizeT NEWTON_SKEW_DENOMINATOR = BIGMATH_NEWTON_SKEW_DENOMINATOR;
  const SizeT NEWTON_RATIO35_B = BIGMATH_NEWTON_RATIO35_B;
  const SizeT NEWTON_RATIO35_NUMERATOR = BIGMATH_NEWTON_RATIO35_NUMERATOR;
  const SizeT NEWTON_RATIO35_DENOMINATOR = BIGMATH_NEWTON_RATIO35_DENOMINATOR;
  const SizeT NEWTON_MID_B = BIGMATH_NEWTON_MID_B;
  const SizeT NEWTON_MID_NUMERATOR = BIGMATH_NEWTON_MID_NUMERATOR;
  const SizeT NEWTON_MID_DENOMINATOR = BIGMATH_NEWTON_MID_DENOMINATOR;
  const SizeT NEWTON_BALANCED_B = BIGMATH_NEWTON_BALANCED_B;
  const SizeT NEWTON_BALANCED_NUMERATOR = BIGMATH_NEWTON_BALANCED_NUMERATOR;
  const SizeT NEWTON_BALANCED_DENOMINATOR = BIGMATH_NEWTON_BALANCED_DENOMINATOR;
  const SizeT QSIZED_MIN_DELTA = BIGMATH_QSIZED_MIN_DELTA;
  const SizeT QSIZED_SMALL_B = BIGMATH_QSIZED_SMALL_B;
  const SizeT QSIZED_SMALL_DELTA_DIV = BIGMATH_QSIZED_SMALL_DELTA_DIV;
  const SizeT NEWTON_HIGH_SKEW_B = BIGMATH_NEWTON_HIGH_SKEW_B;
  const SizeT NEWTON_HIGH_SKEW_NUMERATOR = BIGMATH_NEWTON_HIGH_SKEW_NUMERATOR;
  const SizeT NEWTON_HIGH_SKEW_DENOMINATOR = BIGMATH_NEWTON_HIGH_SKEW_DENOMINATOR;

  std::pair<std::vector<DataT>, std::vector<DataT>> DivideAndRemainder(
      std::vector<DataT> const &a,
      std::vector<DataT> const &b,
      BaseT base,
      bool computeRemainder)
  {
    if (IsZero(b))
      throw std::invalid_argument("Division by zero");

    if (IsZero(a))
    {
      auto q = std::vector<DataT>{0};
      return {q, q};
    }

    Int cmp = Compare(a, b);
    if (cmp == 0)
      return {std::vector<DataT>{1}, std::vector<DataT>{0}};
    if (cmp < 0)
      return {std::vector<DataT>{0}, computeRemainder ? a : std::vector<DataT>()};

    // Newton handles any ratio via blockwise mode; pick when divisor is large enough
    // for reciprocal-setup amortization and skew is in band.
    bool newton_medium_skew =
        b.size() >= NEWTON_MEDIUM_B &&
        NEWTON_SKEW_DENOMINATOR * a.size() >= NEWTON_SKEW_NUMERATOR * b.size();
    // Ratio ≥ 7/2 from 1024: Newton wins ratio ≥ 4 across [1024, 1280).
    bool newton_ratio35 =
        b.size() >= NEWTON_RATIO35_B &&
        NEWTON_RATIO35_DENOMINATOR * a.size() >= NEWTON_RATIO35_NUMERATOR * b.size();
    // Exactly-2.5 band: padded BZ beats Newton below NEWTON_MID_B at ratio
    // 2.5, while the medium band's 14/5 keeps the ratio-3 knife-edge.
    bool newton_mid_skew =
        b.size() >= NEWTON_MID_B &&
        NEWTON_MID_DENOMINATOR * a.size() >= NEWTON_MID_NUMERATOR * b.size();
    // Near-balanced (ratio ≥ 4/3) band: only above NEWTON_BALANCED_B, where BZ's
    // near-balanced path degrades erratically (measured 2×–4.5× slower than
    // Newton at b ≥ 100k limbs); below it BZ wins, so leave it alone.
    bool newton_balanced =
        b.size() >= NEWTON_BALANCED_B &&
        NEWTON_BALANCED_DENOMINATOR * a.size() >= NEWTON_BALANCED_NUMERATOR * b.size();
    bool newton_high_skew =
        b.size() >= NEWTON_HIGH_SKEW_B &&
        NEWTON_HIGH_SKEW_DENOMINATOR * a.size() >= NEWTON_HIGH_SKEW_NUMERATOR * b.size();
    // Ratio-≥8/5 band between medium (3/1 @ 2560) and balanced (4/3 @ 24576).
    bool newton_ratio2 =
        b.size() >= NEWTON_RATIO2_B &&
        NEWTON_RATIO2_DENOMINATOR * a.size() >= NEWTON_RATIO2_NUMERATOR * b.size();
    bool newton_eligible = newton_medium_skew || newton_ratio35 || newton_mid_skew || newton_balanced || newton_high_skew || newton_ratio2;
    if (newton_eligible)
      return NewtonDivision::DivideAndRemainder(a, b, base, computeRemainder);

    // Short-quotient band: large divisors below the Newton balanced band
    // (ratio < 4/3) with a quotient big enough that FastDivision's O(n*delta)
    // loses. BZ's near-balanced path blows up 7-128x on 2^k+1-family divisor
    // sizes here; quotient-sized division scales with the quotient instead.
    bool qsized_main =
        b.size() >= NEWTON_BALANCED_B &&
        a.size() >= b.size() + QSIZED_MIN_DELTA &&
        NEWTON_BALANCED_DENOMINATOR * a.size() < NEWTON_BALANCED_NUMERATOR * b.size();
    // Thin-quotient extension below the balanced floor: delta <= b/8. Generic
    // wins are modest; the point is the 3-14x BZ pathology on 2^k+1-family
    // divisor sizes in this band.
    bool qsized_thin =
        b.size() >= QSIZED_SMALL_B &&
        a.size() >= b.size() + QSIZED_MIN_DELTA &&
        (a.size() - b.size()) * QSIZED_SMALL_DELTA_DIV <= b.size();
    bool qsized_eligible =
        (base == Base2_32 || base == Base2_64) && (qsized_main || qsized_thin);
    if (qsized_eligible)
      return QuotientSizedDivision::DivideAndRemainder(a, b, base, computeRemainder);

    // BZ for large near-balanced divisors and for big-and-skewed cases.
    // The +32-limb quotient-bulk guard in the near-balanced clause excludes
    // degenerate cases where a ≈ b and the quotient is 0-2 limbs — BZ would
    // split a into m = n/2 blocks and run wasted m×m multiplies on mostly-zero
    // high blocks, while FastDivision short-circuits via a single qhat
    // iteration. Regressed 5M×5M balanced 1.45 → 4.48 ms before this guard.
    bool bz_eligible =
        (base == Base2_32 || base == Base2_64) &&
        b.size() > BZ_DIVISOR_THRESHOLD &&
        ((b.size() >= 1024 && a.size() >= b.size() + 32 && a.size() <= 3 * b.size()) ||
         (a.size() > 2048 && a.size() > 3 * b.size()));
    if (bz_eligible)
      return BurnikelZieglerDivision::DivideAndRemainder(a, b, base, computeRemainder);

    return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);
  }

  std::vector<DataT> Divide(std::vector<DataT> const &a,
                            std::vector<DataT> const &b,
                            BaseT base)
  {
    return DivideAndRemainder(a, b, base, false).first;
  }

  std::pair<std::vector<DataT>, std::vector<DataT>> DivideAndRemainder(
      std::vector<DataT> const &a,
      DataT b,
      BaseT base)
  {
    if (b == 0)
      throw std::invalid_argument("Division by zero");

    if (IsZero(a))
    {
      auto q = std::vector<DataT>{0};
      return {q, q};
    }

    Int cmp = Compare(a, b);
    if (cmp == 0)
      return {std::vector<DataT>{1}, std::vector<DataT>{0}};
    if (cmp < 0)
      return {std::vector<DataT>{0}, a};

    return ClassicDivision::DivideAndRemainder(a, b, base);
  }

  std::vector<DataT> Divide(std::vector<DataT> const &a,
                            DataT b,
                            BaseT base)
  {
    if (b == 0)
      throw std::invalid_argument("Division by zero");

    if (IsZero(a))
      return std::vector<DataT>{0};

    Int cmp = Compare(a, b);
    if (cmp == 0)
      return std::vector<DataT>{1};
    if (cmp < 0)
      return std::vector<DataT>{0};

    return ClassicDivision::Divide(a, b, base);
  }
}
