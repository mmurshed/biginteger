/**
 * BigMath: Division dispatcher implementation.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/algorithms/Division.h"

#include <stdexcept>

namespace BigMath
{
  const SizeT NEWTON_LARGE_B = BIGMATH_NEWTON_LARGE_B;
  const SizeT NEWTON_MEDIUM_B = BIGMATH_NEWTON_MEDIUM_B;
  const SizeT BZ_DIVISOR_THRESHOLD = BIGMATH_BZ_DIVISOR_THRESHOLD;
  const SizeT NEWTON_SKEW_NUMERATOR = BIGMATH_NEWTON_SKEW_NUMERATOR;
  const SizeT NEWTON_SKEW_DENOMINATOR = BIGMATH_NEWTON_SKEW_DENOMINATOR;

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
    bool newton_eligible =
        b.size() >= NEWTON_MEDIUM_B &&
        NEWTON_SKEW_DENOMINATOR * a.size() >= NEWTON_SKEW_NUMERATOR * b.size();
    if (newton_eligible)
      return NewtonDivision::DivideAndRemainder(a, b, base, computeRemainder);

    // BZ for large near-balanced divisors and for big-and-skewed cases.
    bool bz_eligible =
        base == Base2_32 &&
        b.size() > BZ_DIVISOR_THRESHOLD &&
        ((b.size() >= 1024 && a.size() <= 3 * b.size()) ||
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
