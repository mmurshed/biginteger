/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef DIVISION
#define DIVISION

#include <stdexcept>
#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "division/BurnikelZieglerDivision.h"
#include "division/ClassicDivision.h"
#include "division/FastDivision.h"
#include "division/NewtonDivision.h"

namespace BigMath
{
#ifndef BIGMATH_NEWTON_LARGE_B
#define BIGMATH_NEWTON_LARGE_B 24576
#endif

#ifndef BIGMATH_NEWTON_MEDIUM_B
#define BIGMATH_NEWTON_MEDIUM_B 16384
#endif

#ifndef BIGMATH_BZ_DIVISOR_THRESHOLD
#define BIGMATH_BZ_DIVISOR_THRESHOLD 512
#endif

  // Newton-Raphson is O(M(n)) vs Knuth's O(m*n). Blockwise mode in NewtonDivision lets it
  // handle ANY ratio na/nb internally — the only thing the dispatcher gates is whether
  // Newton's reciprocal-setup cost amortizes. Two bands:
  //   (a) Divisor large enough that reciprocal multiplies hit NTT (b ≥ 24576 limbs) — Newton
  //       wins even at near-square ratio 4/3.
  //   (b) Smaller divisors (16384 ≤ b < 24576) — Newton only wins at ratio ≥ 2 because
  //       the multiplies inside the reciprocal sit on the Karatsuba/NTT boundary.
  const SizeT NEWTON_LARGE_B = BIGMATH_NEWTON_LARGE_B;
  const SizeT NEWTON_MEDIUM_B = BIGMATH_NEWTON_MEDIUM_B;
  const SizeT BZ_DIVISOR_THRESHOLD = BIGMATH_BZ_DIVISOR_THRESHOLD;

  pair<vector<DataT>, vector<DataT>> DivideAndRemainder(vector<DataT> const &a, vector<DataT> const &b, BaseT base, bool computeRemainder = true)
  {
      // case of 0 divisor
      if (IsZero(b))
      {
        throw invalid_argument("Division by zero");
      }

      // case of 0 dividend
      if (IsZero(a))
      {
        auto q = vector<DataT>{0};
        return {q, q};
      }

      Int cmp = Compare(a, b);
      if (cmp == 0)
      {
        return {vector<DataT>{1}, vector<DataT>{0}}; // case of a/a
      }
      else if (cmp < 0)
      {
        return {vector<DataT>{0}, computeRemainder ? a : vector<DataT>()}; // case of a < b
      }

      // Now: a > b. Newton handles any ratio via blockwise; pick it whenever b is in either
      // of the win bands and the dividend isn't too close to square.
      bool newton_eligible =
          (b.size() >= NEWTON_LARGE_B && 3 * a.size() >= 4 * b.size()) ||
          (b.size() >= NEWTON_MEDIUM_B && a.size() >= 2 * b.size());
      if (newton_eligible)
      {
        return NewtonDivision::DivideAndRemainder(a, b, base, computeRemainder);
      }

      // Use BZ for large even-size divisors where its balanced 2n/n recursion applies.
      // Keep the skewed path too, but avoid tiny divisors where BZ immediately falls back.
      bool bz_eligible =
          base == Base2_32 &&
          b.size() > BZ_DIVISOR_THRESHOLD &&
          b.size() % 2 == 0 &&
          ((b.size() >= 1024 && a.size() <= 3 * b.size()) ||
           (a.size() > 2048 && a.size() > 3 * b.size()));
      if (bz_eligible)
      {
        return BurnikelZieglerDivision::DivideAndRemainder(a, b, base, computeRemainder);
      }

      return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);
  }

  vector<DataT> Divide(vector<DataT> const &a, vector<DataT> const &b, BaseT base)
  {
    return DivideAndRemainder(a, b, base, false).first;
  }

  pair<vector<DataT>, vector<DataT>> DivideAndRemainder(vector<DataT> const &a, DataT b, BaseT base)
  {
      // case of 0 divisor
      if (b == 0)
      {
        throw invalid_argument("Division by zero");
      }

      // case of 0 dividend
      if (IsZero(a))
      {
        auto q = vector<DataT>{0};
        return {q, q};
      }

      Int cmp = Compare(a, b);
      if (cmp == 0)
      {
        return {vector<DataT>{1}, vector<DataT>{0}}; // case of a/a
      }
      else if (cmp < 0)
      {
        return {vector<DataT>{0}, a}; // case of a < b
      }

      // Now: a > b
      return ClassicDivision::DivideAndRemainder(a, b, base);
  }

  vector<DataT> Divide(vector<DataT> const &a, DataT b, BaseT base)
  {
      // case of 0 divisor
      if (b == 0)
      {
        throw invalid_argument("Division by zero");
      }

      // case of 0 dividend
      if (IsZero(a))
      {
        return vector<DataT>{0};
      }

      Int cmp = Compare(a, b);
      if (cmp == 0)
      {
        return vector<DataT>{1}; // case of a/a
      }
      else if (cmp < 0)
      {
        return vector<DataT>{0}; // case of a < b
      }

      // Now: a > b
      return ClassicDivision::Divide(a, b, base);
  }  
}

#endif
