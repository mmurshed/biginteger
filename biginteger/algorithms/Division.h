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
  // Newton-Raphson is O(M(n)) vs Knuth's O(m*n). Wins when the total work product na*nb
  // is large enough that Newton's reciprocal-setup constant (a handful of large multiplies)
  // amortizes, AND the ratio na/nb sits in a band where neither end-case dominates:
  // square (Newton wastes a full reciprocal on a 1-2 limb quotient) or skewed (BZ already wins).
  const ULong NEWTON_DIVISION_WORK_THRESHOLD = (ULong)16384 * 32768;

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

      // Now: a > b
      if (a.size() > 2048 && a.size() > 3 * b.size())
      {
        return BurnikelZieglerDivision::DivideAndRemainder(a, b, base, computeRemainder);
      }

      // Newton sweet spot: ratio na/nb in [4/3, 5/2] and work above threshold.
      if ((ULong)a.size() * b.size() >= NEWTON_DIVISION_WORK_THRESHOLD &&
          3 * a.size() >= 4 * b.size() &&
          2 * a.size() <= 5 * b.size())
      {
        return NewtonDivision::DivideAndRemainder(a, b, base, computeRemainder);
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
