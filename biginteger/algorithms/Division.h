/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef DIVISION
#define DIVISION

#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "division/ClassicDivision.h"
#include "division/KnuthDivision.h"

namespace BigMath
{
  pair<vector<DataT>, vector<DataT>> DivideAndRemainder(vector<DataT> const &a, vector<DataT> const &b, BaseT base, bool computerRemainder = true)
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
        return {vector<DataT>{0}, b}; // case of a < b
      }

      // Now: a > b
      return KnuthDivision::DivideAndRemainder(a, b, base, computerRemainder);
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
        return {vector<DataT>{0}, vector<DataT>{b}}; // case of a < b
      }

      // Now: a > b
      return KnuthDivision::DivideAndRemainder(a, b, base);
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
      return KnuthDivision::Divide(a, b, base);
  }  
}

#endif