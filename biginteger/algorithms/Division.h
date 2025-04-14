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
#include "division/MontgomeryDivision.h"

namespace BigMath
{
  pair<vector<DataT>, vector<DataT>> DivideAndRemainder(vector<DataT> const &a, vector<DataT> const &b, BaseT base)
  {
      // case of 0 divisor
      if (IsZero(b))
      {
        throw invalid_argument("Division by zero");
      }

      // case of 0 dividend
      if (IsZero(b))
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
      auto result = MontgomeryDivision::DivideAndRemainder(a, b, base);

      return result;
  }
}

#endif