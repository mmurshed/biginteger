#ifndef BURNIKELZIEGLER_DIVISION
#define BURNIKELZIEGLER_DIVISION

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
  class BurnikelZieglerDivision
  {
  private:
    static const SizeT BZ_THRESHOLD = 512;

    static vector<DataT> ShiftBase(vector<DataT> const &value, SizeT limbs)
    {
      return ShiftLeft(value, limbs);
    }

    static pair<vector<DataT>, vector<DataT>> DivideRecursive(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base,
        bool computeRemainder)
    {
      if (a.size() <= BZ_THRESHOLD || b.size() <= BZ_THRESHOLD || a.size() <= 2 * b.size())
        return FastDivision::DivideAndRemainder(a, b, base, computeRemainder);

      SizeT lowSize = 2 * (SizeT)b.size();
      vector<DataT> low(a.begin(), a.begin() + lowSize);
      vector<DataT> high(a.begin() + lowSize, a.end());
      TrimZeros(low);
      TrimZeros(high);

      auto highQr = DivideRecursive(high, b, base, true);

      vector<DataT> combined = Add(ShiftBase(highQr.second, lowSize), low, base);
      TrimZeros(combined);

      auto lowQr = FastDivision::DivideAndRemainder(combined, b, base, computeRemainder);

      vector<DataT> quotient = Add(ShiftBase(highQr.first, lowSize), lowQr.first, base);
      TrimZeros(quotient);
      if (quotient.empty())
        quotient.push_back(0);

      return {quotient, computeRemainder ? lowQr.second : vector<DataT>()};
    }

  public:
    static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base,
        bool computeRemainder = true)
    {
      if (IsZero(b))
        throw invalid_argument("Division by zero");

      if (IsZero(a))
        return {vector<DataT>{0}, computeRemainder ? vector<DataT>{0} : vector<DataT>()};

      Int cmp = Compare(a, b);
      if (cmp < 0)
        return {vector<DataT>{0}, computeRemainder ? a : vector<DataT>()};
      if (cmp == 0)
        return {vector<DataT>{1}, computeRemainder ? vector<DataT>{0} : vector<DataT>()};

      return DivideRecursive(a, b, base, computeRemainder);
    }

    static vector<DataT> Divide(vector<DataT> const &a, vector<DataT> const &b, BaseT base)
    {
      return DivideAndRemainder(a, b, base, false).first;
    }
  };
}

#endif
