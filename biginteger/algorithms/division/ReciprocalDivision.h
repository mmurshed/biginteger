#ifndef RECIPROCAL_DIVISION
#define RECIPROCAL_DIVISION

#include <stdexcept>
#include <utility>
#include <vector>
using namespace std;

#include "../../common/Comparator.h"
#include "../../common/Util.h"
#include "BurnikelZieglerDivision.h"
#include "ClassicDivision.h"
#include "FastDivision.h"

namespace BigMath
{
  class ReciprocalDivision
  {
  public:
    class Divider
    {
    private:
      vector<DataT> divisor;
      BaseT base;

    public:
      Divider(vector<DataT> const &b, BaseT radix) : divisor(b), base(radix)
      {
        TrimZeros(divisor);
        if (IsZero(divisor))
          throw invalid_argument("Division by zero");
      }

      pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
          vector<DataT> const &a,
          bool computeRemainder = true) const
      {
        if (IsZero(a))
          return {vector<DataT>{0}, computeRemainder ? vector<DataT>{0} : vector<DataT>()};

        Int cmp = Compare(a, divisor);
        if (cmp < 0)
          return {vector<DataT>{0}, computeRemainder ? a : vector<DataT>()};
        if (cmp == 0)
          return {vector<DataT>{1}, computeRemainder ? vector<DataT>{0} : vector<DataT>()};

        if (divisor.size() == 1)
          return ClassicDivision::DivideAndRemainder(a, divisor[0], base);

        if (a.size() > 2048 && a.size() > 3 * divisor.size())
          return BurnikelZieglerDivision::DivideAndRemainder(a, divisor, base, computeRemainder);

        return FastDivision::DivideAndRemainder(a, divisor, base, computeRemainder);
      }

      vector<DataT> Divide(vector<DataT> const &a) const
      {
        return DivideAndRemainder(a, false).first;
      }

      vector<DataT> const &Divisor() const
      {
        return divisor;
      }
    };

    static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base)
    {
      return Divider(b, base).DivideAndRemainder(a);
    }

    static vector<DataT> Divide(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base)
    {
      return Divider(b, base).Divide(a);
    }
  };
}

#endif
