#ifndef RECIPROCAL_DIVISION
#define RECIPROCAL_DIVISION

#include <stdexcept>
#include <utility>
#include <vector>
using namespace std;

#include "NewtonDivision.h"

namespace BigMath
{
  class ReciprocalDivision
  {
  public:
    class Divider
    {
    private:
      NewtonDivision::Divider divider;

    public:
      Divider(vector<DataT> const &b, BaseT radix) : divider(b, radix) {}

      pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
          vector<DataT> const &a,
          bool computeRemainder = true) const
      {
        return divider.DivideAndRemainder(a, computeRemainder);
      }

      vector<DataT> Divide(vector<DataT> const &a) const
      {
        return divider.Divide(a);
      }

      vector<DataT> const &Divisor() const
      {
        return divider.Divisor();
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
