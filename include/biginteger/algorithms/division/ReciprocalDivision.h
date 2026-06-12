#ifndef RECIPROCAL_DIVISION
#define RECIPROCAL_DIVISION

#include <stdexcept>
#include <utility>
#include <vector>

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
      Divider(std::vector<DataT> const &b, BaseT radix) : divider(b, radix) {}

      std::pair<std::vector<DataT>, std::vector<DataT>> DivideAndRemainder(
          std::vector<DataT> const &a,
          bool computeRemainder = true) const
      {
        return divider.DivideAndRemainder(a, computeRemainder);
      }

      std::vector<DataT> Divide(std::vector<DataT> const &a) const
      {
        return divider.Divide(a);
      }

      std::vector<DataT> const &Divisor() const
      {
        return divider.Divisor();
      }
    };

    static std::pair<std::vector<DataT>, std::vector<DataT>> DivideAndRemainder(
        std::vector<DataT> const &a,
        std::vector<DataT> const &b,
        BaseT base)
    {
      return Divider(b, base).DivideAndRemainder(a);
    }

    static std::vector<DataT> Divide(
        std::vector<DataT> const &a,
        std::vector<DataT> const &b,
        BaseT base)
    {
      return Divider(b, base).Divide(a);
    }
  };
}

#endif
