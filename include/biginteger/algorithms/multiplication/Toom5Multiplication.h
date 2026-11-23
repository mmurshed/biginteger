// Direct-call Toom-5 interface. Public integer multiplication uses its
// platform dispatcher.
#pragma once

#include "ToomHigherMultiplication.h"

#ifndef BIGMATH_TOOM5_THRESHOLD
#define BIGMATH_TOOM5_THRESHOLD 512
#endif

namespace BigMath
{
  class Toom5Multiplication
  {
  public:
    static std::vector<DataT> Multiply(const std::vector<DataT> &a, const std::vector<DataT> &b, BaseT base)
    {
      if (IsZero(a) || IsZero(b))
        return {0};
      if (std::min(a.size(), b.size()) == 1 ||
          std::max(a.size(), b.size()) < BIGMATH_TOOM5_THRESHOLD || base != Base2_64)
        return KaratsubaMultiplication::Multiply(a, b, base);
      return ToomHigherMultiplication::Multiply(a, b, base, ToomHigherMultiplication::Split::Toom5);
    }
  };
}
