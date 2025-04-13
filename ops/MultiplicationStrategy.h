/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_MULTIPLICATION_STRATEGY
#define BIGINTEGER_MULTIPLICATION_STRATEGY

#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "../algorithms/classic/ClassicMultiplication.h"
#include "../algorithms/karatsuba/KaratsubaMultiplication.h"
#include "../algorithms/toomcook/ToomCookMultiplication2.h"
#include "../algorithms/stonehagestrassen/FFTMultiplication.h"

namespace BigMath
{
  class MultiplicationStrategy
  {
    // Karatsuba performs better for over 128 digits for the result
    static const SizeT KARATSUBA_THRESHOLD = 128;
    static const SizeT TOOM_COOK_THRESHOLD = 512;
    static const SizeT FFT_THRESHOLD = 1024;

  public:
    static vector<DataT> MultiplyUnsigned(vector<DataT> const &a, vector<DataT> const &b, BaseT base)
    {
      SizeT size = a.size() + b.size();

      if (size <= KARATSUBA_THRESHOLD)
      {
        return ClassicMultiplication::Multiply(a, b, base);
      }
      else if (size <= TOOM_COOK_THRESHOLD)
      {
        return KaratsubaMultiplication::Multiply(a, b, base);
      }
      else if (size <= FFT_THRESHOLD)
      {
        ToomCookMultiplication2 toomCook;
        return toomCook.Multiply(a, b, base);
      }

      return FFTMultiplication::Multiply(a, b, base);
    }
  };
}

#endif