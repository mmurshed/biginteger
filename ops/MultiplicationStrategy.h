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
#include "../algorithms/multiplication/ClassicMultiplication.h"
#include "../algorithms/multiplication/KaratsubaMultiplication.h"
#include "../algorithms/multiplication/ToomCookMultiplication.h"
#include "../algorithms/multiplication/FFTMultiplication.h"

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
        return ToomCookMultiplication::Multiply(a, b, base);
      }

      return FFTMultiplication::Multiply(a, b, base);
    }
  };
}

#endif