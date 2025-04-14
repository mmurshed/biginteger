/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef MULTIPLICATION
#define MULTIPLICATION

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
  // Karatsuba performs better for over 128 digits for the result
  const SizeT KARATSUBA_THRESHOLD = 128;
  const SizeT TOOM_COOK_THRESHOLD = 512;
  const SizeT FFT_THRESHOLD = 1024;

  vector<DataT> Multiply(vector<DataT> const &a, vector<DataT> const &b, BaseT base)
  {
    if (IsZero(a) || IsZero(b))
      return vector<DataT>{0}; // 0 times anything is zero

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
}

#endif