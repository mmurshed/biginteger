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

#include "../algorithms/multiplication/ClassicMultiplication.h"
#include "../algorithms/multiplication/KaratsubaMultiplication.h"
#include "../algorithms/multiplication/ToomCookMultiplication.h"
#include "../algorithms/multiplication/FFTMultiplication.h"

namespace BigMath
{
  const SizeT FFT_THRESHOLD = 700;

  vector<DataT> Multiply(vector<DataT> const &a, vector<DataT> const &b, BaseT base)
  {
    if (IsZero(a) || IsZero(b))
      return vector<DataT>{0}; // 0 times anything is zero

    // If b is a single digit, use the scalar multiplication
    if(b.size() == 1)
      return ClassicMultiplication::Multiply(a, b[0], base);
    // If a is a single digit, use the scalar multiplication
    if(a.size() == 1)
      return ClassicMultiplication::Multiply(b, a[0], base);

    SizeT size = a.size() + b.size();

    if (size <= FFT_THRESHOLD)
    {
      return ToomCookMultiplication::Multiply(a, b, base);
    }

    return FFTMultiplication::Multiply(a, b, base);
  }

  vector<DataT> Multiply(vector<DataT> const &a, DataT b, BaseT base)
  {
    if (b == 0 || IsZero(a))
      return vector<DataT>{0}; // 0 times anything is zero
    return ClassicMultiplication::Multiply(a, b, base);
  }
}

#endif