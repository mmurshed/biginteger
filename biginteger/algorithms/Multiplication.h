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
#include "../algorithms/multiplication/NTTMultiplication.h"

namespace BigMath
{
#ifndef BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD
#define BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD 32
#endif

#ifndef BIGMATH_NTT_MULTIPLICATION_THRESHOLD
#define BIGMATH_NTT_MULTIPLICATION_THRESHOLD 512
#endif

  const SizeT CLASSIC_MULTIPLICATION_THRESHOLD = BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD;
  const SizeT NTT_MULTIPLICATION_THRESHOLD = BIGMATH_NTT_MULTIPLICATION_THRESHOLD;

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
    SizeT minSize = min(a.size(), b.size());

    if (size <= CLASSIC_MULTIPLICATION_THRESHOLD || minSize <= 32)
      return ClassicMultiplication::Multiply(a, b, base);

    if (size < NTT_MULTIPLICATION_THRESHOLD)
      return KaratsubaMultiplication::Multiply(a, b, base);

    return NTTMultiplication::Multiply(a, b, base);
  }

  vector<DataT> Multiply(vector<DataT> const &a, DataT b, BaseT base)
  {
    if (b == 0 || IsZero(a))
      return vector<DataT>{0}; // 0 times anything is zero
    return ClassicMultiplication::Multiply(a, b, base);
  }
}

#endif
