/**
 * BigInteger Class
 * Top-level square dispatcher: routes to Classic/Karatsuba/NTT square based on size.
 * Mirrors Multiplication.h thresholds.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef SQUARING
#define SQUARING

#include <vector>
using namespace std;

#include "../algorithms/multiplication/ClassicSquare.h"
#include "../algorithms/multiplication/KaratsubaSquare.h"
#include "../algorithms/multiplication/NTTSquare.h"
#include "../algorithms/Multiplication.h"  // for NTT_MULTIPLICATION_THRESHOLD

namespace BigMath
{
#ifndef BIGMATH_NTT_SQUARE_THRESHOLD
  // Default: cross at the same operand-pair size as mult dispatch would,
  // halved because we're checking a.size() not a.size()+b.size().
  // For squaring a×a, both sides equal, so total work scales with a.size() alone.
#define BIGMATH_NTT_SQUARE_THRESHOLD (BIGMATH_NTT_MULTIPLICATION_THRESHOLD / 2)
#endif

  const SizeT NTT_SQUARE_THRESHOLD = BIGMATH_NTT_SQUARE_THRESHOLD;

  inline vector<DataT> Square(vector<DataT> const &a, BaseT base)
  {
    if (IsZero(a))
      return vector<DataT>{0};

    if (a.size() == 1)
      return ClassicSquare::Square(a, base);

    if (a.size() < NTT_SQUARE_THRESHOLD)
      return KaratsubaSquare::Square(a, base);

    return NTTSquare::Square(a, base);
  }
}

#endif
