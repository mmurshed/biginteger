/**
 * BigMath: Squaring dispatcher
 *
 * Top-level square dispatcher: routes to Classic/Karatsuba/NTT square based on
 * operand size. Mirrors Multiplication.h thresholds but uses operand size
 * directly (not the sum), since for a · a both sides are equal.
 *
 * NTT_SQUARE_THRESHOLD defaults to NTT_MULTIPLICATION_THRESHOLD / 2 = 750 limbs.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef SQUARING
#define SQUARING

#include <vector>

#include "../algorithms/multiplication/ClassicSquare.h"
#include "../algorithms/multiplication/KaratsubaSquare.h"
#include "../algorithms/multiplication/NTTSquare.h"
#include "../algorithms/Multiplication.h"  // for BIGMATH_NTT_MULTIPLICATION_THRESHOLD

namespace BigMath
{
#ifndef BIGMATH_NTT_SQUARE_THRESHOLD
#define BIGMATH_NTT_SQUARE_THRESHOLD (BIGMATH_NTT_MULTIPLICATION_THRESHOLD / 2)
#endif

  extern const SizeT NTT_SQUARE_THRESHOLD;

  std::vector<DataT> Square(std::vector<DataT> const &a, BaseT base);
}

#endif
