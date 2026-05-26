/**
 * BigMath: Squaring dispatcher
 *
 * Top-level square dispatcher: routes to Classic/Karatsuba/NTT square based on
 * operand size. Mirrors Multiplication.h thresholds but uses operand size
 * directly (not the sum), since for a · a both sides are equal.
 *
 * NTT_SQUARE_THRESHOLD defaults to 512 limbs based on the square-specific
 * KaratsubaSquare vs NTTSquare crossover.
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
// dispatch_tuner LIMB_64=1: KaratsubaSquare wins through ~1536 limbs; NTT
// takes over at 2048. Under LIMB_64=0 crossover stays at 512.
#if BIGMATH_LIMB_64
#define BIGMATH_NTT_SQUARE_THRESHOLD 2048
#else
#define BIGMATH_NTT_SQUARE_THRESHOLD 512
#endif
#endif

  extern const SizeT NTT_SQUARE_THRESHOLD;

  std::vector<DataT> Square(std::vector<DataT> const &a, BaseT base);
}

#endif
