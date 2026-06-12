/**
 * BigMath: Squaring dispatcher
 *
 * Top-level square dispatcher: routes to Classic/Karatsuba/NTT square based on
 * operand size. Mirrors Multiplication.h thresholds but uses operand size
 * directly (not the sum), since for a · a both sides are equal.
 *
 * NTT_SQUARE_THRESHOLD defaults to 640 limbs (LIMB_64) based on the
 * KaratsubaSquare vs CRT-NTT self-multiply crossover.
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
// CRT-square crossover (2026-06-12): KaratsubaSquare ties the CRT+NEON
// self-multiply at 512 limbs, loses 1.3x+ from 768. Under LIMB_64=0 the
// crossover stays at 512. Canonical copy: build/DispatchThresholds.h.
#if BIGMATH_LIMB_64
#define BIGMATH_NTT_SQUARE_THRESHOLD 640
#else
#define BIGMATH_NTT_SQUARE_THRESHOLD 512
#endif
#endif

  extern const SizeT NTT_SQUARE_THRESHOLD;

  std::vector<DataT> Square(std::vector<DataT> const &a, BaseT base);
}

#endif
