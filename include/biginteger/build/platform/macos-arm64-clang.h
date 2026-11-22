#pragma once
// Survey run 44997161065; retain the cold-square crossover supported by local and CI comparisons.
#if BIGMATH_LIMB_64 == 1
#ifndef BIGMATH_NTT_SQUARE_COLD_THRESHOLD
#define BIGMATH_NTT_SQUARE_COLD_THRESHOLD 512
#endif
// One-thread equal-length products on M1 Max still favor recursive
// multiplication beyond the parked pool's 512-limb crossover.
#ifndef BIGMATH_NTT_SERIAL_BALANCED_MIN_LIMBS
#define BIGMATH_NTT_SERIAL_BALANCED_MIN_LIMBS 704
#endif
#ifndef BIGMATH_DIV_BZ_MIN_B
#define BIGMATH_DIV_BZ_MIN_B 81
#endif
#endif
