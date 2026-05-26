/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_CONSTANTS
#define BIGINTEGER_CONSTANTS

#include <cstdint>

// Limb width selector. When 0 (default), DataT stores 32-bit values in a 64-bit
// container — the historical layout. When 1, DataT stores true 64-bit values
// and Base() switches to Base2_64. Foundation flag for the multi-PR 64-bit-limb
// refactor; subsequent PRs add the 64-bit code paths gated on this macro. The
// 32-bit path remains the default until the refactor is complete and tested.
#ifndef BIGMATH_LIMB_64
#define BIGMATH_LIMB_64 0
#endif

namespace BigMath
{
  typedef int64_t Long;
  typedef uint64_t ULong;

  typedef int32_t Int;
  typedef uint32_t UInt;

  typedef __int128 Long128;
  typedef unsigned __int128 ULong128;

  typedef uint64_t DataT;
  typedef uint32_t SizeT;
  typedef int64_t BaseT;

    // Base 10
    const BaseT Base10 = 10;
    // Base 10^3
    const BaseT Base100 = 100;
    // Base 10^8
    const BaseT Base100M = 100000000lu;
    // Number of digits in `BASE'
    const SizeT Base100_Zeroes = 2;

    // Base 2
    const BaseT Base2 = 2;
    // Base 16
    const BaseT Base16 = 16;
    // Base 2^8
    const BaseT Base2_8 = 256;
    // Base 2^10
    const BaseT Base2_10 = 1024;
    // Base 2^16
    const BaseT Base2_16 = 65536ul;
    // Base 2^32
    const BaseT Base2_32 = 4294967296ul; // 2^32
    // Base 2^31
    const BaseT Base2_31 = 2147483648; // 2^31

    // Sentinel for Base 2^64. The numeric value 2^64 does not fit in BaseT
    // (int64_t), so subsequent PRs use the sentinel 0 in dispatch — 0 is never
    // a legal numeric base. Code paths that need to *use* the base value as a
    // shift/mask (rather than just compare it) take a separate Base2_64 branch
    // and operate on the limb directly via the LIMB_* helpers below.
    const BaseT Base2_64 = 0;

    // Per-limb width helpers. These drive the carry/borrow/shift idioms in
    // Addition, Subtraction, ClassicMultiplication, FastDivision, etc. The
    // helpers compile to either the historical 32-bit ops or the new 64-bit
    // ops depending on BIGMATH_LIMB_64.
#if BIGMATH_LIMB_64
    constexpr SizeT LimbBits = 64;
    constexpr ULong128 LimbBase = (ULong128)1 << 64; // 2^64
    constexpr ULong LimbMask = 0xFFFFFFFFFFFFFFFFULL;
    // Compile-time current Base() value, returned by BigInteger::Base().
    constexpr BaseT CurrentBase = Base2_64;
#else
    constexpr SizeT LimbBits = 32;
    constexpr ULong LimbBase = 1ULL << 32; // 2^32
    constexpr ULong LimbMask = 0xFFFFFFFFULL;
    constexpr BaseT CurrentBase = Base2_32;
#endif
}

#endif
