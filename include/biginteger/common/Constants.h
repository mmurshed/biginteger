/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_CONSTANTS
#define BIGINTEGER_CONSTANTS

#include <cstdint>

// Limb width selector. When 1 (default), DataT stores true 64-bit values and
// Base() == Base2_64. When 0, DataT stores 32-bit values in a 64-bit container
// (the historical layout) and Base() == Base2_32. The 32-bit path is kept
// available for A/B testing and as a fallback; flip via -DBIGMATH_LIMB_64=0.
//
// Bench (M1 Max, -O3 -march=native, bench_vs_gmp): LIMB_64=1 matches or beats
// LIMB_64=0 on every measured op. Wins: small/mid mul -25 to -60%, skewed div
// -36 to -59%, parse -24 to -34%, ToString -19 to -28%. 1M NTT-bound mul is
// flat (NTT-coefficient-bound, not limb-bound).
#ifndef BIGMATH_LIMB_64
#define BIGMATH_LIMB_64 1
#endif

// Internal parallelism for NTT-bound operations. When 1 (default), the
// CRT NTT path dispatches its 6 forwards + 3 inverses as batched work
// units across a small thread pool; pointwise multiply is also chunked
// across workers. When 0, all calls reduce to inline serial loops — zero
// overhead and no pthread linkage pulled into the binary.
//
// Pool size auto-detected at runtime, capped at BIGMATH_MAX_THREADS.
// See docs/THREAD_SAFETY.md for the thread-safety model.
//
// Bench (M1 Max, default LIMB_64=1 + CRT default): 2.3-3.4× speedup on
// large mul / skewed div / parse vs the single-thread path. Opt out via
// -DBIGMATH_USE_THREADS=0 for embedded/header-only-strict consumers.
#ifndef BIGMATH_USE_THREADS
#define BIGMATH_USE_THREADS 1
#endif

#ifndef BIGMATH_MAX_THREADS
#define BIGMATH_MAX_THREADS 8
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
