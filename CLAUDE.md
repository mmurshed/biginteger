# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

C++ BigInteger library in namespace `BigMath`, built with CMake. Public headers live under `include/biginteger/` (plus `include/bigdecimal/`); non-inline implementations (dispatchers, parser, thread pool, ops) live under `src/`. Consumers link `bigmath::bigmath` and `#include "biginteger/BigInteger.h"` plus whatever ops/algorithms they need. The Visual Studio solution `bigtst.sln` / `bigtst.vcxproj` is a Windows leftover; CMake is the build everywhere.

## Build & run

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
cd build && ctest --output-on-failure       # unit tests + correctness + mfa_roundtrip
./build/mult_correctness && ./build/div_correctness   # canonical correctness harnesses
./build/bench_vs_gmp                        # needs libgmp (brew install gmp)
```

Requires a C++20 compiler with `unsigned __int128` (GCC ≥ 10 / Clang ≥ 11; uses `std::span` in division paths). One-off probes compile as `c++ -std=c++20 -O2 -Iinclude probe.cpp build/libbigmath.a` — header-only compilation no longer works; the dispatchers are out-of-line in `src/`.

`tests/algotest.cpp` uses deprecated `<strstream>` and may not compile on modern toolchains — prefer `mult_correctness`/`div_correctness` as the canonical correctness harnesses.

CI (`.github/workflows/qa.yaml`) runs `nanoclaw-task qa-agent` on PRs against `origin/main` and uploads `qa-report.md`.

## Architecture

**Number representation.** `BigInteger` (in `include/biginteger/BigInteger.h`) wraps `vector<DataT>` (little-endian limbs) plus a sign bool. `DataT = uint64_t` and limbs hold true 64-bit values by default (`BIGMATH_LIMB_64=1`, `BigInteger::Base()` returns the `Base2_64` sentinel `0`); carry/borrow chains use `__uint128_t`, except the multiplication basecase loops, which go through the flag-chain kernels in `include/biginteger/common/LimbKernels.h` (`AddN`/`SubN`/`Mul1`/`AddMul1`/`MulBasecase`/`SqrBasecase`). `-DBIGMATH_LIMB_64=0` restores the legacy base-2^32-in-64-bit layout. Type aliases (`DataT`, `SizeT`, `BaseT`, `Long`, etc.) come from `include/biginteger/common/Constants.h`.

**Layer split — `ops/` vs `algorithms/`:**
- `include/biginteger/ops/*.h` (+ `src/ops/*.cpp`) — top-level entry points operating on `BigInteger` (sign handling, IO, `Operations.h` umbrella).
- `include/biginteger/algorithms/*.h` — `vector<DataT>`-level routines. Dispatch headers (`algorithms/Multiplication.h`, `algorithms/Division.h`) declare the strategy pick; the out-of-line dispatchers live in `src/algorithms/*.cpp`. Concrete implementations live under `algorithms/multiplication/` and `algorithms/division/`.

**Multiplication dispatch** (`algorithms/Multiplication.h`):
- single-limb operand → `ClassicMultiplication::Multiply(vec, scalar, base)`
- `a.size()+b.size() <= 56`, or `min <= 24` (`BIGMATH_CLASSIC_MIN_LIMB_THRESHOLD`, the Karatsuba leaf cutoff `BIGMATH_KARATSUBA_THRESHOLD`) → classic schoolbook
- `sum >= 640` and `min >= 64` and `la·lb >= 96000` (`NTT_MULTIPLICATION_THRESHOLD` / `NTT_MIN_LIMB_THRESHOLD` / `NTT_AREA_THRESHOLD`; balanced from ~310 limbs per operand) → `NTTMultiplication` (3-prime CRT effectively always: `BIGMATH_NTT_CRT_THRESHOLD = 256`); transform length is `bit_ceil(coeffs)`, except just above a power of two `H` (coeffs ≤ 1.5·H, 2048 ≤ H ≤ 2^20) where the truncated `H + L` scheme (product mod `x^H+1` and mod `x^L−1`, recombined) replaces the `2H` transform
- otherwise → `KaratsubaMultiplication`, which slices operands that are ≥ 2:1 apart into balanced pieces (the area test, not the sum alone, decides NTT because sliced Karatsuba's cost follows `max · K(min)/min`: a short operand keeps it ahead — 48 × 6144 limbs still favors Karatsuba)
- Between `NTT_AREA_THRESHOLD` and `NTT_COLD_AREA_THRESHOLD = 262144` (≈ 512 × 512) the NTT is taken only while the thread pool is warm (`ParallelWarmOrPrewarm()`): with parked workers a small NTT costs its serial time (~52 µs vs 24 µs warm) and loses to Karatsuba. Cold calls run Karatsuba; a burst prewarms the pool and switches to the NTT. Squares: `NTT_SQUARE_THRESHOLD = 352` warm, `NTT_SQUARE_COLD_THRESHOLD = 608` cold.
- Re-swept 2026-11-20 over a (min, max) shape grid; the old `min <= 64 and max >= 10·min → classic` skew guard is retired (`BIGMATH_CLASSIC_SKEW_MIN_LIMB_THRESHOLD = 0`), NTT entry was `sum >= 1280` before.
- `ToomCookMultiplication` exists as an alternate exercised by tests; its dispatch window `[TOOM3_MULTIPLICATION_THRESHOLD, NTT_MULTIPLICATION_THRESHOLD)` was retired 2026-06-12 by keeping the two equal (CRT+NEON beats it everywhere).

**Why no Schönhage-Strassen.** The default NTT is the 3-prime CRT with 32-bit split and a hard 2^26-coefficient transform ceiling (~640M decimal digits, guarded with a throw); the single-prime Goldilocks path (`P = 2^64 - 2^32 + 1`, 16-bit split, overflow-safe to ~2^31 limbs) remains as the `-DBIGMATH_NTT_CRT=0` fallback. Either way NTT covers every practical size. SSA's asymptotic edge is `log log n`; at any size below ~2^40 limbs the constant-factor cost of mod-2^N+1 arithmetic loses to NTT. The ~700+ lines of SSA aren't justified for this codebase — if you need more big-int wins, look at block-recursive Mulders division or Schönhage's subquadratic GCD instead.

**Division dispatch** (`algorithms/Division.h`; canonical thresholds in `include/biginteger/build/DispatchThresholds.h` with value-identical `#ifndef` fallbacks in `algorithms/Division.h`, snapshotted into `const SizeT` in `src/algorithms/Division.cpp`). 64-bit limbs, re-swept 2026-11-21 on a divisor-length × ratio grid; `nb` = divisor limbs, `δ = na − nb` (quotient has `δ+1` limbs):
- `nb < 96` (`DIV_BZ_MIN_B`) → `FastDivision` (Knuth Algorithm D; a raw-buffer loop on `LimbKernel::SubMul1` with the Möller–Granlund 3/2 quotient digit; stack workspace up to 192 limbs)
- `δ ≥ nb` (ratio ≥ 2): `NewtonDivision` when `nb ≥ 1536` (`DIV_NEWTON_MIN_B`), or `nb ≥ 768` at ratio ≥ 3 (`DIV_NEWTON_SKEW_MIN_B` / `_RATIO`); otherwise `BurnikelZieglerDivision`
- `δ < nb` (quotient shorter than the divisor): `QuotientSizedDivision` when `δ ≥ 3072` (`DIV_QSIZED_MIN_QUOTIENT`; it divides the operand tops — `t = δ+4` limbs of b — with Newton, then one `δ × nb` back-multiply); else `NewtonDivision` when `nb ≥ 3072` and ratio ≥ 7/4 (`DIV_NEWTON_NEAR2_MIN_B`); otherwise `BurnikelZieglerDivision`
- `BurnikelZieglerDivision` on 64-bit limbs is an in-place divide-and-conquer in the shape of GMP's `mpn_dcpi1_div_qr`: it works inside the shifted dividend buffer, splits the divisor into floor/ceil halves by pointer offset (odd lengths need no padding), corrects each half-quotient with one product against the divisor's other half, and bottoms out in `FastDivision::DivideNormalized64` at `BIGMATH_DC_DIV_THRESHOLD = 64`. A quotient shorter than the divisor recurses on the divisor's top limbs only, so thin quotients cost what quotient-sized division costs.
- `NewtonDivision`: reciprocal via a guarded precision schedule (`ApproxReciprocal`: exact base of ≤ 32 limbs, each earlier precision `ceil((p+2)/2)`, unpadded step `R·B^e + floor(R·E'/B^(2cur))`; contract `0 ≤ floor(B^(2n)/D) − R ≤ 3`, no refinement iteration), then blockwise chunks of ≤ 2n limbs with wrapped (cyclic) remainders. A leading block of ≤ 32 quotient limbs (`BIGMATH_NEWTON_SCHOOLBOOK_HEAD`; typically the one limb the normalization shift pushes out of a `k·n`-limb dividend) goes to `FastDivision` instead of a full-size chunk. `NewtonDivision::FallbackCount()` counts silent drops to schoolbook; tests assert it stays 0.
- Quotient-only calls (`Divide`, `computeRemainder = false`) with `δ + 4 < nb` skip the rule above: the quotient comes from a `(2δ+4)`-by-`(δ+4)` division of the operands' top limbs (through the dispatcher), exact unless the top remainder is below `B^(δ+2)`, in which case one `δ × nb` product decides between `q'` and `q' − 1`.
- single-limb divisor → `ClassicDivision`
- `KnuthDivision` and `ReciprocalDivision` are alternates used by correctness tests for cross-checking.
- The legacy 32-in-64 layout (`-DBIGMATH_LIMB_64=0`) keeps the pre-2026-11 dispatch bands (`NEWTON_*`, `QSIZED_*`, `BZ_DIVISOR_THRESHOLD`), the vector-based zero-padding BZ (`BIGMATH_BZ_RECURSION_THRESHOLD = 128`) and Newton's exact-doubling reciprocal with its refinement iteration.

When adding a new algorithm, slot the implementation under `algorithms/<op>/<Name>.h`, then update the dispatch in `algorithms/<Op>.h` and add the canonical threshold `#define` to `include/biginteger/build/DispatchThresholds.h` (the `#ifndef` fallbacks in algorithm headers must stay value-identical to it).

**Tests cross-check by construction.** Correctness tests (`mult_correctness.cpp`, `div_correctness.cpp`) run every algorithm against classic/BurnikelZiegler and compare results limb-for-limb via `Compare(...)`; the division test also verifies the identity `q*b + r == a` and `r < b`. When changing an algorithm, run these — passing tests means the new code agrees with the reference path, not just that it ran.

**Calculator app** (`calculator/`): REPL frontend using `ExpressionEvaluator.h` over `BigInteger`. Independent of the test binaries.
