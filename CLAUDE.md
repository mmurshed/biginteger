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

**Number representation.** `BigInteger` (in `include/biginteger/BigInteger.h`) wraps `vector<DataT>` (little-endian limbs) plus a sign bool. `DataT = uint64_t` and limbs hold true 64-bit values by default (`BIGMATH_LIMB_64=1`, `BigInteger::Base()` returns the `Base2_64` sentinel `0`); carry/borrow chains use `__uint128_t`. `-DBIGMATH_LIMB_64=0` restores the legacy base-2^32-in-64-bit layout. Type aliases (`DataT`, `SizeT`, `BaseT`, `Long`, etc.) come from `include/biginteger/common/Constants.h`.

**Layer split — `ops/` vs `algorithms/`:**
- `include/biginteger/ops/*.h` (+ `src/ops/*.cpp`) — top-level entry points operating on `BigInteger` (sign handling, IO, `Operations.h` umbrella).
- `include/biginteger/algorithms/*.h` — `vector<DataT>`-level routines. Dispatch headers (`algorithms/Multiplication.h`, `algorithms/Division.h`) declare the strategy pick; the out-of-line dispatchers live in `src/algorithms/*.cpp`. Concrete implementations live under `algorithms/multiplication/` and `algorithms/division/`.

**Multiplication dispatch** (`algorithms/Multiplication.h`):
- single-limb operand → `ClassicMultiplication::Multiply(vec, scalar, base)`
- `a.size()+b.size() <= 96`, or heavy skew (`min <= 64` and `max >= 10·min`) → classic schoolbook
- below `NTT_MULTIPLICATION_THRESHOLD = 1280` total limbs → `KaratsubaMultiplication` (dropped from 5120 on 2026-06-12 after the NEON Shoup butterflies made CRT NTT ~2.5× faster in the small-mid band)
- otherwise → `NTTMultiplication` (3-prime CRT effectively always: `BIGMATH_NTT_CRT_THRESHOLD = 256`)
- `ToomCookMultiplication` exists as an alternate exercised by tests; its dispatch window was retired 2026-06-12 (CRT+NEON beats it everywhere).

**Why no Schönhage-Strassen.** The default NTT is the 3-prime CRT with 32-bit split and a hard 2^26-coefficient transform ceiling (~640M decimal digits, guarded with a throw); the single-prime Goldilocks path (`P = 2^64 - 2^32 + 1`, 16-bit split, overflow-safe to ~2^31 limbs) remains as the `-DBIGMATH_NTT_CRT=0` fallback. Either way NTT covers every practical size. SSA's asymptotic edge is `log log n`; at any size below ~2^40 limbs the constant-factor cost of mod-2^N+1 arithmetic loses to NTT. The ~700+ lines of SSA aren't justified for this codebase — if you need more big-int wins, look at block-recursive Mulders division or Schönhage's subquadratic GCD instead.

**Division dispatch** (`algorithms/Division.h`; threshold macros live in `include/biginteger/build/DispatchThresholds.h` with `#ifndef` fallbacks in `algorithms/Division.h`, snapshotted into `const SizeT` in `src/algorithms/Division.cpp`):
- `NewtonDivision` (Newton-Raphson reciprocal, O(M(n)); handles arbitrary `na/nb` via blockwise mode — top chunk in [n+1, 2n], slide down by n, thread the remainder) when any skew band holds:
  - `b ≥ 896` at ratio ≥ 8 (`NEWTON_HIGH_SKEW` 8/1), or
  - `b ≥ 1280` at ratio ≥ 7/2 (`NEWTON_RATIO35` — 7/2 keeps the ratio-4.0000±1-limb knife on the Newton side), or
  - `b ≥ 1792` at ratio ≥ 14/5 (`NEWTON_SKEW` — 2.8 likewise for ratio 3.0000), or
  - `b ≥ 2560` at ratio ≥ 5/2 (`NEWTON_MID`), or
  - `b ≥ 4096` at ratio ≥ 2 (`NEWTON_RATIO20`), or
  - `b ≥ 8192` at ratio ≥ 8/5 (`NEWTON_RATIO2`), or
  - `b ≥ 131072` at ratio ≥ 4/3 (`NEWTON_BALANCED` — floor raised from 24576 on 2026-06-12: padded BZ with the 128-limb basecase wins ratio 1.4-1.5 through ~98k limbs).
  (Frontier re-swept twice 2026-06-12: once after the BZ odd-size padding fix, again after the BZ basecase retune 512→128 — both changed every Newton/BZ crossover; `include/biginteger/build/DispatchThresholds.h` is the canonical place these live — its `#define`s win over the `#ifndef` fallbacks in `algorithms/*.h`.)
- `QuotientSizedDivision` when (`b ≥ 24576` (`QSIZED_MAIN_B`, decoupled from the balanced floor), `a ≥ b + 64`, ratio < 4/3) OR (thin-quotient: `b ≥ 8192`, `64 ≤ Δ ≤ b/8`): divides the operand TOPS (`t = Δ+4` limbs of b, `Δ+t` of a) for the (Δ+1)-limb quotient, then one Δ×nb back-multiply for the remainder — cost scales with the quotient, not the divisor.
- else `BurnikelZieglerDivision` for power-of-two base when `b > 512` and the BZ band fits (near-balanced `b ≥ 1024, b+32 ≤ a ≤ 3b`, or big-and-skewed `a > 2048 && a > 3b`). BZ pads both operands with bottom zero limbs so the divisor size divides by 2 at every recursion level (quotient unchanged, remainder's bottom pad limbs are zero); before 2026-06-12 any odd size on the way down fell back to FastDivision at full size — the root cause of the "2^k+1 family" pathology and 1.3-10× losses on odd-limb-count divisors at ratio < 2.8. Recursion basecase `BIGMATH_BZ_RECURSION_THRESHOLD = 128` (was 512; the retune is worth another 1.2-1.5× across the BZ band — Karatsuba-backed multiplies beat 512-limb Knuth-D basecase calls).
- otherwise multi-limb → `FastDivision` (Knuth Algorithm D variant)
- single-limb divisor → `ClassicDivision`
- `KnuthDivision` and `ReciprocalDivision` are alternates used by correctness tests for cross-checking.

The balanced band exists because BZ's recursive 2n/n halving lands intermediate NTT multiplies just over power-of-2 transform-length boundaries for non-power-of-2 divisor sizes, blowing up 5–60× vs Newton (worst at `n = 2^k+1`); Newton pads once and stays flat. The former ratio ∈ (1, 4/3) residual is closed by `QuotientSizedDivision` (2026-06-11): the `2^k+1`-family blowup shapes went from seconds to tens of ms, and it beats BZ even at BZ's exact-power-of-2 best case.

When adding a new algorithm, slot the implementation under `algorithms/<op>/<Name>.h`, then update the dispatch in `algorithms/<Op>.h` and add the canonical threshold `#define` to `include/biginteger/build/DispatchThresholds.h` (the `#ifndef` fallbacks in algorithm headers must stay value-identical to it).

**Tests cross-check by construction.** Correctness tests (`mult_correctness.cpp`, `div_correctness.cpp`) run every algorithm against classic/BurnikelZiegler and compare results limb-for-limb via `Compare(...)`; the division test also verifies the identity `q*b + r == a` and `r < b`. When changing an algorithm, run these — passing tests means the new code agrees with the reference path, not just that it ran.

**Calculator app** (`calculator/`): REPL frontend using `ExpressionEvaluator.h` over `BigInteger`. Independent of the test binaries.
