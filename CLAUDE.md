# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Header-only C++ BigInteger library in namespace `BigMath`. All implementation lives under `biginteger/` as `.h` files; consumers `#include "biginteger/BigInteger.h"` plus whatever ops/algorithms they need. No build system on macOS/Linux — Visual Studio solution `bigtst.sln` / `bigtst.vcxproj` exists for Windows, everything else compiles standalone per-test.

## Build & run

Each test in `tests/` is a self-contained translation unit with its own `main()`. Compile with C++20 (uses `std::span` in division paths; also uses `unsigned __int128` everywhere — GCC/Clang required):

```
c++ -std=c++20 -O2 tests/mult_correctness.cpp -o mult_correctness && ./mult_correctness
c++ -std=c++20 -O2 tests/div_correctness.cpp -o div_correctness && ./div_correctness
c++ -std=c++20 -O2 tests/multperf_simple.cpp -o multperf_simple && ./multperf_simple
c++ -std=c++20 -O2 tests/divperf_simple.cpp -o divperf_simple && ./divperf_simple
c++ -std=c++20 -O2 tests/stress/stresstest.cpp -o stresstest   # needs input.txt/answer.txt
c++ -std=c++20 -O2 calculator/calculator.cpp -o calculator && ./calculator
```

`tests/algotest.cpp` uses deprecated `<strstream>` and may not compile on modern toolchains — prefer `mult_correctness`/`div_correctness` as the canonical correctness harnesses.

CI (`.github/workflows/qa.yaml`) runs `nanoclaw-task qa-agent` on PRs against `origin/main` and uploads `qa-report.md`. No CMake/Make to run locally.

## Architecture

**Number representation.** `BigInteger` (in `biginteger/BigInteger.h`) wraps `vector<DataT>` (little-endian limbs) plus a sign bool. `DataT = uint64_t` but limbs hold base-2^32 values (`BigInteger::Base()` returns `Base2_32`); the upper 32 bits are headroom for carries during arithmetic. Type aliases (`DataT`, `SizeT`, `BaseT`, `Long`, etc.) come from `biginteger/common/Constants.h`.

**Layer split — `ops/` vs `algorithms/`:**
- `biginteger/ops/*.h` — top-level entry points operating on `BigInteger` (sign handling, IO, `Operations.h` umbrella).
- `biginteger/algorithms/*.h` — `vector<DataT>`-level routines. Dispatch headers (`algorithms/Multiplication.h`, `algorithms/Division.h`) pick a strategy by operand size; concrete implementations live under `algorithms/multiplication/` and `algorithms/division/`.

**Multiplication dispatch** (`algorithms/Multiplication.h`):
- single-limb operand → `ClassicMultiplication::Multiply(vec, scalar, base)`
- `a.size()+b.size() <= 96` or `min(a,b) <= 32` → classic schoolbook
- below `NTT_MULTIPLICATION_THRESHOLD = 1280` total limbs → `KaratsubaMultiplication` (dropped from 5120 on 2026-06-12 after the NEON Shoup butterflies made CRT NTT ~2.5× faster in the small-mid band)
- otherwise → `NTTMultiplication` (3-prime CRT effectively always: `BIGMATH_NTT_CRT_THRESHOLD = 256`)
- `ToomCookMultiplication` exists as an alternate exercised by tests; its dispatch window was retired 2026-06-12 (CRT+NEON beats it everywhere).

**Why no Schönhage-Strassen.** NTT here uses the Goldilocks prime `P = 2^64 - 2^32 + 1` with 16-bit input split. Convolution-accumulation overflow only at ~2^31 limbs (≈ 8 GB operands), so NTT covers every practical size. SSA's asymptotic edge is `log log n`; at any size below ~2^40 limbs the constant-factor cost of mod-2^N+1 arithmetic loses to NTT. The ~700+ lines of SSA aren't justified for this codebase — if you need more big-int wins, look at block-recursive Mulders division or Schönhage's subquadratic GCD instead.

**Division dispatch** (`algorithms/Division.h`, thresholds defined in `src/algorithms/Division.cpp`):
- `NewtonDivision` (Newton-Raphson reciprocal, O(M(n)); handles arbitrary `na/nb` via blockwise mode — top chunk in [n+1, 2n], slide down by n, thread the remainder) when any skew band holds:
  - `b ≥ 640` at ratio ≥ 8 (`NEWTON_HIGH_SKEW` 8/1), or
  - `b ≥ 1024` at ratio ≥ 7/2 (`NEWTON_RATIO35`), or
  - `b ≥ 1280` at ratio ≥ 14/5 (`NEWTON_SKEW` — 2.8 keeps the ratio-3.0000±1-limb knife on the Newton side), or
  - `b ≥ 1792` at ratio ≥ 5/2 (`NEWTON_MID`), or
  - `b ≥ 4096` at ratio ≥ 8/5 (`NEWTON_RATIO2`), or
  - `b ≥ 24576` at ratio ≥ 4/3 (`NEWTON_BALANCED` — near-balanced band; ratio lowered from 2/1 and floor from 98304 on 2026-06-11).
  (Frontier re-swept 2026-06-12 after the BZ odd-size padding fix — see below — changed every Newton/BZ crossover; `include/biginteger/build/DispatchThresholds.h` is the canonical place these live — its `#define`s win over the `#ifndef` fallbacks in `algorithms/*.h`.)
- `QuotientSizedDivision` when (`b ≥ 24576`, `a ≥ b + 64`, ratio < 4/3) OR (thin-quotient: `b ≥ 8192`, `64 ≤ Δ ≤ b/8`): divides the operand TOPS (`t = Δ+4` limbs of b, `Δ+t` of a) for the (Δ+1)-limb quotient, then one Δ×nb back-multiply for the remainder — cost scales with the quotient, not the divisor.
- else `BurnikelZieglerDivision` for power-of-two base when `b > 512` and the BZ band fits (near-balanced `b ≥ 1024, b+32 ≤ a ≤ 3b`, or big-and-skewed `a > 2048 && a > 3b`). BZ pads both operands with bottom zero limbs so the divisor size divides by 2 at every recursion level (quotient unchanged, remainder's bottom pad limbs are zero); before 2026-06-12 any odd size on the way down fell back to FastDivision at full size — the root cause of the "2^k+1 family" pathology and 1.3-10× losses on odd-limb-count divisors at ratio < 2.8.
- otherwise multi-limb → `FastDivision` (Knuth Algorithm D variant)
- single-limb divisor → `ClassicDivision`
- `KnuthDivision` and `ReciprocalDivision` are alternates used by correctness tests for cross-checking.

The balanced band exists because BZ's recursive 2n/n halving lands intermediate NTT multiplies just over power-of-2 transform-length boundaries for non-power-of-2 divisor sizes, blowing up 5–60× vs Newton (worst at `n = 2^k+1`); Newton pads once and stays flat. The former ratio ∈ (1, 4/3) residual is closed by `QuotientSizedDivision` (2026-06-11): the `2^k+1`-family blowup shapes went from seconds to tens of ms, and it beats BZ even at BZ's exact-power-of-2 best case.

When adding a new algorithm, slot the implementation under `algorithms/<op>/<Name>.h`, then update the dispatch in `algorithms/<Op>.h` — the thresholds there are the only place size cutoffs live.

**Tests cross-check by construction.** Correctness tests (`mult_correctness.cpp`, `div_correctness.cpp`) run every algorithm against classic/BurnikelZiegler and compare results limb-for-limb via `Compare(...)`; the division test also verifies the identity `q*b + r == a` and `r < b`. When changing an algorithm, run these — passing tests means the new code agrees with the reference path, not just that it ran.

**Calculator app** (`calculator/`): REPL frontend using `ExpressionEvaluator.h` over `BigInteger`. Independent of the test binaries.
