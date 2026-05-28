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
- below `NTT_MULTIPLICATION_THRESHOLD = 32768` → `KaratsubaMultiplication`
- otherwise → `NTTMultiplication`
- `ToomCookMultiplication` exists as an alternate exercised by tests but is not in the default dispatch path.

**Why no Schönhage-Strassen.** NTT here uses the Goldilocks prime `P = 2^64 - 2^32 + 1` with 16-bit input split. Convolution-accumulation overflow only at ~2^31 limbs (≈ 8 GB operands), so NTT covers every practical size. SSA's asymptotic edge is `log log n`; at any size below ~2^40 limbs the constant-factor cost of mod-2^N+1 arithmetic loses to NTT. The ~700+ lines of SSA aren't justified for this codebase — if you need more big-int wins, look at block-recursive Mulders division or Schönhage's subquadratic GCD instead.

**Division dispatch** (`algorithms/Division.h`):
- big-and-skewed (`a.size() > 2048 && a.size() > 3*b.size()`) → `BurnikelZieglerDivision`
- Newton band → `NewtonDivision` — Newton-Raphson reciprocal, O(M(n)). Handles arbitrary `na/nb` via internal blockwise mode (top chunk in [n+1, 2n], slide down by n, thread the remainder). Two dispatch bands: `b ≥ 24576` at ratio ≥ 4/3, or `b ≥ 16384` at ratio ≥ 2. Wins range from 4× (ratio 2) to 20× (ratio 4 at large sizes) vs FD/BZ.
- otherwise multi-limb → `FastDivision` (Knuth Algorithm D variant)
- single-limb divisor → `ClassicDivision`
- `KnuthDivision` and `ReciprocalDivision` are alternates used by correctness tests for cross-checking.

When adding a new algorithm, slot the implementation under `algorithms/<op>/<Name>.h`, then update the dispatch in `algorithms/<Op>.h` — the thresholds there are the only place size cutoffs live.

**Tests cross-check by construction.** Correctness tests (`mult_correctness.cpp`, `div_correctness.cpp`) run every algorithm against classic/BurnikelZiegler and compare results limb-for-limb via `Compare(...)`; the division test also verifies the identity `q*b + r == a` and `r < b`. When changing an algorithm, run these — passing tests means the new code agrees with the reference path, not just that it ran.

**Calculator app** (`calculator/`): REPL frontend using `ExpressionEvaluator.h` over `BigInteger`. Independent of the test binaries.
