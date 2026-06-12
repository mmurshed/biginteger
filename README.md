# BigMath

Arbitrary-precision integer library in C++20. Header-light with a thin static-library shell. 64-bit limbs, 3-prime CRT NTT multiplication (multithreaded by default), Newton–Raphson division with cached reciprocals, divide-and-conquer base conversion. Includes a REPL calculator built on top.

Target: portable C++ (with optional NEON intrinsics on aarch64, scalar fallbacks everywhere else) that reaches GMP parity or better across the large-operand bands — no hand-written assembly.

## Performance

BigMath vs GMP 6.3, Apple M1 Max, paired same-run measurements (canonical run 2026-06-12 at v13.0). Below the dashed line BigMath is faster than GMP:

![BigMath vs GMP ratio by operand size](docs/images/bigmath_vs_gmp.png)

- **Multiplication beats GMP at every balanced size from 500k digits up** — peaking at 3.2× faster (10M×10M: 64 vs 209 ms) and staying at or below GMP through 200M digits warm (0.89–0.97×) — and at **every measured skewed shape ≥500k×50k** (0.40–0.73×).
- **Division beats GMP from 20M-digit dividends (0.45–0.82×)** and reaches parity at 5M (200M÷40M: 3.1 s vs GMP's 4.5 s).
- **Decimal I/O beats GMP from 500k digits in both directions**: parse 0.38–0.69× (20M digits: 306 vs 815 ms, 2.7× faster), warm ToString 0.49–0.62× (1M digits: 28 vs 50 ms).
- Sub-NTT sizes (≲ 13k digits) remain 2–3× behind GMP's hand-tuned basecase — the documented portable-C++ wall.

The margins come from the 2026 optimization campaign — NEON Shoup NTT butterflies, fused-MFA transforms, the wraparound-Newton division family, the Burnikel–Ziegler odd-size padding fix, parallel decimal I/O fan-outs, and CRT-NTT squaring — followed by a full-codebase correctness audit. The condensed history with per-step wins is in [BENCHMARK.md](BENCHMARK.md#optimization-history-condensed); releases in [CHANGELOG.md](CHANGELOG.md).

![Session before/after](docs/images/session_2026_06_11.png)

Full tables, methodology, and per-PR history: [BENCHMARK.md](BENCHMARK.md).

---

## Features

- **64-bit limb representation** (`BIGMATH_LIMB_64=1`, default). `DataT` stores true 64-bit values; every carry/borrow chain uses `__uint128_t` accumulators. Halves loop iteration counts in scalar paths vs the legacy 32-bit-in-64-bit layout. Opt out via `-DBIGMATH_LIMB_64=0`.
- **Multiplication:** Classic schoolbook → Karatsuba (64-bit-hybrid leaf) → 3-prime CRT NTT (2013265921 / 469762049 / 1811939329, 32-bit coefficient splitting) from 1280 total limbs (~13k digits). Toom-3 and the single-prime Goldilocks NTT remain as cross-check alternates; the NEON-accelerated CRT path beats both at every measured size.
- **Radix-4 + radix-8 fused NTT butterflies** (PRs #59, #60): 3× fewer memory passes vs radix-2; ~1.6× wall-clock at ≥2M limbs.
- **NEON Shoup butterflies on aarch64** (PRs #96–#99, `BIGMATH_NEON=1` default on Apple Silicon): the radix-8/4/2 layers and 1/n scaling run 4 lanes wide using Shoup's precomputed-reciprocal multiplication — `(w·x) mod P` from one widening multiply, one low multiply, one conditional subtract, with interleaved `(twiddle, w′)` tables keeping the strided gathers at one cache line per pair. Bit-exact with the scalar path; ~4× per butterfly pass, −26…40% end-to-end on NTT-bound ops. Scalar fallback everywhere else.
- **MFA / Bailey 6-step CRT NTT** (`BIGMATH_NTT_MFA_THRESHOLD=2^20`, default). Large CRT transforms switch to a cache-friendly matrix Fourier layout with fully fused stages: the transpose, the operand packing (straight from the limb arrays), and the pointwise product all ride the tiled stage passes, and the forward/inverse stages run as row-chunked `ParallelDo(6)` work units. The gate sat at `2^24` until the fused stages flipped the break-even (PR #109); cyclic wrap-around products use the same pipeline up to the CRT ceiling `2^26` (PR #110).
- **Multithreaded NTT** (`BIGMATH_USE_THREADS=1`, default). Small thread pool (size `min(hw_concurrency, BIGMATH_MAX_THREADS=8)`) parallelizes the CRT path: in the fused MFA window each transform stage runs as six (operand×prime or prime×row-range) work units; the non-MFA path batches 6 forwards + 3 inverses. 2.3-3.4× speedup on large mul / skewed div / parse. Opt out via `-DBIGMATH_USE_THREADS=0` to drop pthread linkage.
- **Division:** Classic short division → Knuth Algorithm D (`FastDivision`, Möller-Granlund qhat, bit-shift normalization) → Burnikel–Ziegler (small near-balanced) → Newton–Raphson with cached reciprocals and wrap-around cyclic products (half-length transforms for the remainder, quotient-estimate, and reciprocal-iteration steps — GMP `mu_div`/`invertappr` style) → quotient-sized division for short-quotient shapes (cost scales with the quotient, not the divisor). Burnikel–Ziegler pads operands with bottom zero limbs so its recursion never hits an odd split (PR #116 — before that, odd divisor limb counts fell back to Knuth D at full size, 1.3–10× slower on ~half of real shapes). Newton dispatch bands form a ratio/size frontier (8/1 @ 896 limbs, 7/2 @ 1280, 14/5 @ 1792, 5/2 @ 2560, 2/1 @ 4096, 8/5 @ 8192, 4/3 @ 131072) with fractional ratios chosen to keep integer digit-ratio knife-edges (2.0000, 3.0000, 4.0000 ± 1 limb) off the band boundaries. Identity `q·b + r == a` is cross-checked in `tests/div_correctness.cpp`.
- **Squaring:** Specialized Classic / Karatsuba squarers below 640 limbs; above, the CRT NTT self-multiply (skips the duplicate forward transforms — 2.5–5.8× over the former Goldilocks `NTTSquare`, which remains the `-DBIGMATH_NTT_CRT=0` fallback).
- **String I/O:** Linear chunked parser/formatter for small inputs, divide-and-conquer with cached Newton reciprocals at scale; from 100k digits both directions fan their D&C subtrees out over the thread pool (the subtree outputs are fixed-width fields at precomputable offsets). Asymptotic `O(M(L) · log L)` both directions.
- **BigDecimal:** Java-style fixed-point decimal (unscaled BigInteger + int scale) with exact +, −, \*; rounded division taking 8 rounding modes; parse/format covering plain and scientific notation.
- **Calculator REPL:** Variables, hex/bin/dec output, multi-line continuation, comments, `:help :quit :vars :reset :base :digits :time :load :save` directives. Built as a separate executable.
- **Thread-safe by construction** for concurrent use of distinct objects from distinct threads. See [docs/THREAD_SAFETY.md](docs/THREAD_SAFETY.md).

## Layout

```
include/biginteger/  public BigInteger headers (algorithms, ops, common)
include/bigdecimal/  public BigDecimal headers
src/                 non-inline BigInteger implementations
bigdecimal/          BigDecimal implementation
tests/               correctness, performance, unit tests
calculator/          REPL frontend
docs/                technical references (see below)
```

## Requirements

- CMake ≥ 3.16
- C++20 compiler with `__uint128_t` + `__builtin_*_overflow` — GCC ≥ 10, Clang ≥ 11, or Apple Clang. MSVC is not supported.
- (Optional) `libgmp` for the GMP comparison benchmark.

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Outputs:
- `build/libbigmath.a` — static library
- `build/calculator` — REPL
- `build/{mult_correctness, div_correctness, unit_tests}` — test binaries
- `build/{multperf_simple, divperf_simple, bench_vs_gmp}` — performance harnesses

See [BUILDING.md](BUILDING.md) for shared-library mode, threshold overrides, and installation.

## Test

```sh
cd build && ctest --output-on-failure
```

`mult_correctness` is fast; `div_correctness` runs a broad matrix and takes a few minutes; `unit_tests` is sub-second.

The standalone perf harnesses `multperf_simple` and `divperf_simple` (BigMath-internal, no GMP) take an optional run-count `k` (default 3); each size is timed `k` times and the mean reported:

```sh
./build/multperf_simple 5   # average of 5 runs per size
./build/divperf_simple 5
```

See [BENCHMARK.md](BENCHMARK.md#in-tree-simple-harnesses-k-averaged-no-gmp) for the latest k=5 numbers.

## Use

```cpp
#include "biginteger/BigInteger.h"
#include "biginteger/common/Builder.h"
#include "biginteger/common/Parser.h"
#include "biginteger/ops/Operations.h"

using namespace BigMath;

BigInteger a = BigIntegerBuilder::From("123456789012345678901234567890");
BigInteger b = BigIntegerBuilder::From("987654321098765432109876543210");
BigInteger c = a * b;
std::cout << ToString(c) << '\n';
```

CMake consumers:

```cmake
find_package(bigmath REQUIRED)
target_link_libraries(my_app PRIVATE bigmath::bigmath)
```

## Calculator

```sh
./build/calculator
> x = 2^256 - 1
> x * x
> :base 16
> x
> :quit
```

Operators: `+ - * / % ^` (with `^` right-associative). Literals: decimal, `0x…` hex, `0b…` binary, underscore separators (`1_000_000`). `#` to end of line is a comment. Trailing `\` continues a line.

## Documentation

- [docs/BASE.md](docs/BASE.md) — number representation, 64-bit limbs, why little-endian
- [docs/MULTIPLICATION.md](docs/MULTIPLICATION.md) — Classic / Karatsuba / Toom-3 / NTT (Goldilocks + multi-prime CRT) and their tradeoffs
- [docs/DIVISION.md](docs/DIVISION.md) — Classic / Fast (Knuth D + Möller-Granlund qhat) / Burnikel–Ziegler / Newton / Reciprocal-cached
- [docs/STRING_CONVERSION.md](docs/STRING_CONVERSION.md) — chunked decimal I/O, D&C parse/format, Newton-divider chain
- [docs/BIGDECIMAL.md](docs/BIGDECIMAL.md) — fixed-point decimal model, rounding modes, performance
- [docs/THREAD_SAFETY.md](docs/THREAD_SAFETY.md) — concurrency model, opt-in internal parallelism
- [BUILDING.md](BUILDING.md) — CMake build, install, build flags, threshold tuning

Each doc covers algorithms, dispatch, benchmark numbers vs GMP, optimizations that landed, and approaches that were tried and rejected with reasons.

## Performance numbers

Apple M1 Max, vs GMP 6.3.0, `-O3 -march=native`, full default stack (`BIGMATH_LIMB_64=1` + `BIGMATH_NTT_CRT=1` + `BIGMATH_USE_THREADS=1`, 8-thread pool). Canonical `bench_vs_gmp` run, 2026-06-12 at v13.0; rows marked *warm* are steady-state best-of-3 (see [BENCHMARK.md](BENCHMARK.md) methodology):

| operation | size (digits) | BigMath | GMP | ratio |
|---|---|---:|---:|---:|
| mul | 100 000 × 100 000 | 0.80 ms | 0.62 ms | 1.29× |
| mul | **1 000 000 × 1 000 000** | **7.9 ms** | **9.3 ms** | **0.85×** ← BigMath faster |
| mul | **5 000 000 × 5 000 000** | **32.9 ms** | **63.5 ms** | **0.52×** ← 1.9× faster |
| mul | **10 000 000 × 10 000 000** | **64 ms** | **209 ms** | **0.31×** ← 3.2× faster |
| mul | **20 000 000 × 20 000 000** | **129 ms** | **275 ms** | **0.47×** ← 2.1× faster |
| mul | **200 000 000 × 200 000 000** | **2 809 ms (warm)** | **2 907 ms** | **0.97×** |
| mul (skewed) | **1 000 000 / 100 000** | **2.8 ms** | **4.6 ms** | **0.61×** ← BigMath faster |
| mul (skewed) | **20 000 000 / 2 000 000** | **94 ms** | **161 ms** | **0.58×** |
| mul (skewed) | **50 000 000 / 5 000 000** | **272 ms** | **677 ms** | **0.40×** ← 2.5× faster |
| div (skewed) | 500 000 / 100 000 | 7.7 ms | 4.8 ms | 1.59× |
| div (skewed) | **5 000 000 / 1 000 000** | **69 ms** | **69 ms** | **1.00×** ← parity |
| div (skewed) | **50 000 000 / 10 000 000** | **590 ms** | **1 318 ms** | **0.45×** ← 2.2× faster |
| div (skewed) | **200 000 000 / 40 000 000** | **3 119 ms** | **4 544 ms** | **0.69×** ← BigMath faster |
| parse | **1 000 000** | **12.5 ms** | **20.9 ms** | **0.60×** ← BigMath faster |
| parse | **20 000 000** | **306 ms** | **815 ms** | **0.38×** ← 2.7× faster |
| ToString | 100 000 | 3.4 ms (warm) | 2.3 ms | 1.46× |
| ToString | **1 000 000** | **27.7 ms (warm)** | **50.3 ms** | **0.55×** ← BigMath faster |
| ToString | **10 000 000** | **488 ms (warm)** | **901 ms** | **0.54×** |

What remains behind GMP — division below ~200k-digit dividends and everything below ~13k digits — is the hand-tuned-assembly basecase band, documented as the accepted portable-C++ wall in the subsystem docs. See [BENCHMARK.md](BENCHMARK.md) for full tables, warm/cold methodology, and the condensed optimization history.

Opt-out flags (`-DBIGMATH_USE_THREADS=0` / `-DBIGMATH_NTT_CRT=0` / `-DBIGMATH_LIMB_64=0`) revert any subset of the defaults — useful for embedded targets, header-only-strict consumers, or A/B comparison.

## License

See [LICENSE](LICENSE).

## Author

S. M. Mahbub Murshed (murshed@gmail.com)
