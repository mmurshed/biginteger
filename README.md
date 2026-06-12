# BigMath

Arbitrary-precision integer library in C++20. Header-light with a thin static-library shell. 64-bit limbs, 3-prime CRT NTT multiplication (multithreaded by default), Newton–Raphson division with cached reciprocals, divide-and-conquer base conversion. Includes a REPL calculator built on top.

Target: portable C++ (with optional NEON intrinsics on aarch64, scalar fallbacks everywhere else) that reaches GMP parity or better across the large-operand bands — no hand-written assembly.

## Performance

BigMath vs GMP 6.3, Apple M1 Max, paired same-run measurements (2026-06-12). Below the dashed line BigMath is faster than GMP:

![BigMath vs GMP ratio by operand size](docs/images/bigmath_vs_gmp.png)

- **Multiplication beats GMP at every balanced size from 500k digits up** — peaking at 3.2× faster (10M×10M: 65 vs 206 ms) and staying at or below GMP through 200M digits (warm steady state) — and at **every measured skewed shape ≥500k×50k** (0.40–0.77×).
- **Division beats GMP from 20M-digit dividends (0.45–0.82×)** and reaches parity from 5M (200M÷40M: 3.0 s vs GMP's 4.5 s).
- **Decimal I/O beats GMP from 500k digits in both directions** after the parallel D&C fan-outs (PRs #118/#119): warm ToString 0.50–0.69× (1M digits: 28 vs 49 ms), parse 0.38–0.71× (20M digits: 304 vs 799 ms, 2.6× faster).
- Sub-NTT sizes (≲ 13k digits) remain 2–3× behind GMP's hand-tuned basecase — the documented portable-C++ wall.

Two 2026-06 optimization runs produced the current margins — PRs #82–#99 (wraparound Newton division, dispatch band retunes, quotient-sized division, NEON Shoup NTT butterflies) and PRs #107–#112 (fused-MFA pass fusion, row-chunked stage parallelism, MFA gate 2^24→2^20, cyclic products on the fused pipeline, on-the-fly operand packing). PR #116 then fixed a long-standing Burnikel–Ziegler blind spot: odd divisor limb counts silently fell back to quadratic Knuth D, costing 1.3–10× on ~half of real division shapes in the 1k–25k-limb band (12289-limb divisor at ratio 1.5: 95 → 8 ms) — the prior benchmark tables, whose digit-derived shapes mostly landed on even sizes, never sampled it.

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
- **Squaring:** Specialized Classic / Karatsuba / NTT squarers (1.4–1.6× over `Multiply(a,a)`).
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

## Performance

Apple M1 Max, vs GMP 6.3.0, `-O3 -march=native`, full default stack (`BIGMATH_LIMB_64=1` + `BIGMATH_NTT_CRT=1` + `BIGMATH_USE_THREADS=1`, 8-thread pool). Canonical `bench_vs_gmp` run, 2026-06-12 evening (post PR #107–#112):

| operation | size | BigMath | GMP | ratio |
|---|---|---:|---:|---:|
| mul | 100 000 × 100 000 | 0.66 ms | 0.61 ms | 1.08× |
| mul | **1 000 000 × 1 000 000** | **7.8 ms** | **9.6 ms** | **0.81×** ← BigMath faster |
| mul | **5 000 000 × 5 000 000** | **32.6 ms** | **64.4 ms** | **0.51×** ← BigMath 2× faster |
| mul | **10 000 000 × 10 000 000** | **65 ms** | **206 ms** | **0.32×** ← BigMath 3.2× faster |
| mul | **20 000 000 × 20 000 000** | **129 ms** | **276 ms** | **0.47×** ← BigMath 2.1× faster |
| mul | **200 000 000 × 200 000 000** | **3 453 ms** | **3 307 ms** | **1.04×** ← 0.97× warm |
| mul (skewed) | **1 000 000 / 100 000** | **2.5 ms** | **4.6 ms** | **0.54×** ← BigMath faster |
| mul (skewed) | **20 000 000 / 2 000 000** | **94 ms** | **162 ms** | **0.58×** ← BigMath faster |
| mul (skewed) | **50 000 000 / 5 000 000** | **267 ms** | **671 ms** | **0.40×** ← BigMath 2.5× faster |
| div (skewed) | 500 000 / 100 000 | 7.6 ms | 4.6 ms | 1.65× |
| div (skewed) | **5 000 000 / 1 000 000** | **70 ms** | **70 ms** | **1.01×** ← parity |
| div (skewed) | **50 000 000 / 10 000 000** | **587 ms** | **1 307 ms** | **0.45×** ← BigMath 2.2× faster |
| div (skewed) | **200 000 000 / 40 000 000** | **3 010 ms** | **4 492 ms** | **0.67×** ← BigMath faster |
| parse | 1 000 000 digits | 32 ms | 21 ms | 1.58× |
| parse | **20 000 000 digits** | **814 ms** | **804 ms** | **1.01×** ← parity |
| ToString | 100 000 digits | 3.0 ms | 2.4 ms | 1.25× |
| ToString | **1 000 000 digits** | **28 ms** | **49 ms** | **0.57×** ← BigMath faster |
| ToString | **10 000 000 digits** | **481 ms (warm)** | **908 ms** | **0.53×** ← BigMath faster |

**BigMath beats GMP on balanced multiplication at every size from 500k digits up** — the former ≥50M-digit losses ("GMP SSA recovers") closed once the fused-MFA pipeline landed: 50M–200M digits run at 0.90–0.97× warm. **Skewed division flipped from a 1.9–2.8× loss to a 0.45–0.67× win at ≥20M-digit dividends** as Newton inherits the fused multiplies and runs its wrap-around products on the same pipeline. **ToString flipped from the largest remaining gap to a 0.50–0.69× win at ≥500k digits** (PR #118: the D&C formatter's subtrees are fixed-width fields at precomputable offsets — one ParallelDo over 8 of them; `tostring_chain_plan.md` records the outcome and the deferred cold-chain lever). The 10k–2M-digit division plan (`smallskew_div_plan.md`) is executed: its sweep surfaced the Burnikel–Ziegler odd-size fallback fixed in PR #116, and what remains of the sub-50k-digit division gap is the documented portable-C++ basecase wall, accepted there. See [BENCHMARK.md](BENCHMARK.md) for full tables, warm/cold methodology, and per-PR history.

Opt-out flags (`-DBIGMATH_USE_THREADS=0` / `-DBIGMATH_NTT_CRT=0` / `-DBIGMATH_LIMB_64=0`) revert any subset of the defaults — useful for embedded targets, header-only-strict consumers, or A/B comparison.

## License

See [LICENSE](LICENSE).

## Author

S. M. Mahbub Murshed (murshed@gmail.com)
