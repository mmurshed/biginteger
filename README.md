# BigMath

Arbitrary-precision integer library in C++20. Header-light with a thin static-library shell. 64-bit limbs, 3-prime CRT NTT multiplication (multithreaded by default), Newton–Raphson division with cached reciprocals, divide-and-conquer base conversion. Includes a REPL calculator built on top.

Target: portable C++ (with optional NEON intrinsics on aarch64, scalar fallbacks everywhere else) that reaches GMP parity or better across the large-operand bands — no hand-written assembly.

## Performance

BigMath vs GMP 6.3, Apple M1 Max, paired same-run measurements (canonical runs: multiplication 2026-11-20, division and decimal I/O 2026-11-21). GMP is single-threaded; BigMath uses its 8-thread pool on NTT-sized work.

- **Multiplication beats GMP from 10k digits up** (was 500k before the 2026-11 pass) — 0.21–0.63× from 50k to 20M digits, peaking at 4.8× faster (10M×10M: 43 vs 207 ms), and 0.63–0.72× warm through 200M digits. Skewed shapes: 0.27–0.51×.
- **Division beats GMP from 200k-digit dividends (0.19–0.67×)** — 50M÷10M digits: 244 ms vs GMP's 1.31 s; 200M÷40M: 1.45 s vs 4.47 s. It was 500k before the 2026-11 division pass, which sped up the 12 published skewed division cases by 1.14–1.58×. The historical table compares BigMath quotient-only with GMP quotient and remainder.
- **Decimal I/O**: parse beats GMP from 100k digits (0.23–0.54×; 20M digits: 183 vs 804 ms), ToString from 200k digits cold (0.38–0.74×) and from 100k warm (0.30–0.57×).
- Below ~10k digits multiplication is 1.4–1.5× behind GMP (was 2.3–3.2×). Still trailing: division below ~200k-digit dividends (1.4–2.0×, was 2–2.9×) and small decimal I/O.

The margins come from the 2026 optimization campaign — NEON Shoup NTT butterflies, fused-MFA transforms, the wraparound-Newton division family, the Burnikel–Ziegler odd-size fix — and the 2026-11 multiplication pass: flag-chain limb kernels, all-NEON transform tails, a rewritten thread pool, a parallel finalize pass, truncated transform lengths and a re-swept dispatcher; and the 2026-11 division pass: a raw-buffer basecase, in-place Burnikel–Ziegler, a guarded-precision Newton reciprocal, a re-swept division dispatcher, quotient-only division and lock-free job pickup in the thread pool.

![BigMath vs GMP ratio by operand size](docs/images/bigmath_vs_gmp.png)

*Figure uses the canonical runs: multiplication 2026-11-20; division and decimal I/O 2026-11-21. Multiplication at ≥50M digits uses warm timings; ToString cold and measured warm timings are shown separately. See [BENCHMARK.md](BENCHMARK.md) for operand shapes and timing methodology.*

Full tables, methodology, and per-PR history: [BENCHMARK.md](BENCHMARK.md).

---

## Features

- **64-bit limb representation** (`BIGMATH_LIMB_64=1`, default). `DataT` stores true 64-bit values; every carry/borrow chain uses `__uint128_t` accumulators. Halves loop iteration counts in scalar paths vs the legacy 32-bit-in-64-bit layout. Opt out via `-DBIGMATH_LIMB_64=0`.
- **Multiplication:** Classic schoolbook → Karatsuba (64-bit-hybrid leaf) → 3-prime CRT NTT (2013265921 / 469762049 / 1811939329, 32-bit coefficient splitting) from 320 limbs per operand for balanced products with a warm pool, or 512 each with a cold pool. Toom-3 and the single-prime Goldilocks NTT remain as cross-check alternates; the NEON-accelerated CRT path beats both at every measured size.
- **Radix-4 + radix-8 fused NTT butterflies** (PRs #59, #60): 3× fewer memory passes vs radix-2; ~1.6× wall-clock at ≥2M limbs.
- **NEON Shoup butterflies on aarch64** (PRs #96–#99, `BIGMATH_NEON=1` default on Apple Silicon): the radix-8/4/2 layers and 1/n scaling run 4 lanes wide using Shoup's precomputed-reciprocal multiplication — `(w·x) mod P` from one doubling multiply-high (`vqdmulhq_s32`), a multiply-subtract, and one conditional subtract, with interleaved `(twiddle, w′)` tables keeping the strided gathers at one cache line per pair. Bit-exact with the scalar path; ~4× per butterfly pass, −26…40% end-to-end on NTT-bound ops. Scalar fallback everywhere else.
- **MFA / Bailey 6-step CRT NTT** (`BIGMATH_NTT_MFA_THRESHOLD=2^20`, default). Large CRT transforms switch to a cache-friendly matrix Fourier layout with fully fused stages: the transpose, the operand packing (straight from the limb arrays), and the pointwise product all ride the tiled stage passes, and the forward/inverse stages run as row-chunked `ParallelDo(6)` work units. The gate sat at `2^24` until the fused stages flipped the break-even (PR #109); cyclic wrap-around products use the same pipeline up to the CRT ceiling `2^26` (PR #110).
- **Multithreaded NTT** (`BIGMATH_USE_THREADS=1`, default). Small thread pool (size `min(hw_concurrency, BIGMATH_MAX_THREADS=8)`) parallelizes the CRT path: in the fused MFA window each transform stage runs as six (operand×prime or prime×row-range) work units; the non-MFA path batches 6 forwards + 3 inverses. 2.3-3.4× speedup on large mul / skewed div / parse. Opt out via `-DBIGMATH_USE_THREADS=0` to drop pthread linkage.
- **Division:** single-limb kernel with a normalized running remainder → Knuth Algorithm D on raw buffers (`FastDivision`: flag-chain `submul_1`, Möller–Granlund 3/2 quotient digit) below 96-limb divisors → in-place divide-and-conquer Burnikel–Ziegler (GMP `dcpi1` shape, no padding) → Newton–Raphson with a guarded-precision reciprocal and wrap-around cyclic products from 1536-limb divisors (768 at ratio ≥ 3) → quotient-sized division when the quotient is shorter than the divisor but ≥ 3072 limbs. Quotient-only calls (`operator/`) with a quotient shorter than the divisor are computed from the operands' top limbs. One rule per regime, re-swept 2026-11-21 on a divisor × ratio grid (`build/DispatchThresholds.h`). Identity `q·b + r == a` is cross-checked in `tests/div_correctness.cpp`.
- **Squaring:** Specialized Classic / Karatsuba squarers below 352 limbs; above, the CRT NTT self-multiply (skips the duplicate forward transforms — 2.5–5.8× over the former Goldilocks `NTTSquare`, which remains the `-DBIGMATH_NTT_CRT=0` fallback).
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

Apple M1 Max, vs GMP 6.3.0, `-O3 -march=native`, full default stack (`BIGMATH_LIMB_64=1` + `BIGMATH_NTT_CRT=1` + `BIGMATH_USE_THREADS=1`, 8-thread pool). Canonical `bench_vs_gmp` runs (multiplication 2026-11-20; division and decimal I/O 2026-11-21); rows marked *warm* discard the first call.

| operation | size (digits) | BigMath | GMP | ratio |
|---|---|---:|---:|---:|
| mul | 1 000 × 1 000 | 1.3 µs | 0.9 µs | 1.39× |
| mul | **10 000 × 10 000** | **27 µs** | **30 µs** | **0.90×** ← BigMath faster |
| mul | **100 000 × 100 000** | **0.25 ms** | **0.60 ms** | **0.42×** |
| mul | **1 000 000 × 1 000 000** | **4.8 ms** | **9.6 ms** | **0.50×** |
| mul | **10 000 000 × 10 000 000** | **43 ms** | **207 ms** | **0.21×** ← 4.8× faster |
| mul | **20 000 000 × 20 000 000** | **93 ms** | **276 ms** | **0.34×** |
| mul | **200 000 000 × 200 000 000** | **2 121 ms (warm)** | **2 957 ms** | **0.72×** |
| mul (skewed) | **1 000 000 / 100 000** | **1.4 ms** | **4.7 ms** | **0.29×** |
| mul (skewed) | **50 000 000 / 5 000 000** | **204 ms** | **677 ms** | **0.30×** |
| div (skewed) | **200 000 / 50 000** | **1.2 ms** | **1.7 ms** | **0.67×** ← BigMath faster |
| div (skewed) | **1 000 000 / 200 000** | **5.0 ms** | **10.2 ms** | **0.49×** |
| div (skewed) | **50 000 000 / 10 000 000** | **244 ms** | **1 307 ms** | **0.19×** ← 5.4× faster |
| div (skewed) | **200 000 000 / 40 000 000** | **1 450 ms** | **4 472 ms** | **0.32×** |
| parse | **1 000 000** | **7.1 ms** | **20.8 ms** | **0.34×** |
| parse | **20 000 000** | **183 ms** | **804 ms** | **0.23×** ← 4.4× faster |
| ToString | **100 000** | **1.4 ms (warm)** | **2.4 ms** | **0.57×** |
| ToString | **1 000 000** | **14.8 ms (warm)** | **49.6 ms** | **0.30×** |
| ToString | **10 000 000** | **362 ms (cold)** | **901 ms** | **0.40×** |

What remains behind GMP: division below ~200k-digit dividends (1.4–2.0×; worst for divisors of a few hundred limbs at high ratio, where GMP's divide-and-conquer sits on an assembly basecase and Toom multiplies), multiplication below ~10k digits (1.4–1.5×, GMP's assembly basecase plus Toom tiers), small decimal I/O, and 1–4-limb products where the result allocation dominates. See [BENCHMARK.md](BENCHMARK.md) for full tables, warm/cold methodology, and the condensed optimization history.

Opt-out flags (`-DBIGMATH_USE_THREADS=0` / `-DBIGMATH_NTT_CRT=0` / `-DBIGMATH_LIMB_64=0`) revert any subset of the defaults — useful for embedded targets, header-only-strict consumers, or A/B comparison.

## License

See [LICENSE](LICENSE).

## Author

S. M. Mahbub Murshed (murshed@gmail.com)
