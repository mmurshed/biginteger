# BigMath

Arbitrary-precision integer library in C++20. Header-light with a thin static-library shell. 64-bit limbs, 3-prime CRT NTT multiplication (multithreaded by default), Newton–Raphson division with cached reciprocals, divide-and-conquer base conversion. Includes a REPL calculator built on top.

Target: portable C++ (with optional NEON intrinsics on aarch64, scalar fallbacks everywhere else) that reaches GMP parity or better across the large-operand bands — no hand-written assembly.

## Performance

Both runs used the same operands on an Apple M1 Max and share one GMP 6.3.0 pass, which used one thread. Values below 1× mean BigMath was faster. The arithmetic panels show warm calls, and the decimal conversion panels show first and warm calls. Most large cases had one warm measurement per library.

### Single-threaded BigMath and GMP

BigMath’s worker pool was disabled at compile time. The timings below are from warm calls.

- **Multiplication.** GMP was faster at 13 of the 14 equal-length sizes. BigMath won at 10 million digits per operand (137 ms versus 202 ms). At 200 million digits per operand, GMP took 2.87 s and BigMath took 10.75 s.
- **Division.** BigMath was faster in both output modes at 50 million ÷ 10 million and 100 million ÷ 20 million digits. At 200 million ÷ 40 million, BigMath took 1.46× GMP’s time for the quotient and 1.37× for quotient and remainder.
- **Decimal conversion.** GMP was faster at every sampled size for both parsing and formatting. At 20 million digits, BigMath took 938 ms to parse and 2.87 s to format, versus 800 ms and 2.10 s for GMP.

![Single-thread BigMath time divided by single-thread GMP time for multiplication, division, decimal parsing and formatting](docs/images/bigmath_vs_gmp_single_thread.png)

### Multithreaded BigMath (up to eight threads) and single-threaded GMP

![Default-build BigMath time divided by single-thread GMP time for multiplication, division, decimal parsing and formatting](docs/images/bigmath_vs_gmp.png)

### Small operands (1–16 limbs)

Values of up to 4 limbs (256 bits) are stored inside the `BigInteger` object, and operands of up to 16 limbs are combined in a stack buffer. Before and after that change, single-threaded, against GMP writing into a reused output: a 1-limb product went from 5.7× to 2.1× GMP's time, a 1-limb remainder from 3.9× to 1.1×, and copying a small value from 5.3× to 0.12×.

![BigMath time divided by GMP time on 1 to 16 limb operands, before and after inline storage](docs/images/bigmath_vs_gmp_small.png)

[Single-thread data and methodology](BENCHMARK.md#single-thread-comparison) · [Multithreaded data and methodology](BENCHMARK.md#default-build-with-up-to-eight-bigmath-threads) · [Small-operand data](BENCHMARK.md#small-operands)

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

The internal multiplication table reports the mean of five calls per algorithm. [Internal algorithm crossovers](BENCHMARK.md#internal-algorithm-crossovers)

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

Warm wall-clock times from the November 21, 2026 Apple M1 Max run. The first ratio compares the two single-thread timings shown. The last column gives the ratio for the default build (up to eight BigMath threads) against the same GMP timings. Ratios use the unrounded times. [Single-thread CSV](docs/benchmarks/2026-11-21-sbo-m1-max-single-thread.csv) · [Default-build CSV](docs/benchmarks/2026-11-21-sbo-m1-max.csv)

| Operation | Decimal digits | BigMath 1 thread (ms) | GMP 1 thread (ms) | 1-thread ratio | Default-build ratio |
|---|---:|---:|---:|---:|---:|
| Multiply | 10,000 × 10,000 | 0.04971 | 0.02792 | 1.78× | 0.88× |
| Multiply | 1,000,000 × 1,000,000 | 10.2 | 8.977 | 1.14× | 0.31× |
| Multiply | 5,000,000 × 5,000,000 | 62.33 | 62.29 | 1.00× | 0.20× |
| Multiply | 10,000,000 × 10,000,000 | 137.1 | 201.8 | 0.68× | 0.14× |
| Multiply | 20,000,000 × 20,000,000 | 312 | 267 | 1.17× | 0.27× |
| Multiply | 200,000,000 × 200,000,000 | 10,746 | 2,865 | 3.75× | 0.76× |
| Multiply | 5,000,000 × 500,000 | 23.53 | 30.79 | 0.76× | 0.27× |
| Multiply | 200,000,000 × 20,000,000 | 4,923 | 1,911 | 2.58× | 0.53× |
| Divide, quotient | 200,000 ÷ 50,000 | 2.552 | 1.512 | 1.69× | 0.70× |
| Divide, quotient | 10,000,000 ÷ 2,000,000 | 160.3 | 145 | 1.11× | 0.36× |
| Divide, quotient | 50,000,000 ÷ 10,000,000 | 993.7 | 1,220 | 0.81× | 0.20× |
| Divide, quotient | 200,000,000 ÷ 40,000,000 | 6,095 | 4,189 | 1.46× | 0.32× |
| Divide, quotient and remainder | 200,000 ÷ 50,000 | 2.549 | 1.684 | 1.51× | 0.68× |
| Divide, quotient and remainder | 10,000,000 ÷ 2,000,000 | 157.2 | 154.6 | 1.02× | 0.34× |
| Divide, quotient and remainder | 50,000,000 ÷ 10,000,000 | 991.6 | 1,297 | 0.76× | 0.18× |
| Divide, quotient and remainder | 200,000,000 ÷ 40,000,000 | 6,100 | 4,441 | 1.37× | 0.30× |
| Parse | 100,000 | 1.625 | 1.058 | 1.54× | 0.97× |
| Parse | 20,000,000 | 937.5 | 799.8 | 1.17× | 0.22× |
| ToString | 100,000 | 3.858 | 2.271 | 1.70× | 0.59× |
| ToString | 20,000,000 | 2,874 | 2,096 | 1.37× | 0.29× |

Opt-out flags (`-DBIGMATH_USE_THREADS=0` / `-DBIGMATH_NTT_CRT=0` / `-DBIGMATH_LIMB_64=0`) revert any subset of the defaults — useful for embedded targets, header-only-strict consumers, or A/B comparison.

## License

See [LICENSE](LICENSE).

## Author

S. M. Mahbub Murshed (murshed@gmail.com)
