# BigMath

Arbitrary-precision integer library in C++20. Header-light with a thin static-library shell. 64-bit limbs, 3-prime CRT NTT multiplication (multithreaded by default), Newton–Raphson division with cached reciprocals, divide-and-conquer base conversion. Includes a REPL calculator built on top.

Target: pure C++ that matches GMP on skewed multiplication and gets within 1.4-5× on the rest, without dropping to platform-specific assembly.

---

## Features

- **64-bit limb representation** (`BIGMATH_LIMB_64=1`, default). `DataT` stores true 64-bit values; every carry/borrow chain uses `__uint128_t` accumulators. Halves loop iteration counts in scalar paths vs the legacy 32-bit-in-64-bit layout. Opt out via `-DBIGMATH_LIMB_64=0`.
- **Multiplication:** Classic schoolbook → Karatsuba (64-bit-hybrid leaf) → 3-prime CRT NTT (998244353 / 985661441 / 754974721, 30-bit primes, 32-bit coefficient splitting) gated above 5000 limbs sum; Goldilocks NTT for smaller cases. Single dispatcher in `algorithms/Multiplication.h` picks by operand size.
- **Radix-4 + radix-8 fused NTT butterflies** (PRs #59, #60). Adjacent radix-2 layers collapse into single load/store butterflies — radix-4 fuses 2 layers (4 elements), radix-8 fuses 3 layers (8 elements). Same modular op count; 3× fewer memory passes vs radix-2. ~1.6× wall-clock at ≥2M limbs, widens the BigMath-wins band against GMP from 5M-10M to **2M-20M**.
- **Multithreaded NTT** (`BIGMATH_USE_THREADS=1`, default). Small thread pool (size `min(hw_concurrency, BIGMATH_MAX_THREADS=8)`) parallelizes the CRT path: 6 forwards + 3 inverses as batched work units, one `ParallelDo` dispatch per phase. 2.3-3.4× speedup on large mul / skewed div / parse. Opt out via `-DBIGMATH_USE_THREADS=0` to drop pthread linkage.
- **Division:** Classic short division → Knuth Algorithm D (`FastDivision` with Möller-Granlund 3-by-2 qhat for Base2_32) → Burnikel–Ziegler (balanced large) → Newton–Raphson with reciprocal caching (skewed large). Identity `q·b + r == a` is cross-checked in `tests/div_correctness.cpp`.
- **Squaring:** Specialized Classic / Karatsuba / NTT squarers (1.4–1.6× over `Multiply(a,a)`).
- **String I/O:** Linear chunked parser/formatter for small inputs, divide-and-conquer with cached Newton reciprocals at scale. Asymptotic `O(M(L) · log L)` both directions.
- **BigDecimal:** Java-style fixed-point decimal (unscaled BigInteger + int scale) with exact +, −, \*; rounded division taking 8 rounding modes; parse/format covering plain and scientific notation.
- **Calculator REPL:** Variables, hex/bin/dec output, multi-line continuation, comments, `:help :quit :vars :reset :base :digits :time :load :save` directives. Built as a separate executable.
- **Thread-safe by construction** for concurrent use of distinct objects from distinct threads. See [docs/THREAD_SAFETY.md](docs/THREAD_SAFETY.md).

## Layout

```
include/biginteger/    public BigInteger headers (algorithms, ops, common)
include/bigdecimal/    public BigDecimal headers
src/                   non-inline BigInteger implementations
bigdecimal/            BigDecimal implementation
include/biginteger/    public headers (algorithms, ops, common)
src/                   non-inline implementations
tests/                 correctness, performance, unit tests
calculator/            REPL frontend
docs/                  technical references (see below)
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

Apple M1 Max, vs GMP 6.3.0, `-O3 -march=native`, full default stack (`BIGMATH_LIMB_64=1` + `BIGMATH_NTT_CRT=1` + `BIGMATH_USE_THREADS=1`, 8-thread pool). Min of 5 runs:

| operation | size | BigMath | GMP | ratio |
|---|---|---:|---:|---:|
| mul | 100 000 × 100 000 | 1.20 ms | 0.76 ms | 1.57× |
| mul | 1 000 000 × 1 000 000 | 9.84 ms | 8.89 ms | 1.11× |
| mul | **2 000 000 × 2 000 000** | **19.8 ms** | **20.8 ms** | **0.95×** ← BigMath faster |
| mul | **5 000 000 × 5 000 000** | **45.8 ms** | **64.0 ms** | **0.72×** ← BigMath faster |
| mul | **10 000 000 × 10 000 000** | **103 ms** | **213 ms** | **0.48×** ← BigMath 2.08× faster |
| mul | **20 000 000 × 20 000 000** | **276 ms** | **280 ms** | **0.99×** ← parity |
| mul | 50 000 000 × 50 000 000 | 1 393 ms | 681 ms | 2.05× ← GMP SSA recovers |
| mul | 100 000 000 × 100 000 000 | 3 221 ms | 1 449 ms | 2.22× |
| mul (skewed) | **1 000 000 / 100 000** | **4.27 ms** | **4.48 ms** | **0.95×** ← BigMath faster |
| mul (skewed) | **2 000 000 / 200 000** | **8.51 ms** | **9.40 ms** | **0.91×** ← BigMath faster |
| mul (skewed) | 10 000 000 / 1 000 000 | 84.2 ms | 70.7 ms | 1.19× |
| mul (skewed) | **50 000 000 / 5 000 000** | **625 ms** | **666 ms** | **0.94×** ← BigMath faster |
| div (skewed) | 500 000 / 100 000 | 18.2 ms | 4.76 ms | 3.84× |
| div (skewed) | 10 000 000 / 2 000 000 | 412 ms | 157 ms | 2.63× |
| div (skewed) | 50 000 000 / 10 000 000 | 2 356 ms | 1 312 ms | **1.79×** |
| parse | 1 000 000 digits | 46.5 ms | 20.5 ms | 2.27× |
| parse | 20 000 000 digits | 1 243 ms | 800 ms | **1.55×** |
| ToString | 100 000 digits | 19.5 ms | 2.36 ms | 8.24× |
| ToString | 1 000 000 digits | 228 ms | 48.7 ms | 4.67× |
| ToString | 20 000 000 digits | 5 238 ms | 2 115 ms | **2.48×** |

**BigMath beats GMP on balanced multiplication across the 2M–20M digit band** — the NTT asymptotic edge plus radix-4+8 fused butterflies (PRs #59, #60) widen the win window and push the **10M peak to 2.08× faster than GMP** (103 ms vs 213 ms). At ≥50M GMP's Schönhage-Strassen reverses the lead but the loss factor halved (3.58× → 2.05× at 50M). Skewed mults win four bands (500k×50k through 50M×5M, 0.91-0.97×). Skewed division at 50M×10M dropped to **1.79×** as Newton inherits the NTT speedup via its per-chunk muls. ToString narrows from 8.2× at 100k to 2.48× at 20M; parse to 1.55× at 20M. See [BENCHMARK.md](BENCHMARK.md) for the full table or the per-doc ratio tables for the breakdown.

Opt-out flags (`-DBIGMATH_USE_THREADS=0` / `-DBIGMATH_NTT_CRT=0` / `-DBIGMATH_LIMB_64=0`) revert any subset of the defaults — useful for embedded targets, header-only-strict consumers, or A/B comparison.

## License

See [LICENSE](LICENSE).

## Author

S. M. Mahbub Murshed (murshed@gmail.com)
