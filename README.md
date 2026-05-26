# BigMath

Arbitrary-precision integer library in C++20. Header-light with a thin static-library shell. Goldilocks-prime NTT multiplication, Newton–Raphson division, divide-and-conquer base conversion. Includes a REPL calculator built on top.

Target: pure C++ that gets within 2–6× of GMP across the typical workload, without dropping to platform-specific assembly.

---

## Features

- **Multiplication:** Classic schoolbook → Karatsuba (64-bit hybrid leaf) → NTT with the Goldilocks prime `2⁶⁴ − 2³² + 1`. Single dispatcher in `algorithms/Multiplication.h` picks by operand size.
- **Division:** Classic short division → Knuth Algorithm D (`FastDivision`) → Burnikel–Ziegler (balanced large) → Newton–Raphson with reciprocal caching (skewed large). Identity `q·b + r == a` is cross-checked in `tests/div_correctness.cpp`.
- **Squaring:** Specialized Classic / Karatsuba / NTT squarers (1.4–1.6× over `Multiply(a,a)`).
- **String I/O:** Linear chunked parser/formatter for small inputs, divide-and-conquer with cached Newton reciprocals at scale. Asymptotic `O(M(L) · log L)` both directions.
- **BigDecimal:** Java-style fixed-point decimal (unscaled BigInteger + int scale) with exact +, −, \*; rounded division taking 8 rounding modes; parse/format covering plain and scientific notation.
- **Calculator REPL:** Variables, hex/bin/dec output, multi-line continuation, comments, `:help :quit :vars :reset :base :digits :time :load :save` directives. Built as a separate executable.

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

- [docs/BASE.md](docs/BASE.md) — number representation, base-2³² limbs, why little-endian
- [docs/MULTIPLICATION.md](docs/MULTIPLICATION.md) — Classic / Karatsuba / Toom-3 / NTT and their tradeoffs
- [docs/DIVISION.md](docs/DIVISION.md) — Classic / Fast (Knuth D) / Burnikel–Ziegler / Newton / Reciprocal-cached
- [docs/STRING_CONVERSION.md](docs/STRING_CONVERSION.md) — chunked decimal I/O, D&C parse/format, Newton-divider chain
- [docs/BIGDECIMAL.md](docs/BIGDECIMAL.md) — fixed-point decimal model, rounding modes, performance
- [BUILDING.md](BUILDING.md) — CMake build, install, threshold tuning

Each doc covers algorithms, dispatch, benchmark numbers vs GMP, optimizations that landed, and approaches that were tried and rejected with reasons.

## Performance

Apple M1 Max, vs GMP 6.3.0, `-O3 -march=native`:

| operation | size | BigMath | GMP | ratio |
|---|---|---:|---:|---:|
| mul | 5 000 × 5 000 digits | 0.043 ms | 0.012 ms | 3.5× |
| mul | 1 000 000 × 1 000 000 | 37.5 ms | 8.9 ms | 4.2× |
| div | 500 000 / 100 000 | 55 ms | 4.6 ms | 12× |
| parse | 100 000 digits | 7.4 ms | 1.1 ms | 7.1× |
| ToString | 100 000 digits | 41 ms | 2.3 ms | 18× |

The residual gap is the structural cost of pure C++ vs GMP's hand-tuned ARM64 assembly, plus the 32-bit-limb representation that keeps every inner-loop accumulator in a single register. NTT-bound operations hit the Goldilocks 16-bit coefficient floor regardless of representation. See the docs for full ratio tables and per-algorithm tuning notes.

## License

See [LICENSE](LICENSE).

## Author

S. M. Mahbub Murshed (murshed@gmail.com)
