# Building BigMath

Production-ready CMake build, replacing the prior header-only `c++ -std=c++20 …` per-test recipes.

## Layout

```
biginteger/                     repository root
├── CMakeLists.txt              top-level build file
├── include/biginteger/         public headers (consumed via -Iinclude)
│   ├── BigInteger.h
│   ├── common/
│   ├── algorithms/
│   │   ├── multiplication/
│   │   └── division/
│   └── ops/
├── src/                        non-inline implementations
│   ├── common/
│   ├── algorithms/
│   └── ops/
├── tests/                      correctness + perf harnesses
├── docs/                       MULTIPLICATION.md, DIVISION.md, …
└── build/                      out-of-source build output (gitignored)
```

## Requirements

- CMake ≥ 3.16
- C++20 compiler with `__uint128_t` and `__builtin_*_overflow`. GCC ≥ 10 or Clang ≥ 11. Apple Clang on macOS works.
- For the GMP comparison bench: `libgmp` (Homebrew: `brew install gmp`).

## Quick build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
```

The build emits:

- `build/libbigmath.a` — static library, the primary deliverable.
- `build/{mult_correctness, div_correctness, multperf_simple, divperf_simple, bench_vs_gmp}` — test and benchmark binaries linked against `bigmath`.

`-DBIGMATH_BUILD_SHARED=ON` also produces `libbigmath.{so|dylib}` for shared-library distribution.

## Testing

```sh
cd build && ctest --output-on-failure
```

`mult_correctness` runs in <1 s. `div_correctness` covers a broader matrix and takes ~2 min on Apple M1 Max.

## Benchmark

```sh
./build/bench_vs_gmp
```

Compares this library to GMP across multiplication, division, parse, and ToString at sizes 1 000 – 1 000 000 decimal digits. The bench takes ~30 s.

Results discussed in [docs/MULTIPLICATION.md](docs/MULTIPLICATION.md), [docs/DIVISION.md](docs/DIVISION.md), [docs/STRING_CONVERSION.md](docs/STRING_CONVERSION.md).

## Tuning thresholds at build time

All dispatch thresholds and per-algorithm cutoffs are macro-overridable:

```sh
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_FLAGS_RELEASE="-O3 -march=native \
        -DBIGMATH_NTT_MULTIPLICATION_THRESHOLD=1500 \
        -DBIGMATH_KARATSUBA_THRESHOLD=48 \
        -DBIGMATH_NEWTON_MEDIUM_B=1024"
```

See `tests/performance/dispatch_tuner.cpp` for an empirical sweep that recommends values for the current machine.

## Consuming the library from another CMake project

After `cmake --install build --prefix /some/where`, downstream projects can do:

```cmake
find_package(bigmath REQUIRED)        # or via add_subdirectory(...)
target_link_libraries(my_app PRIVATE bigmath::bigmath)
```

The header path resolves as `<biginteger/BigInteger.h>`.

## Sanity check the refactor

```sh
# Build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8

# Correctness
cd build && ctest --output-on-failure && cd ..

# Bench: any ratio in the ±5% band of the documented numbers is "within noise"
./build/bench_vs_gmp
```

Expected: `mul 1000×1000` ≈ 3.5× vs GMP; `ToString 100 000 digits` ≈ 17× vs GMP. See [docs/MULTIPLICATION.md](docs/MULTIPLICATION.md) and [docs/DIVISION.md](docs/DIVISION.md) for the full reference table.
