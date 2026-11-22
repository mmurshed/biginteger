# Building BigMath

For prebuilt libraries and the GitHub release workflow, see the [release guide](docs/RELEASING.md).

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
- C++20 compiler with `__uint128_t` and `__builtin_*_overflow`. The release matrix tests GCC on Ubuntu 22.04, Clang on Ubuntu 24.04, Apple Clang, and MSYS2 UCRT64 Clang.
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

Compares this library to GMP across multiplication, division, parse, and ToString (1k digits up through the multi-million-digit rows; see BENCHMARK.md for the full canonical run). The bench takes ~30 s.

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

See the [platform tuning guide](docs/PLATFORM_TUNING.md) for running the dispatch tuner and validating profiles on the current machine or GitHub Actions.

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

Expected: see the current ratio table in README.md (post-#116-#119 BigMath beats GMP at ≥500k digits on mul, div, parse, and ToString; e.g. `ToString 100 000 digits` ≈ 1.25× vs GMP). See [docs/MULTIPLICATION.md](docs/MULTIPLICATION.md) and [docs/DIVISION.md](docs/DIVISION.md) for the full reference table.
