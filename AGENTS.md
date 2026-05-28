# Repository Guidelines

## Project Structure & Module Organization

This is a header-only C++ BigInteger library in namespace `BigMath`.

- `biginteger/BigInteger.h` defines the main `BigInteger` wrapper.
- `biginteger/common/` contains shared types, constants, parsing, comparison, and utilities.
- `biginteger/ops/` contains public `BigInteger` operators and sign-aware wrappers.
- `biginteger/algorithms/` contains limb-vector algorithms and dispatchers.
- `biginteger/algorithms/multiplication/` and `biginteger/algorithms/division/` contain concrete arithmetic implementations.
- `tests/` contains standalone correctness and performance programs.
- `calculator/` contains the interactive calculator frontend.
- `bigtst.sln` and `bigtst.vcxproj` are Visual Studio project files.

## Build, Test, and Development Commands

There is no Makefile or CMake project. Compile tests directly:

```sh
c++ -std=c++20 -O2 tests/mult_correctness.cpp -o /tmp/mult_correctness && /tmp/mult_correctness
c++ -std=c++20 -O2 tests/div_correctness.cpp -o /tmp/div_correctness && /tmp/div_correctness
c++ -std=c++20 -O2 tests/multperf_simple.cpp -o /tmp/multperf_simple && /tmp/multperf_simple
c++ -std=c++20 -O2 tests/divperf_simple.cpp -o /tmp/divperf_simple && /tmp/divperf_simple
c++ -std=c++20 -O2 calculator/calculator.cpp -o /tmp/calculator && /tmp/calculator
```

Use GCC or Clang. Several optimized paths rely on `unsigned __int128`.

## Coding Style & Naming Conventions

Follow the existing header-only style: classes use `PascalCase`, methods use `PascalCase`, and type aliases use names like `DataT`, `SizeT`, and `BaseT`. Keep implementations inside `namespace BigMath`. Prefer small helper methods near the algorithm class they support. Avoid unrelated formatting churn, especially in older CRLF-formatted files.

## Testing Guidelines

Correctness tests are executable `.cpp` files, not a unit-test framework. When changing multiplication or division, run the relevant correctness test and at least one simple performance harness. Division results must satisfy `q * b + r == a` and `r < b`. Multiplication algorithms should compare limb-for-limb against the reference implementation.

## Commit & Pull Request Guidelines

Recent commits use short imperative or descriptive subjects, for example `Branchless reduction` and `Blockwise Newton`. Keep commits focused on one algorithm or bug fix. Pull requests should describe the changed algorithm, list benchmark/correctness commands run, and call out any dispatch-threshold changes.

## Architecture Notes

Numbers are little-endian vectors of `DataT`. The current internal radix is `Base2_32`, even though `DataT` is 64-bit, leaving upper bits for carry headroom. Public operators live in `ops/`; algorithm selection belongs in `algorithms/Multiplication.h` and `algorithms/Division.h`.
