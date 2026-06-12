# Repository Guidelines

## Project Structure & Module Organization

This is a CMake-built C++ BigInteger library in namespace `BigMath` (`include/` + `src/` split; see BUILDING.md).

- `include/biginteger/BigInteger.h` defines the main `BigInteger` wrapper.
- `include/biginteger/common/` contains shared types, constants, parsing, comparison, and utilities.
- `include/biginteger/ops/` (+ `src/ops/`) contains public `BigInteger` operators and sign-aware wrappers.
- `include/biginteger/algorithms/` (+ `src/algorithms/`) contains limb-vector algorithms and dispatchers.
- `include/biginteger/algorithms/multiplication/` and `.../division/` contain concrete arithmetic implementations.
- `include/bigdecimal/` + `bigdecimal/` contain the fixed-point decimal layer.
- `tests/` contains unit tests plus standalone correctness and performance programs.
- `calculator/` contains the interactive calculator frontend.
- `bigtst.sln` and `bigtst.vcxproj` are legacy Visual Studio project files.

## Build, Test, and Development Commands

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
cd build && ctest --output-on-failure
./build/mult_correctness && ./build/div_correctness
```

Use GCC or Clang. Several optimized paths rely on `unsigned __int128`.

## Coding Style & Naming Conventions

Follow the existing header-only style: classes use `PascalCase`, methods use `PascalCase`, and type aliases use names like `DataT`, `SizeT`, and `BaseT`. Keep implementations inside `namespace BigMath`. Prefer small helper methods near the algorithm class they support. Avoid unrelated formatting churn, especially in older CRLF-formatted files.

## Testing Guidelines

Correctness tests are executable `.cpp` files, not a unit-test framework. When changing multiplication or division, run the relevant correctness test and at least one simple performance harness. Division results must satisfy `q * b + r == a` and `r < b`. Multiplication algorithms should compare limb-for-limb against the reference implementation.

## Commit & Pull Request Guidelines

Recent commits use short imperative or descriptive subjects, for example `Branchless reduction` and `Blockwise Newton`. Keep commits focused on one algorithm or bug fix. Pull requests should describe the changed algorithm, list benchmark/correctness commands run, and call out any dispatch-threshold changes.

## Architecture Notes

Numbers are little-endian vectors of `DataT`. The default internal radix is `Base2_64` (64-bit limbs; `-DBIGMATH_LIMB_64=0` for the legacy `Base2_32` layout), even though `DataT` is 64-bit, leaving upper bits for carry headroom. Public operators live in `ops/`; algorithm selection belongs in `algorithms/Multiplication.h` and `algorithms/Division.h`.
