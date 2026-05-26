# Multiplication Optimization Plan

## Summary

Optimize the portable C++ multiplication path only. Do not add ARM64 assembly, do not re-enable Toom-Cook dispatch, and do not start a full 64-bit limb refactor. Focus on NTT overhead, square/multiply shared planning, and benchmark-driven thresholds.

## Key Changes

- Extract duplicated NTT internals from `NTTMultiplication` and `NTTSquare` into one internal NTT core helper.
- Replace separate forward-root, inverse-root, and per-call `n_inv` work with one thread-local transform plan keyed by transform size: forward roots, inverse roots, and inverse length.
- Rework Base2_32 NTT finalization:
  - carry only across the exact convolution coefficient count, not all padded transform slots.
  - pack 16-bit coefficients directly into 32-bit limbs during carry propagation.
  - remove the intermediate normalized `c16` then `c32` pass.
- Add reusable thread-local NTT work buffers for `fa` and `fb` to reduce repeated vector allocation churn in repeated same-size multiplication and Newton-driven workloads.
- Apply the same shared-plan and direct-finalization path to `NTTSquare`, preserving its single-forward-transform advantage.
- Add `dispatch_tuner` as a CMake performance target if it is still not built by CMake, but do not register it as a `ctest`.

## Dispatch Tuning

- Re-run thresholds after NTT changes using Release `-O3 -march=native`.
- Set `BIGMATH_NTT_MULTIPLICATION_THRESHOLD` to the first equal-size total-limb point where NTT beats Karatsuba by at least 5% for two consecutive measured sizes.
- Set `BIGMATH_CLASSIC_MIN_LIMB_THRESHOLD` to the highest small operand size where Classic beats both Karatsuba and NTT by at least 5% in skewed `large x small` cases.
- Keep `BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD = 0` unless Classic beats Karatsuba by at least 5% in the equal-size sweep.
- Add or extend square tuning so `BIGMATH_NTT_SQUARE_THRESHOLD` is chosen from measured `KaratsubaSquare` vs `NTTSquare`, using the same two-consecutive-size and 5% rule.

## Test Plan

Correctness:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8
cd build && ctest --output-on-failure
./mult_correctness
./unit_tests
```

Performance:

```sh
./build/dispatch_tuner --full
./build/multperf_simple
./build/bench_vs_gmp
```

Run `bench_vs_gmp` only when GMP is available.

Benchmark focus cases:

- equal-size multiply around 512, 750, 1024, 1500, 2048, and 4096 limbs.
- skewed multiply such as `4096 x 4`, `4096 x 16`, `4096 x 48`, and `4096 x 64`.
- square around 256, 512, 750, 1024, and 4096 limbs.
- at least one large decimal `ToString` or parse scenario because `Pow10` and Newton division consume multiplication heavily.

## Assumptions

- Priority is portable C++ performance.
- Public BigInteger APIs remain unchanged.
- Toom-Cook remains a correctness cross-check implementation, not a dispatcher candidate.
- Large NTT multiplication is already near its pure-C++ structural limit, so expected wins are modest: mostly allocation, cold-plan, finalization, and threshold improvements.
