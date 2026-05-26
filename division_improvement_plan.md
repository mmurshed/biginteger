# Division Improvement Plan

## Diagnosis

BigMath skewed division remains roughly 3-6x behind GMP in `BENCHMARK.md`, narrowing to about 3x at the largest tested sizes. The dominant path is Newton division. `docs/DIVISION.md` profiles the cost as mostly `NewtonDivision::DivideNormalizedWithReciprocal` plus `ApproxReciprocal`, both driven by NTT multiplication.

Small allocation-only changes are unlikely to move the result. The next useful work should focus on dispatch, recursive division shape, and reducing Newton's multiplication cost.

## 1. Audit Burnikel-Ziegler Dispatch Under Base2_64

Status: completed.

`src/algorithms/Division.cpp` currently gates BZ with `base == Base2_32`, while the default build uses `BIGMATH_LIMB_64=1` and `CurrentBase == Base2_64`. Verify whether this is intentional or a stale guard.

Benchmark direct calls to `FastDivision`, `BurnikelZieglerDivision`, `NewtonDivision`, and dispatcher for:

```text
1024/512, 2048/1024, 4096/2048, 8192/4096,
16384/8192, 32768/16384
```

If BZ is correct and faster for Base2_64, enable it in dispatch and update docs.

Finding: BZ is correct under Base2_64 and is much faster than FastDivision on near-balanced shapes. Dispatch now allows BZ for `Base2_64`.

## 2. Add Shape-Focused Division Benchmarks

Status: completed.

Current benchmark coverage emphasizes equal-size short-circuit cases and skewed Newton cases. Add a benchmark focused on meaningful division shapes:

```text
2n/n, 3n/2n, 5n/2n, 3n/n, 5n/n, 10n/n
```

Record direct algorithm timings and dispatcher choice. Use this before retuning thresholds.

Added `tests/performance/division_shape_bench.cpp`.

## 3. Retune Dispatch Thresholds

Status: completed.

After benchmark coverage exists, regenerate thresholds for:

- Newton divisor lower bound
- Newton skew ratio
- BZ divisor threshold
- BZ recursion threshold
- Base2_64-specific BZ eligibility

Avoid tightening Newton back to `4*b`; prior docs show that caused severe boundary regressions.

Finding: Newton was entering too early for 1024-2048 limb divisors. The retune uses:

- medium-skew Newton: `b >= 4096` and `a >= 3b`
- high-skew Newton: `b >= 2048` and `a >= 8b`
- BZ: enabled for both `Base2_32` and `Base2_64`

## 4. Prototype GMP-Style Pre-Inverted D&C Division

Status: prototyped and rejected for now.

The largest structural gap is Newton chunked division versus GMP's recursive pre-inverted division. Prototype a new recursive path for large `2n/n` and near-balanced inputs, similar in spirit to GMP's `mpn_dcpi1_div_q`.

Target:

- precompute inverse once
- compute quotient blocks recursively
- avoid Newton's repeated full reciprocal chunk multiplications where possible
- benchmark against BZ and Newton

This is the highest-ceiling improvement, but also the highest-risk change.

Experiment: implemented a precomputed BZ divisor tree that cached divisor splits and high-half recursive dividers across BZ blocks. It measured mixed/flat against the existing BZ implementation, so it was reverted. A true GMP-style `mpn_dcpi1_div_q` equivalent still requires a deeper algorithmic implementation, not just cached BZ structure.

## 5. Investigate True Truncated Products

Status: investigated, not implemented.

Newton often only needs the high half of products such as `chunk * reciprocal`. The previous Mulders-style high product was rejected because it decomposed one NTT into two NTTs and measured flat.

Only revisit this if the implementation reduces actual NTT work, for example with a true high-window convolution or truncated NTT API.

Finding: the current multiplication API always computes full products, and the prior exact `MulHigh` experiment decomposed one NTT into two NTTs and measured flat. No production change was made.

## 6. Improve Medium-Band FastDivision

FastDivision still matters below Newton/BZ thresholds and for small divisors. Benchmark and consider:

- quotient-only path when remainder is not needed
- tighter Base2_64 subtract-multiply pointer loops
- small divisor bands: `m/2`, `m/4`, `m/8`, `m/16`, `m/32`, `m/64`
- avoiding temporary resizing in hot loops

This will not close the large Newton gap, but can reduce medium-size regressions.

## Avoid For Now

- Reattempting the previous Mulders `MulHigh`; it measured flat.
- Division-layer scratch buffers without output-buffer multiplication/subtraction; a ToString/Newton scratch experiment regressed.
- Lowering BZ recursion threshold blindly below 512 limbs.
- Removing `NewtonDivision::Divider`; cached reciprocal division is critical for repeated-divisor workloads and ToString.

## Recommended Order

1. Audit and benchmark BZ with Base2_64.
2. Add shape-focused division benchmarks.
3. Retune dispatch thresholds from measurements.
4. Prototype pre-inverted recursive division for `2n/n`.
5. Investigate true truncated NTT/high-product support only after the above.
