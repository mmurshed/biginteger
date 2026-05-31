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

## 3b. Near-Balanced Newton Band + BZ Non-Power-of-2 Blowup (PR #79)

Status: completed.

Profiling near-balanced (ratio ≈ 2) division at large `b` found two compounding bugs:

1. **Newton single-block pathology.** The single-block path (`na ≤ 2n+1`) ran a `2n+1`-limb chunk through the truncated `(chunk·R) >> 2n` estimate; the error scales with `chunk / B^(2n)` and reaches ~`B` once the chunk exceeds `2n` limbs (the `+1` comes from the Knuth normalize shift on `a ≈ 2n`), overflowing the fixup cap and bailing to quadratic FastDivision (23× spike at `nb = 50000`). Fixed by routing `na > 2n` through the blockwise path.

2. **BZ non-power-of-2 blowup.** BZ's recursive 2n/n halving lands its intermediate NTT multiplies just over power-of-2 transform-length boundaries for non-power-of-2 divisor sizes; the FFT length doubles and the constant factor compounds across recursion depth into a **5–60× slowdown vs Newton**, worst at `n = 2^k+1` (≈90 s for a 262145-limb divisor vs Newton ~0.9 s).

Fix: new Newton balanced band — `b >= NEWTON_BALANCED_B (98304)` and `a >= 2b` → Newton. Measured 8.8× at 500k/250k, 5.6× at 2M/1M, 97.7× at the 2¹⁸+1 worst case; exact-pow2 sizes (BZ best case) regress ~4%. Harness: `tests/performance/division_balanced_bench.cpp`.

**Residual / next step:** ratio ∈ (1, 2) at large `b` still routes to BZ and hits the same blowup (~2.7× slower than Newton at ratio 1.5). Extending the balanced band below ratio 2 needs a quotient-bulk lower bound so genuinely tiny-quotient `a ≈ b` cases (where a full reciprocal is wasteful) stay on FastDivision.

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
