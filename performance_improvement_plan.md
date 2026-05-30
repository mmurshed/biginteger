# Performance Improvement Plan

## Diagnosis

BigMath has already closed a meaningful part of the GMP gap in the large balanced multiplication bands, and the recent dispatch work fixed one real skewed-multiplication miss. The remaining room for improvement is now less about choosing a different asymptotic algorithm and more about removing overhead in the hot paths that still run on every large operation.

The strongest current candidates are:

- allocation and temporary-vector churn inside multiplication and Newton division
- repeated setup work that can be cached or prepared once
- machine-specific threshold drift between builds and architectures
- the remaining skewed multiplication band where a blockwise strategy might outperform one full symmetric multiply

The plan below is ordered by likely return on effort.

## 1. Reuse Scratch Buffers In Hot Paths

Status: proposed.

The first improvement to pursue is reusable scratch space in the code that currently allocates and frees temporary vectors repeatedly:

- `NTTMultiplication`
- `NTTMultiplicationCrt`
- `NewtonDivision`
- any helper that repeatedly builds short-lived coefficient buffers or carry buffers

The goal is to reduce `std::vector` construction, resize, and trimming overhead in the paths that dominate large multiplication, skewed multiplication, division, parse, and ToString.

What to measure:

- balanced large multiplication
- skewed multiplication in the pre-NTT range
- Newton division on large divisors
- parse and ToString on very large inputs

Success criteria:

- lower wall-clock time without changing dispatch behavior
- no correctness regressions in `mult_correctness`, `div_correctness`, and the relevant unit tests
- no increase in code complexity that leaks through the public API

## 2. Extend Prepared and Cached APIs

Status: completed.

The library already benefits from internal caching in several places. The next step is to broaden that idea where the same operand or setup is reused many times.

Good candidates:

- prepared NTT multiplication for repeated products against one fixed operand
- reusable reciprocal state for repeated division by the same divisor
- optional cache-aware APIs for internal callers such as `ToString` and `BigDecimal`

This will not help one-off `a * b` calls as much as the scratch-buffer work, but it can remove repeated setup cost in real workloads that reuse the same large operand.

Implemented as public wrappers in `biginteger/ops/Multiplication.h` and `biginteger/ops/Division.h`:

- `PreparedMultiplication` for repeated `BigInteger` multiplies against one prepared operand
- `CachedDivision` for repeated `BigInteger` division against one cached divisor

Both wrappers preserve the existing sign handling and are covered by unit tests.

What to measure:

- repeated large multiplies with one fixed operand
- repeated division by the same divisor
- warm-cache versus cold-cache behavior

Success criteria:

- measurable improvement for repeated workloads
- no regression for one-off calls
- an API that stays explicit about when reuse is valid

## 3. Auto-Tune Dispatch Thresholds Per Build

Status: completed.

The current thresholds are reasonable on the machine they were measured on, but the best cutoffs can move with compiler, CPU, and limb layout. The library already has benchmark harnesses that make this feasible.

The next step is to make threshold selection more systematic:

- collect shape-aware benchmark data during a local tuning run
- emit suggested compile-time overrides for multiplication and division dispatch
- keep the defaults conservative when no tuning pass is run

This matters most for:

- `CLASSIC_MULTIPLICATION_THRESHOLD`
- `TOOM3_MULTIPLICATION_THRESHOLD`
- `NTT_MULTIPLICATION_THRESHOLD`
- `NTT_SQUARE_THRESHOLD`
- `BZ_DIVISOR_THRESHOLD`
- Newton skew thresholds

Success criteria:

- threshold suggestions are derived from measurements, not hand-tuned guesswork
- the defaults remain sensible for general users
- local builds on different machines can be tuned without code changes

Implemented as a build-local threshold profile in `include/biginteger/build/DispatchThresholds.h`, plus an emission mode in `tests/performance/dispatch_tuner.cpp`:

- the dispatcher now includes a single shared threshold header before the algorithm-specific fallbacks
- `dispatch_tuner --emit-header [path]` writes a build-specific override file after measuring the current machine
- the checked-in header keeps conservative defaults, while local tuning can regenerate it for the active build

## 4. Revisit Extreme Skew Multiplication

Status: research candidate.

The recent skew gate fixed a real dispatch bug: Toom-3 is good in the balanced band, but not for 2:1+ skewed inputs there. The next open question is whether a blockwise skew strategy can be made to win consistently for very large ratios without regressing the middle band.

This should only be pursued if the benchmark shape shows a clear gap:

- very large `n/m` ratios
- large enough smaller operand to make one-shot Classic unattractive
- a regime where Karatsuba and one full NTT are both leaving performance on the table

The previous prototype already showed that a naive blockwise strategy can win in narrow cases and then fall off quickly. Any revival should be benchmark-driven and gated narrowly, not promoted by instinct.

Success criteria:

- measurable win in a clearly defined skewed band
- no regression in balanced or moderately skewed cases
- dispatcher rules that remain easy to explain

## 5. Optimize Finalization Only If It Still Shows Up

Status: conditional.

If profiling shows finalization or limb normalization still taking a significant share after scratch reuse and dispatch cleanup, then optimize the NTT output and carry-propagation stages next.

Possible targets:

- tighter carry-propagation loops
- fewer trimming passes
- pointer-based finalization helpers
- narrower cleanup paths for common coefficient counts

This is worth doing only if the benchmark data still points there after the higher-return work above.

## Recommended Order

1. Reuse scratch buffers in multiplication and Newton division.
2. Extend prepared and cached APIs for repeated workloads.
3. Add threshold auto-tuning driven by shape-aware benchmarks.
4. Revisit extreme skew multiplication only if a clear band remains.
5. Optimize finalization if profiling still identifies it as a hot spot.

## What Not To Do Yet

- Do not add new asymptotic multiplication algorithms unless a benchmarked production band exists.
- Do not lower thresholds based on isolated microbenchmarks that ignore dispatcher overhead.
- Do not add platform-specific assembly before the portable hot paths have been squeezed first.
- Do not broaden the public API until the repeated-workload use case is clear and measured.

## Suggested Verification

Use the existing benchmark and correctness harnesses after each change:

- `tests/mult_correctness.cpp`
- `tests/div_correctness.cpp`
- `tests/performance/dispatch_tuner.cpp`
- `tests/performance/multiplication_shape_bench.cpp`
- `tests/performance/division_shape_bench.cpp`

The benchmark data should answer one question for each change: did this reduce total work in a real production band, or did it only move cost around?
