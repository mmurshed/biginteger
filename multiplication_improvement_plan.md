# Multiplication Improvement Plan

## Diagnosis

BigMath multiplication is already strong at very large balanced sizes: `BENCHMARK.md` shows it beating GMP at 5M+ digits. The remaining gaps are concentrated in smaller scalar/Karatsuba bands, the NTT crossover band, and very large skewed multiplication where the full NTT path may do too much symmetric work.

The next useful work should be benchmark-driven. Avoid adding asymptotically attractive algorithms without a measured production band.

## 1. Add Shape-Focused Multiplication Benchmarks

Add a focused benchmark that measures direct algorithms and dispatcher behavior across balanced and skewed limb shapes.

Suggested cases:

```text
balanced: 64, 96, 128, 256, 512, 1024, 2048, 4096, 8192 limbs
skewed:   1024x64, 2048x128, 4096x256, 8192x512, 16384x1024
large:    10n/n, 20n/n, 50n/n
```

Measure:

- Classic
- Karatsuba
- NTT Goldilocks
- CRT NTT
- dispatcher

Use end-to-end `Multiply()` timings, not only direct microbenchmarks.

## 2. Retune Dispatch by Shape

Current dispatch is mostly threshold-based. Total limb count is not enough for skewed operands.

Use the shape benchmark to tune:

- Karatsuba vs NTT threshold
- CRT NTT threshold
- serial vs threaded threshold
- skew-specific fallback bands

Goal: keep NTT for balanced and moderate-skew cases, but avoid full-size NTT when a blockwise or Karatsuba path wins for high skew.

## 3. Prototype Blockwise Skewed Multiplication

For `A` much larger than `B`, split `A` into chunks near the size of `B`:

```text
A = A0 + A1 * B^k + A2 * B^(2k) + ...
A * B = sum(Ai * B) shifted
```

Benchmark this against one full NTT. This may help large skewed cases such as `5M x 500k` and `10M x 1M`, where BigMath falls behind GMP again.

Reject the prototype if repeated smaller multiplications lose to one large NTT.

## 4. Optimize NTT Finalization

`FinalizeBase2_64` shows up as a meaningful cost in NTT-heavy profiles. Improvements here affect multiplication, division, parse, and ToString.

Investigate:

- tighter carry-propagation loops
- fewer cleanup/trimming passes
- specialized paths for common coefficient counts
- pointer-based finalization
- separate microbenchmark for finalization

Expected gain is likely single-digit percent, but broad.

## 5. Tune CRT NTT and Thread Thresholds Separately

The ideal thresholds may differ for:

- Goldilocks vs CRT
- serial vs threaded
- balanced vs skewed operands

Recheck `BIGMATH_NTT_CRT_THRESHOLD` and thread activation thresholds using shape-aware benchmarks. Avoid assuming the balanced optimum is also best for skewed multiplication.

## 6. Investigate Faster NTT Kernel

The butterfly and modular multiplication loops dominate large products.

Potential work:

- stronger loop unrolling
- pointer-based butterfly loops
- ARM64-specific assembly kernel
- better scheduling around `MUL` / `UMULH`
- reduced loads/stores in pointwise multiply

This is high effort and platform-sensitive. Do it only after benchmarks confirm the NTT kernel remains the bottleneck.

## 7. Consider Prepared NTT Operands

For repeated multiplication by the same large operand, add an internal or public prepared-operand API that caches:

- coefficient split
- transform size
- forward NTT spectra
- CRT per-prime spectra

This does not help one-off `a * b`, but can help exponentiation, repeated products, and some division/string-conversion internals.

## 8. Audit Squaring Use

Squaring already has specialized implementations. Ensure callers use `Square(x)` instead of `Multiply(x, x)` where applicable.

Known good use:

- `Pow10(d)` even branch

Future audit targets:

- exponentiation if added
- benchmark/test helpers
- internal algorithm code that computes self-products

## Avoid For Now

- Do not add Toom-3 or Toom-5 to dispatch; both were measured slower in the production band.
- Do not lower the NTT threshold based only on direct microbenchmarks; dispatch overhead changes the result.
- Do not reattempt Mulders-style high multiplication unless it reduces actual NTT work.
- Do not implement Schonhage-Strassen or Furer-style multiplication for this input range.

## Recommended Order

1. Add shape-focused multiplication benchmark.
2. Retune dispatch from measured balanced and skewed shapes.
3. Prototype blockwise skewed multiplication.
4. Optimize `FinalizeBase2_64`.
5. Retune CRT/thread thresholds by shape.
6. Investigate NTT kernel improvements.
7. Consider prepared NTT operands for repeated-workload APIs.

## 2026-05-26 Implementation Findings

### Completed

- Added `tests/performance/multiplication_shape_bench.cpp`, a shape-aware benchmark that measures Classic, Karatsuba, NTT, dispatcher, and an experimental blockwise skew prototype across balanced and skewed limb ratios.
- Retuned default dispatch for very small high-skew shapes: when `minSize <= 64` and `maxSize >= 10 * minSize`, dispatch now uses Classic instead of Karatsuba.
- Simplified Base2_64 NTT finalization loops in both Goldilocks and CRT NTT paths by packing coefficients directly into output limbs instead of using slot/flush lambdas.

### Measured Wins

The clearest win is the skew dispatch guard. Before the retune, `20n/n` at `1280x64` limbs dispatched through Karatsuba and took about `0.1522 ms`, while Classic took `0.0821 ms`. After the retune, dispatcher time is about `0.0848 ms`. The `50n/n` at `3200x64` case improved from about `0.4469 ms` to roughly `0.21-0.25 ms`, depending on run noise.

The finalization-loop changes were correctness-neutral and small/mixed in end-to-end timings. They remain because they remove avoidable branch/lambda overhead in hot NTT output code without changing algorithm selection.

### Prototyped But Not Promoted

- Blockwise skew multiplication can beat one-shot Classic in narrow `min=64`, high-ratio cases inside the benchmark harness, but it regressed when promoted into production dispatch and loses badly once the smaller operand reaches 512+ limbs. It remains a benchmark prototype only.
- Barrett-style modular reduction for CRT NTT prime multiplication was tested and rejected. It was consistently slower than the compiler's `% P` code for these fixed 32-bit NTT primes.
- Prepared/reusable NTT operands were considered but not added. The benefit requires a repeated-same-operand API because transform size depends on both operands; one-off `Multiply(a, b)` cannot reuse a prepared spectrum safely without a larger API contract.
- Assembly was not attempted. The portable work did not leave a newly isolated bottleneck that justifies platform-specific code in this pass.
