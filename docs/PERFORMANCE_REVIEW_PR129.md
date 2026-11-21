# Performance review of PRs #129–#150

## Implemented changes and measured results

The follow-up uses quotient-only division for decimal truncation, shared carry and borrow kernels for whole-vector integer arithmetic, and a dedicated path for small self-products. Scalar division moves its outputs, and scalar modulo computes only the remainder. The Karatsuba recursion reads its two partial products from the output buffer. Burnikel–Ziegler (BZ) writes mid-sized correction products directly into its scratch space.

The table compares #150 with the combined changes on an Apple M1 Max using Release builds, full-width deterministic inputs and the paired process order described below. Times are microseconds per public call. For the scale reduction, the signed input has 192 integer digits and 10,000 fractional digits. Both builds truncate it to the same 192-digit integer. The 2:1 decimal division provides a control where the short-quotient shortcut is ineligible.

| Operation | #150 µs | Follow-up µs | Speedup |
|---|---:|---:|---:|
| Add, 1,024 + 1,024 limbs | 1.394 | 0.693 | 2.01× |
| Subtract, 1,024 − 1,024 limbs | 1.796 | 0.694 | 2.59× |
| Decimal divide DOWN, 1,034 ÷ 1,024 limbs | 7.988 | 0.347 | 23.00× |
| Decimal divide DOWN, 1,024 ÷ 512 limbs | 98.983 | 94.050 | 1.05× |
| Decimal SetScale DOWN, drop 10,000 places | 4.367 | 0.329 | 13.26× |
| Scalar remainder, 64 limbs | 0.489 | 0.311 | 1.57× |
| Self-product, 256 limbs | 16.698 | 11.891 | 1.40× |

**Where the square path ends.** On the same build, forcing the square dispatcher at 319 limbs improved bursts by 1.39× and isolated calls by 1.36×. At 384 limbs it improved isolated calls by 1.45× but made bursts 21% slower in throughput (24.12 → 30.46 µs per call). The default therefore sends identical operands of up to 319 limbs to the square dispatcher. Above that ceiling, multiplication keeps its existing warm NTT entry. When that entry chooses serial Karatsuba because the pool is cold, identical operands use Karatsuba squaring. Equal-valued, distinct operands follow ordinary multiplication. The sweep used 13 samples per process, with 48 products per burst or a 2 ms pause before each isolated call.

**Ordinary products and division.** Small differences moved between runs. In the initial combined run, distinct-operand multiplication ranged from 0.90× to 1.03× across five shapes. A dedicated repeat of the same inputs, taking the median of four processes per build in mirrored order, measured 1.01–1.04×. BZ-heavy QR ranged from 0.99× to 1.05× on the repeat. At 5,120 ÷ 512 limbs the repeat measured 876.4 → 840.4 µs. These copy reductions have a smaller effect than the public-operator changes.

### Matched division and the remaining gap

At **100,000 ÷ 10,000 decimal digits**, the canonical runner now times the same outputs on both sides and checks each against GMP outside the timer. Each first call and the best of three subsequent calls on the same operands are shown below. The second mode can reuse process state from the first.

| Operation | BigMath first ms | GMP first ms | BigMath warm ms | GMP warm ms | Warm BigMath / GMP |
|---|---:|---:|---:|---:|---:|
| Quotient only | 0.912 | 0.480 | 0.890 | 0.451 | 1.97× |
| Quotient and remainder | 0.891 | 0.469 | 0.864 | 0.459 | 1.88× |

Separate instrumentation of the eight BZ calls on those same operands attributed **79.6%** of elapsed BZ time to 2,040 correction products, **16.4%** to 1,088 schoolbook leaves, and **0.8%** to 744 addback iterations. The balance includes normalization, subtraction and dispatch. These shares include the instrumentation overhead. The correction products dominate the outstanding gap on this input.

The half-size reciprocal suggested in #144 now has an optional GMP-checked probe. Its two guarded top limbs support division in half-divisor-sized blocks. The probe measures reciprocal construction separately and includes full back-products and correction in the end-to-end QR timing. It beat the dispatcher in **5 of 30** normalized shapes, including **2 of 10** at ratio 2:1. At that ratio its median runtime was 1.15× the dispatch time. For 4,096 ÷ 2,048 limbs, construction fell from 437.0 to 139.6 µs while the complete operation took 521.4 µs versus 480.0 µs through dispatch and 501.4 µs in GMP. Production keeps the existing Newton algorithm. The earlier fused-MFA and narrow Toom-3 proposals remain documented in the historical review.

### Benchmark and correctness checks

- Both Release layouts passed all four CTest targets, including 339 unit cases with 64-bit limbs and 319 with 32-bit limbs. The 32-bit baseline had 16 failing tests that supplied noncanonical limbs. The repaired fixtures generate data for the active radix. AddressSanitizer and UndefinedBehaviorSanitizer also passed the default multiplication, division and unit tests.
- Regression coverage includes signed scalar quotients and remainders wider than 32 bits, decimal scale shifts in both directions, truncation to zero, identical and distinct operands around the square cutoff, and a calculator power. Both simple performance programs agreed with their arithmetic references.
- The split driver rebuilds through CMake. Changing a library source rebuilt BigMath and reused the GMP timings. Changing only a cached QR remainder hash triggered a correctness failure. Both benchmark harnesses report the requested division mode, first-call and subsequent timings, and repetition counts.
- The full tuner generated a 64-bit profile using the configured CRT square backend and full-width limbs. A separate build with that header passed CTest and checked 28 held-out division shapes against the schoolbook reference. Its median dispatch time across those shapes was 1.01× the generic profile's time, with individual ratios of 0.60–1.86×. This host profile stays experimental. The README and benchmark notes give the measured 12-case division speedup and the actual balanced NTT entry floors.

The following commands replay the division comparison, the optional reciprocal experiment and the held-out shape checks from a Release build with GMP installed. The original prototype replay and its complete data appear later in this document.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target bench_vs_gmp division_half_reciprocal_probe division_shape_bench -j8
./build/bench_vs_gmp --div 100000 10000
./build/division_half_reciprocal_probe
./build/division_shape_bench 101 777 1553 3091
```

## Historical findings at #150

**Review date: 2026-11-21.** The largest missed opportunity I measured was in a caller of the new division code. For a decimal division that rounds toward zero, BigMath still computes a remainder and discards it. Using the quotient-only path from PR #147 made a 1,034-by-1,024-limb case **7.86× faster** on the M1 Max. The broadest opportunity was simpler: using PR #136’s existing carry kernels for ordinary BigInteger addition and subtraction made the 1,024-limb public operators **1.95× and 2.51× faster**, respectively. These are comparisons of local prototypes against merged BigMath. Each describes a particular operation and operand shape.

## What the pass already delivered

I reviewed the 22 merged PRs numbered #129–#150, using the end of #128 (`7dfb87f`) as the baseline and #150 (`0e9b3c3`) as the final source. #130 merged before #129, and #137 before #136. The commits in this range carry Claude co-author trailers. PR #151 updates the figure and was open at review time. The performance experiments below use clean archives of #150, independent of that branch and the local division benchmark work. [Full source comparison](https://github.com/mmurshed/biginteger/compare/7dfb87f...0e9b3c3)

| PRs | Changes relevant to this review |
|---|---|
| #129–#132 | Moved multiplication results into their wrappers, retuned the schoolbook cutoff, fixed the Linux test build, and extended and simplified the NEON transform butterflies. |
| #133–#139 | Revised worker scheduling, parallelized CRT reconstruction and inverse work, added carry kernels, repaired skewed Karatsuba workspace handling, truncated selected transforms, and retuned multiplication by operand shape. |
| #140–#141 | Published the multiplication results, then added a second NTT entry threshold for calls that find the thread pool cold. |
| #142–#145 | Rebuilt the division basecase and Burnikel–Ziegler recursion, reduced Newton reciprocal work, and swept the division dispatcher over divisor size and operand ratio. |
| #146–#150 | Improved single-limb division, short quotient-only division, worker pickup, and Newton’s leading block, then published the division results and comparison caveats. |

The canonical tables show balanced multiplication crossing GMP at the 10k-digit sample and skewed division at the 200k-digit dividend sample on this machine. The latter comparison uses different GMP and BigMath output modes. PR #150 also reports a matched quotient-and-remainder comparison that wins at its 200k sample (0.77× GMP time). These are measured sizes on a threaded BigMath build. [Canonical tables and methodology at #150](https://github.com/mmurshed/biginteger/blob/0e9b3c3/BENCHMARK.md)

## Where another patch would help

Times below are microseconds per public call. A speedup above 1 means the prototype took less time. Operand sizes in this table are counts of **64-bit limbs**. Rows with similar or worse times are included to show where an optimization stops helping.

| Operation and shape | Merged | Prototype | Speedup |
|---|---:|---:|---:|
| BigInteger add, 1,024 + 1,024 limbs | 1.363 | 0.701 | 1.95× |
| BigInteger subtract, 1,024 − 1,024 limbs | 1.743 | 0.695 | 2.51× |
| BigDecimal divide, round DOWN, 522 ÷ 512 limbs | 4.173 | 0.757 | 5.51× |
| BigDecimal divide, round DOWN, 1,034 ÷ 1,024 limbs | 7.737 | 0.985 | 7.86× |
| BigDecimal divide, round DOWN, 1,024 ÷ 512 limbs | 97.135 | 97.292 | 1.00× |
| BigInteger multiplied by itself, 256 limbs | 16.671 | 11.641 | 1.43× |
| BigInteger multiplied by itself, 384 limbs | 26.081 | 31.882 | **0.82×** |
| Scalar remainder, 64-limb dividend | 0.489 | 0.350 | 1.40× |
| Scalar quotient, 4,096-limb dividend | 20.122 | 20.207 | 1.00× |

**1. Use quotient-only division for decimal truncation.** `BigDecimal::Divide` rescales the magnitudes, calls `DivideAndRemainder`, then consults the rounding mode. With `DOWN`, the rounding decision is always to keep the truncated quotient. The prototype used `num / den` after the existing scale adjustment, applied the result sign, and preserved the requested scale. That lets #147 divide using only the high limbs when the quotient is short. At 2:1 that shortcut is ineligible, which explains the flat control row. The benchmark used scale zero and all four input-sign combinations. Before shipping, cover nonzero scale shifts in both directions, exact multiples, a negative result that rounds to zero, and the existing divide-by-zero behavior. Reducing scale in `SetScale(..., DOWN)` has the same opportunity and needs a separate measurement. [Decimal implementation](../bigdecimal/BigDecimal.cpp) · [Quotient dispatcher](../src/algorithms/Division.cpp) · [PR #147](https://github.com/mmurshed/biginteger/pull/147)

**2. Reuse the carry kernels in whole-vector addition and subtraction.** #136 sped up the add/subtract passes *inside* Karatsuba. The public whole-vector operations still call a range routine that checks bounds on each limb. The prototype used `AddExt` for the full sum and `SubN` plus `SubBorrow` for the magnitude difference, preserving the existing result allocation, trimming, and public sign handling. At 4,096 limbs the gains were 2.04× and 2.59×. Start with the returning, full-vector overloads. The range overloads have distinct contracts for aliased operands, short destination windows and propagation into an existing destination. Test unequal lengths, opposite signs, long carry/borrow runs and both limb layouts before extending the change. [Addition](../src/algorithms/Addition.cpp) · [Subtraction](../src/algorithms/Subtraction.cpp) · [Existing kernels](../include/biginteger/common/LimbKernels.h)

**3. Select specialized squaring for small self-products.** The calculator evaluates a power by repeatedly executing `b = b * b`. That expression enters ordinary multiplication. Large CRT products already detect identical operand references and skip duplicate forward transforms. At smaller sizes a distinct square dispatcher can exploit symmetry earlier. An identity check on the public operands improved 16-, 64- and 256-limb self-products by 1.24×, 1.40× and 1.43×. Routing every size through `Square` regressed at 384 limbs in both paired runs, and the 4,096-limb gain was about 1%. Sweep the transition with isolated calls and bursts before selecting an upper bound. Check an actual exponentiation, equal values stored in different objects, negative bases and the cost of the identity check on ordinary products. [Public multiply](../src/ops/Multiplication.cpp) · [Square dispatcher](../src/algorithms/Squaring.cpp) · [Calculator caller](../calculator/ExpressionEvaluator.h)

**4. Finish result moves in scalar division.** #142 moved results on the BigInteger-divisor path. The `DataT` overload still constructs a quotient from an lvalue vector. Its quotient-and-remainder path also copies both vectors into a pair, both into BigIntegers, then both BigIntegers into another pair. Moving at those boundaries helped the 64-limb scalar remainder by 1.40×, falling to 1.08× at 4,096 limbs. The tested divisor was 1,000,000,007. A remainder-only kernel could also avoid allocating and writing an entire quotient for `a % d`. That is a separate, unmeasured change. Check the sign of negative remainders, scalar zero, and both radix layouts before shipping. [Scalar operators](../src/ops/ScalarDivision.cpp) · [Single-limb implementation](../include/biginteger/algorithms/division/ClassicDivision.h)

## Two plausible core changes with small returns

I tested each of these independently against #150 using the public product or quotient-and-remainder operation. Both directly target work still present after the pass.

| Isolated prototype | Example timings (merged → prototype, µs) | Decision |
|---|---|---|
| Karatsuba: subtract the two partial products directly from their existing output slices before recombination, eliminating two copies at each balanced recursive node. | 64 × 64: 1.712 → 1.722. 256 × 256: 16.625 → 16.262. 1,024 × 48: 23.601 → 23.054. | Roughly 0.99–1.03× across five product shapes. Defer as a standalone speed project. [Karatsuba recombination](../include/biginteger/algorithms/multiplication/KaratsubaMultiplication.h) |
| Burnikel–Ziegler: write mid-sized Karatsuba correction products directly into the division scratch buffer and reuse recursion workspace. This avoids staging both operands and copying the product back in `MulInto`. | 1,024 ÷ 512: 95.596 → 96.384. 5,120 ÷ 512: 861.465 → 839.477. 1,024 ÷ 256: 93.208 → 88.629. | About 0.96–1.05× across five division shapes. This change alone cannot explain or close the documented ~2× gap at 100k ÷ 10k decimal digits. [BZ product path](../include/biginteger/algorithms/division/BurnikelZieglerDivision.h) |

For the remaining **100k ÷ 10k decimal-digit** case, first collect a matched Q/Q and QR/QR timing at that exact shape. The nearby 5,120 ÷ 512 **limb** probe gained only 2.6% from bypassing BZ product staging. Time the schoolbook leaves, correction products and addbacks separately on the exact operands before choosing another kernel to rewrite.

Two larger ideas were already identified in the PRs. #138 left the large fused-MFA transform staircase after timing the adjacent size class. #144 proposed a half-size reciprocal and two quotient blocks near a 2:1 division. A pointer-based Toom-3 prototype had also been tried, with about a 4% average advantage in a narrow pre-NTT window. For further work specifically on core multiply/divide, the half-size reciprocal is a defined experiment: time construction and complete 2:1 QR separately, include adjacent ratios and powers of two ±1, and compare against the current dispatcher and GMP QR. Its gain remains unknown. [PR #138](https://github.com/mmurshed/biginteger/pull/138) · [PR #144](https://github.com/mmurshed/biginteger/pull/144) · [Toom-3 measurements](MULTIPLICATION.md#toom-3-on-the-limb-kernels--re-evaluated-2026-11-20-still-out)

## Fix the measurements before the next threshold sweep

- **The documented tuner times a different square backend.** `TuneSquaring` benchmarks `NTTSquare::Square` (Goldilocks), while the default dispatcher uses CRT self-multiplication. The tuner also fills only the low 32 bits of every limb and sets bit 31 of the top limb even in a 64-bit build. Its division sweep starts Newton measurements at 8,192 limbs and emits legacy division thresholds, while the default 64-bit dispatcher can enter Newton at 768 or 1,536. Generate inputs for the active radix, time the selected backend, measure warm and cold entry separately, and replay the emitted profile against held-out shapes on the target CPU. [Tuner](../tests/performance/dispatch_tuner.cpp) · [Active thresholds](../include/biginteger/build/DispatchThresholds.h) · [Tuning instructions](PLATFORM_TUNING.md)
- **The split benchmark can silently reuse an older BigMath executable.** Its shell driver compiles library sources into `run_bigmath` but decides whether to rebuild using only the runner source and `bench_common.h`. A change to a kernel, library source, workload header, or build flags alone can leave the old binary in use. Build through CMake or key the cache on all inputs. Verify by changing a library source or header between consecutive benchmark invocations. This affects the split driver. The canonical November tables identify the separate `bench_vs_gmp` executable as their source. [Split benchmark driver](../tests/performance/benchsuite/run_benchmark.sh)
- **Compare matched division operations in both harnesses.** #150 already documented that the canonical test times BigMath `/` against GMP `mpz_tdiv_qr`. The split runner repeats the mismatch and hashes only the quotient. Publish distinct Q/Q and QR/QR rows, verify each requested output outside the timer, and label first-call and steady-state measurements. The benchmark methodology also says up to 1,000+ iterations, while the integer harness uses 7, 5, 3 or 1 and caps ToString at one from 100k digits. Record the actual repetitions with each row. [Canonical runner](../tests/performance/bench_vs_gmp.cpp) · [BigMath split runner](../tests/performance/benchsuite/run_bigmath.cpp) · [GMP split runner](../tests/performance/benchsuite/run_gmp.cpp) · [GMP division APIs](https://gmplib.org/manual/Integer-Division)
- **Correct two published descriptions.** Dividing the November 20 BigMath times by November 21 times for all 12 corresponding skewed-division rows gives **1.14–1.58×**. The endpoints are 817.6 → 715.2 ms at 100M ÷ 20M and 87.825 → 55.673 ms at 10M ÷ 2M. The README’s upper bound of 3.2× needs a separately identified workload or baseline. Its multiplication feature list also gives the old 1,280-total-limb NTT entry. Under the default 64-bit rules, balanced operands first satisfy all warm entry conditions at **320 limbs each** (640 total), and reach the unconditional cold threshold at **512 each**. The area condition alone would suggest 310 each and misses the total-size floor. [November 20 tables](https://github.com/mmurshed/biginteger/blob/40599dd/BENCHMARK.md) · [November 21 tables at #150](https://github.com/mmurshed/biginteger/blob/0e9b3c3/BENCHMARK.md) · [README](../README.md) · [Multiply dispatcher](../src/algorithms/Multiplication.cpp)

## How the new numbers were obtained

Apple M1 Max, Apple Clang 21.0.0, CMake Release, `-O3 -march=native -DNDEBUG`, default 64-bit limbs and thread pool. The four public-operator edits shared a prototype build. The Karatsuba and BZ experiments each changed one kernel in a separate build. Each executable used the same four deterministic full-width operand pairs as its comparison build at each shape (seed 129, top bit set). Each timed batch consumed the low result limb, included result construction and destruction, and excluded operand generation and GMP checks. The decimal timing also copied the unscaled quotient out of the returned decimal on both sides. I took the median of seven batches in each process, ran merged / prototype / prototype / merged, and report the mean of the two process medians per build. The two independent core variants ran in mirrored order around the baseline. These warm local timings do not estimate cold latency or performance on another CPU.

All timed inputs passed exact GMP value checks. Division probes checked the requested quotient or both quotient and remainder, and decimal checks also verified the scale. Selected repository tests passed on the corresponding prototypes: 97 for the public-operator changes, 27 for Karatsuba, and 4 for BZ/dispatch. These counts cover the selected suites. The experiments used one machine, and targeted random shapes do not establish production readiness across limb layouts or prove a library-wide speedup.

I also ran the replay script below from fresh source archives. All GMP checks passed again. It measured **1.99×** for 1,024-limb addition, **2.56×** for subtraction, and **8.34×** for the 1,034 ÷ 1,024 decimal case. The 384-limb self-square was close to even on replay (0.98× versus 0.82× originally), reinforcing the need to measure that transition separately.

<details>
<summary>Original timings for every measured shape (µs per call)</summary>

Each cell gives the two process medians in execution order for that build. The speedups in the review use their arithmetic means.

| Operation | Limbs A, B | Merged | Public-operator prototype |
|---|---|---:|---:|
| add | 16, 16 | 0.04348, 0.04321 | 0.03833, 0.03094 |
| sub | 16, 16 | 0.05100, 0.05128 | 0.03756, 0.03171 |
| scalar-q | 16, 1 | 0.17650, 0.10776 | 0.08510, 0.09168 |
| scalar-r | 16, 1 | 0.19697, 0.19460 | 0.10624, 0.10378 |
| self-mul | 16, 16 | 0.17372, 0.17517 | 0.14217, 0.13994 |
| add | 64, 64 | 0.12415, 0.11940 | 0.06737, 0.06812 |
| sub | 64, 64 | 0.21842, 0.14945 | 0.06992, 0.07003 |
| scalar-q | 64, 1 | 0.37283, 0.36577 | 0.34287, 0.34075 |
| scalar-r | 64, 1 | 0.49050, 0.48792 | 0.34710, 0.35298 |
| self-mul | 64, 64 | 1.68449, 1.68405 | 1.20655, 1.20165 |
| add | 256, 256 | 0.35400, 0.37040 | 0.19107, 0.19140 |
| sub | 256, 256 | 0.46547, 0.46620 | 0.18620, 0.19253 |
| scalar-q | 256, 1 | 1.31560, 1.33387 | 1.25820, 1.24713 |
| scalar-r | 256, 1 | 1.49327, 1.47273 | 1.27480, 1.27193 |
| self-mul | 256, 256 | 16.87320, 16.46920 | 11.63767, 11.64402 |
| add | 384, 384 | 0.51002, 0.50761 | 0.26763, 0.26673 |
| sub | 384, 384 | 0.68059, 0.67809 | 0.26853, 0.26973 |
| scalar-q | 384, 1 | 1.95002, 1.93740 | 1.86218, 1.92168 |
| scalar-r | 384, 1 | 2.15335, 2.16036 | 1.89223, 1.89213 |
| self-mul | 384, 384 | 26.87365, 25.28765 | 35.58065, 28.18413 |
| add | 1024, 1024 | 1.37046, 1.35550 | 0.68189, 0.71955 |
| sub | 1024, 1024 | 1.74279, 1.74279 | 0.68029, 0.70887 |
| scalar-q | 1024, 1 | 5.24733, 5.14503 | 5.01683, 4.95833 |
| scalar-r | 1024, 1 | 5.60443, 5.44204 | 4.89957, 5.00053 |
| self-mul | 1024, 1024 | 39.29944, 45.52869 | 45.76562, 44.29688 |
| add | 4096, 4096 | 5.26250, 5.48958 | 2.61875, 2.66040 |
| sub | 4096, 4096 | 6.81980, 6.82083 | 2.61875, 2.65417 |
| scalar-q | 4096, 1 | 20.28855, 19.95520 | 20.22918, 20.18437 |
| scalar-r | 4096, 1 | 22.43645, 21.24790 | 20.13543, 20.28125 |
| self-mul | 4096, 4096 | 155.69531, 155.50262 | 154.51044, 154.05469 |
| down | 512, 512 | 1.04297, 1.00911 | 0.56445, 0.53189 |
| down | 522, 512 | 4.09506, 4.25131 | 0.71614, 0.79817 |
| down | 1034, 1024 | 7.76302, 7.71094 | 0.93294, 1.03645 |
| down | 4106, 4096 | 30.78842, 31.62175 | 4.98112, 4.84766 |
| down | 1024, 512 | 96.89714, 97.37239 | 97.48959, 97.09375 |

| Operation | Limbs A, B | Merged | Karatsuba prototype | BZ prototype |
|---|---|---:|---:|---:|
| mul | 64, 64 | 1.71311, 1.71050 | 1.68707, 1.75608 | 1.67839, 1.77560 |
| mul | 129, 127 | 5.47743, 5.51302 | 5.33854, 5.55816 | 5.44098, 5.67751 |
| mul | 256, 256 | 16.44661, 16.80382 | 15.98351, 16.53950 | 17.13889, 16.18229 |
| mul | 256, 63 | 7.21136, 7.23741 | 6.95703, 7.04383 | 7.01867, 6.87674 |
| mul | 1024, 48 | 23.42014, 23.78255 | 22.66016, 23.44878 | 23.36110, 24.04904 |
| qr | 256, 128 | 9.75000, 9.80208 | 10.11458, 10.01383 | 9.54517, 9.93058 |
| qr | 1024, 512 | 95.52433, 95.66667 | 94.41325, 96.21183 | 96.15975, 96.60767 |
| qr | 5120, 512 | 860.37500, 862.55558 | 844.83333, 844.91317 | 839.35075, 839.60417 |
| qr | 1538, 769 | 187.11808, 191.90275 | 200.47567, 192.61117 | 207.21525, 188.00700 |
| qr | 1024, 256 | 92.65625, 93.76042 | 90.57983, 90.57992 | 87.67367, 89.58333 |

</details>

<details>
<summary>Reproduce the prototype builds, GMP checks and timing batches</summary>

Save the Python block below as a file and run it from the repository root. It builds four source trees in a new temporary directory. Use a C++20 compiler and an installed GMP development library. Set `GMP_PREFIX` if GMP is outside the compiler's default search path. On the review machine the matching Xcode SDK was selected with `SDKROOT=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk`. The fixed thresholds in the BZ prototype restrict the experiment to the measured default 64-bit Karatsuba band.

```python
"""Replay PR #129–150 audit from the repository root. Requires Python 3, CMake, C++20 and GMP."""
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

repo = Path.cwd()
work = Path(tempfile.mkdtemp(prefix="bigmath-perf-review-"))
print("Writing sources and CSV files to", work, flush=True)
os.chdir(work)
archive = work / "source.tar"
subprocess.run(["git", "-C", str(repo), "archive", "-o", str(archive),
                "0e9b3c3e44c506348b74b1bed3bcfe51bc12617e"], check=True)
Path("src").mkdir()
subprocess.run(["tar", "-xf", str(archive), "-C", "src"], check=True)
shutil.copytree("src", "alt")

def replace(file,old,new):
    p=Path("alt")/file
    data=p.read_text()
    assert data.count(old)==1, file
    p.write_text(data.replace(old,new))
add="src/algorithms/Addition.cpp"
sub="src/algorithms/Subtraction.cpp"
for f in (add,sub):
    replace(f,'#include <algorithm>','#include <algorithm>\n#include "biginteger/common/LimbKernels.h"')
replace(add, '    result[size - 1] = 0;', '''    if (base == Base2_64 && !a.empty() && !b.empty()) {
      LimbKernel::AddExt(result.data(), size, a.data(), a.size(), b.data(), b.size());
      TrimZeros(result);
      return result;
    }
    result[size - 1] = 0;''')
replace(sub, '    std::vector<DataT> result(size);', '''    std::vector<DataT> result(size);
    if (base == Base2_64 && !b.empty() && a.size() >= b.size()) {
      DataT bo = LimbKernel::SubN(result.data(), a.data(), b.data(), b.size());
      LimbKernel::SubBorrow(result.data()+b.size(), a.data()+b.size(), a.size()-b.size(), bo);
      TrimZeros(result);
      return result;
    }''')
f="bigdecimal/BigDecimal.cpp"
replace(f,'    auto [q, r] = DivideAndRemainder(num, den);','''    if (mode == RoundingMode::DOWN) {
      auto q = num / den;
      q.SetSign(resultNeg);
      return BigDecimal(std::move(q), newScale);
    }
    auto [q, r] = DivideAndRemainder(num, den);''')
f="src/ops/Multiplication.cpp"
replace(f,'#include <utility>','#include <utility>\n#include "biginteger/algorithms/Squaring.h"')
replace(f,'    auto result = Multiply(','''    if (&a == &b) return BigInteger(Square(a.GetInteger(), BigInteger::Base()), false);
    auto result = Multiply(''')
f="src/ops/ScalarDivision.cpp"
for a,b in [('result(q, false)','result(std::move(q), false)'),('q(result.first, false)','q(std::move(result.first), false)'),('r(result.second, false)','r(std::move(result.second), false)'),('return {q, r};','return {std::move(q), std::move(r)};')]: replace(f,a,b)
replace('include/biginteger/algorithms/division/ClassicDivision.h','return {w, v};','return {std::move(w), std::move(v)};')

for variant in ['kara','bz']:
    shutil.copytree('src',variant)
    kara=Path(variant+'/include/biginteger/algorithms/multiplication/KaratsubaMultiplication.h')
    data=kara.read_text()
    if variant=='kara':
        begin=data.index('            // Now safely copy T1 and T2')
        end=data.index('            // c[m..',begin)
        data=data[:begin]+'''            SubtractFromPtr(t3, lenT3, c, lenT1, base);
            SubtractFromPtr(t3, lenT3, c + 2 * m, lenT2, base);

'''+data[end:]
    else:
        start=data.index('    public:\n')
        data=data[:start]+'''    public:
        static void ProbeMultiplyInto(const DataT* a, SizeT la, const DataT* b, SizeT lb, DataT* dest)
        {
            static thread_local std::vector<DataT> scratch;
            size_t need = 16 * (size_t)std::max(la, lb) + 1024;
            if (scratch.size() < need) scratch.resize(need);
            MultiplyRecursive(a, la, b, lb, dest, scratch.data(), Base2_64);
        }
'''+data[start:]
        bz=Path('bz/include/biginteger/algorithms/division/BurnikelZieglerDivision.h')
        s=bz.read_text()
        loc=s.index('      static thread_local std::vector<DataT> av, bv;')
        bz.write_text(s[:loc]+'''      if (la + lb > 56 && la + lb < 640)
      {
        KaratsubaMultiplication::ProbeMultiplyInto(a, la, b, lb, tp);
        return;
      }
'''+s[loc:])
    kara.write_text(data)

common = r'''#include "biginteger/ops/Operations.h"
#include "biginteger/common/Builder.h"
#include "bigdecimal/BigDecimal.h"
#include <gmp.h>
#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
using namespace BigMath;
using namespace std;

static mt19937_64 rng(129);
static volatile DataT sink;
static BigInteger operand(size_t n) {
    vector<DataT> v(n);
    for (auto &x : v) x = rng();
    v.back() |= DataT(1) << 63;
    return BigInteger(move(v), false);
}
static void import_int(mpz_t z, const BigInteger &x) {
    mpz_import(z, x.size(), -1, sizeof(DataT), 0, 0, x.GetInteger().data());
    if (x.IsNegative()) mpz_neg(z, z);
}
static void check(const BigInteger &x, mpz_srcptr expected) {
    mpz_t got;
    mpz_init(got);
    import_int(got, x);
    if (mpz_cmp(got, expected)) abort();
    mpz_clear(got);
}
template<class F> void measure(const char *op, int n, int m, int reps, F fn) {
    using clock = chrono::steady_clock;
    array<double,7> us;
    for (int i=0; i<4; ++i) { auto v=fn(i); sink=v[0]; }
    for (double &sample : us) {
        auto start = clock::now();
        for(int i=0;i<reps;++i) {auto v=fn(i%4); sink=v[0];}
        sample=chrono::duration<double,micro>(clock::now()-start).count()/reps;
    }
    sort(us.begin(),us.end());
    printf("%s,%d,%d,%.5f\n",op,n,m,us[3]);
}
'''
Path('probe.cpp').write_text(common + r'''int main(){
  mpz_t a,b,result;
  mpz_inits(a,b,result,nullptr);
  constexpr DataT scalar=1000000007;
  for(int n : {16,64,256,384,1024,4096}) {
    array<BigInteger,4> x,y;
    for(int i=0;i<4;++i){
      x[i]=operand(n); y[i]=operand(n);
      if(x[i].CompareTo(y[i])<0) swap(x[i],y[i]);
      import_int(a,x[i]); import_int(b,y[i]);
      mpz_add(result,a,b);check(x[i]+y[i],result);
      mpz_sub(result,a,b);check(x[i]-y[i],result);
      mpz_mul(result,a,a);check(x[i]*x[i],result);
      mpz_tdiv_q_ui(result,a,scalar);check(x[i]/scalar,result);
      mpz_tdiv_r_ui(result,a,scalar);check(x[i]%scalar,result);
    }
    int reps=max(40,160000/n);
    measure("add",n,n,reps,[&](int i){return x[i]+y[i];});
    measure("sub",n,n,reps,[&](int i){return x[i]-y[i];});
    measure("scalar-q",n,1,reps,[&](int i){return x[i]/scalar;});
    measure("scalar-r",n,1,reps,[&](int i){return x[i]%scalar;});
    measure("self-mul",n,n,max(16,12000/n),[&](int i){return x[i]*x[i];});
  }
  for(auto [n,m]: {pair{512,512},pair{522,512},pair{1034,1024},pair{4106,4096},pair{1024,512}}) {
    array<BigDecimal,4> x,y;
    for(int i=0;i<4;++i) {
       auto v=operand(n), w=operand(m);
       if(n==m && v.CompareTo(w)<0) swap(v,w);
       v.SetSign(i&1); w.SetSign(i&2);
       import_int(a,v); import_int(b,w);
       x[i]=BigDecimal(v,0);y[i]=BigDecimal(w,0);
       mpz_tdiv_q(result,a,b);
       auto z=x[i].Divide(y[i],0,RoundingMode::DOWN);
       if(z.Scale()!=0) abort();
       check(z.UnscaledValue(),result);
    }
    measure("down",n,m,64,[&](int i){
       auto z=x[i].Divide(y[i],0,RoundingMode::DOWN);
       return z.UnscaledValue();
    });
  }
  mpz_clears(a,b,result,nullptr);
}
''')
Path('core.cpp').write_text(common + r'''
int main() {
    mpz_t a,b,ref,rem;
    mpz_inits(a,b,ref,rem,nullptr);
    for(auto [n,m]: {pair{64,64},pair{129,127},pair{256,256},pair{256,63},pair{1024,48}}) {
        array<BigInteger,4> x,y;
        for(int i=0;i<4;++i) {
            x[i]=operand(n); y[i]=operand(m);
            import_int(a,x[i]); import_int(b,y[i]);
            mpz_mul(ref,a,b); check(x[i]*y[i],ref);
        }
        measure("mul",n,m,96,[&](int i){return x[i]*y[i];});
    }
    for(auto [n,m]: {pair{256,128},pair{1024,512},pair{5120,512},pair{1538,769},pair{1024,256}}) {
        array<BigInteger,4> x,y;
        for(int i=0;i<4;++i) {
            x[i]=operand(n); y[i]=operand(m);
            import_int(a,x[i]); import_int(b,y[i]);
            mpz_tdiv_qr(ref,rem,a,b);
            auto qr=DivideAndRemainder(x[i],y[i]);
            check(qr.first,ref); check(qr.second,rem);
        }
        measure("qr",n,m,12,[&](int i){auto qr=DivideAndRemainder(x[i],y[i]);return std::move(qr.first);});
    }
    mpz_clears(a,b,ref,rem,nullptr);
}
''')

sdk = os.environ.get("SDKROOT")
cmake_flags = ["-DCMAKE_OSX_SYSROOT=" + sdk] if sdk else []
compile_flags = ["-isysroot", sdk] if sdk else []
gmp = os.environ.get("GMP_PREFIX")
if gmp:
    compile_flags += ["-I" + gmp + "/include", "-L" + gmp + "/lib"]

for tree, build, drivers in [
    ("src", "build", [("probe", "before"), ("core", "ref-core")]),
    ("alt", "candidate", [("probe", "after")]),
    ("kara", "karabuild", [("core", "kara-core")]),
    ("bz", "bzbuild", [("core", "bz-core")]),
]:
    subprocess.run(["cmake", "-S", tree, "-B", build, "-DCMAKE_BUILD_TYPE=Release",
                    "-DBIGMATH_BUILD_TESTS=OFF", *cmake_flags], check=True)
    subprocess.run(["cmake", "--build", build, "-j8"], check=True)
    for source, exe in drivers:
        subprocess.run([os.environ.get("CXX", "c++"), "-std=c++20", "-O3", "-march=native",
                        "-DNDEBUG", *compile_flags, "-I" + tree + "/include", source + ".cpp",
                        build + "/libbigmath.a", "-lgmp", "-pthread", "-o", exe], check=True)

for sequence in [
    ["before", "after", "after", "before"],
    ["ref-core", "kara-core", "bz-core", "bz-core", "kara-core", "ref-core"],
]:
    for i, exe in enumerate(sequence):
        with open(f"{exe}-{i}.csv", "w") as out:
            subprocess.run(["./" + exe], stdout=out, check=True)
print("CSV columns: operation, limbs A, limbs B (1 for scalar), median microseconds")
print("Results:", work)

```

</details>
