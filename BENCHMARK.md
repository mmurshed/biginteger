# BigMath vs GMP — full benchmark suite

Measured 2026-05-30 (fresh local rerun). Single-run snapshot capturing the state of `master` after the 2026-05 optimization pass (LIMB_64 default, multi-prime CRT NTT default, multithreaded NTT default, M-G div2by1 in `ClassicDivision`, M-G 3/2 qhat in `FastDivision` Base2_64, BZ Knuth-normalize fix, radix-4 + radix-8 fused NTT butterflies (PRs #59, #60), **Matrix Fourier Algorithm / Bailey 6-step CRT NTT retuned to n ≥ 2^24 coefficients**).

For per-subsystem deep dives — algorithms, dispatch, optimization history, rejected approaches — see:

- [docs/MULTIPLICATION.md](docs/MULTIPLICATION.md)
- [docs/DIVISION.md](docs/DIVISION.md)
- [docs/STRING_CONVERSION.md](docs/STRING_CONVERSION.md)
- [docs/BASE.md](docs/BASE.md)

This file is the **raw numbers**, all five subsystems in one place, kept for headline-level reference.

---

## Methodology

- **Hardware:** Apple M1 Max (10-core, 8 performance + 2 efficiency). 64 GB RAM.
- **OS:** macOS (Darwin 25.5.0).
- **Compilers:** Apple Clang via `c++ -std=c++20 -O3 -march=native`.
- **Reference:** GMP 6.3.0 (Homebrew).
- **Default stack:** `BIGMATH_LIMB_64=1`, `BIGMATH_NTT_CRT=1`, `BIGMATH_USE_THREADS=1` (8-thread pool, `min(hardware_concurrency, BIGMATH_MAX_THREADS)`).
- **Timing:** `min` over iteration counts ranging 1–20 depending on operand size. Sub-millisecond cases use 1000+ iters; multi-second cases use 1 iter.
- **Random inputs:** mt19937_64 seeded per case, full decimal digit range.

Harnesses:

- `tests/performance/bench_vs_gmp.cpp` — canonical bench. Mul covers ≤200M digits (balanced + skewed), div ≤200M×40M skewed, parse ≤50M, ToString ≤20M, plus BigDecimal arithmetic + I/O. Single-iter timing for digits ≥100k; multi-iter best-of below that. Per-row `[op size] setup...` log to stderr so long setup phases don't look stuck.

---

## Multiplication

Balanced (`a.size() == b.size()`):

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 1 000 × 1 000 | 0.002 | 0.001 | 2.27× |
| 5 000 × 5 000 | 0.030 | 0.011 | 2.64× |
| 10 000 × 10 000 | 0.101 | 0.046 | 2.20× |
| 50 000 × 50 000 | 0.698 | 0.398 | 1.76× |
| 100 000 × 100 000 | 1.144 | 0.709 | 1.61× |
| 500 000 × 500 000 | 5.120 | 4.322 | 1.18× |
| 1 000 000 × 1 000 000 | 10.176 | 8.845 | 1.15× |
| **2 000 000 × 2 000 000** | **21.570** | **20.994** | **1.03×** ← near parity |
| **5 000 000 × 5 000 000** | **48.947** | **65.245** | **0.75×** ← BigMath faster |
| **10 000 000 × 10 000 000** | **106.856** | **209.888** | **0.51×** ← BigMath 1.96× faster |
| **20 000 000 × 20 000 000** | **275.252** | **278.298** | **0.99×** ← parity |
| **50 000 000 × 50 000 000** | **1 248.374** | **674.945** | **1.85×** ← GMP faster |
| **100 000 000 × 100 000 000** | **2 617.801** | **1 453.331** | **1.80×** ← GMP faster |
| **200 000 000 × 200 000 000** | **6 578.614** | **3 049.051** | **2.16×** ← GMP faster |

Skewed (`a.size() >> b.size()`):

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 100 000 × 10 000 | 0.640 | 0.303 | 2.11× |
| 500 000 × 50 000 | 2.168 | 2.135 | 1.02× |
| **1 000 000 × 100 000** | **4.541** | **4.576** | **0.99×** ← BigMath faster |
| **2 000 000 × 200 000** | **9.238** | **9.479** | **0.97×** ← BigMath faster |
| 5 000 000 × 500 000 | 39.149 | 30.829 | 1.27× |
| 10 000 000 × 1 000 000 | 90.079 | 70.716 | 1.27× |
| 20 000 000 × 2 000 000 | 239.615 | 162.845 | 1.47× |
| **50 000 000 × 5 000 000** | **615.801** | **684.121** | **0.90×** ← BigMath faster |
| 100 000 000 × 10 000 000 | 1 159.229 | 1 011.270 | 1.15× |
| 200 000 000 × 20 000 000 | 2 446.578 | 1 952.674 | 1.25× |

**Observations:**

- **BigMath beats GMP on balanced multiplication across the 5M-10M band and is near parity at 20M.** Radix-4 + radix-8 fused NTT butterflies (PRs #59, #60) added 1.5-1.6× wall-clock vs prior. The 2026-05-27 MFA retune moved the gate from `2^21` to `2^24`, and the current 10M balanced row sits at **106.856 ms**.
- **MFA / Bailey 6-step CRT NTT (PR #65) is now reserved for the very-large regime.** The default gate is `2^24` transform coefficients. Focused limb benchmarks show this avoids the 300k-2M limb regression band while preserving MFA wins at 3M+ limbs.
- Below 500k, GMP's hand-tuned basecase keeps a 1.5-3.3× lead.
- **Skewed mults: BigMath is around parity at 500k×50k and 1M×100k, with a slight BigMath lead at 2M×200k.** BigMath falls back behind GMP at 50M×5M and 100M×10M, and is 1.25× at the new 200M×20M row.

### MFA focused threshold check

`mul_xl_bench` was run after raising the default gate to `2^24`. These are limb counts, not decimal digits. For Base2_64 balanced multiplication, `L` limbs per operand map to transform length `bit_ceil(4L - 1)`.

| limbs per operand | transform length | active path | ms |
|---:|---:|---|---:|
| 200 000 | 2^20 | non-MFA | 39.942 |
| 300 000 | 2^21 | non-MFA | 90.660 |
| 500 000 | 2^21 | non-MFA | 92.290 |
| 1 000 000 | 2^22 | non-MFA | 255.832 |
| 2 000 000 | 2^23 | non-MFA | 704.094 |
| **3 000 000** | **2^24** | **MFA** | **1 269.304** |
| **5 000 000** | **2^25** | **MFA** | **2 620.899** |

The retuned gate keeps the old MFA-off path through `2^23` and enables MFA at `2^24+`, matching the measured break-even region.

### Multithreaded NTT check

`BIGMATH_USE_THREADS=1` is the default. A focused `mul_xl_bench` run on 2026-05-27 compared the default build against `-DBIGMATH_USE_THREADS=0`:

| limb size | serial ms | threaded ms | speedup |
|---:|---:|---:|---:|
| 100 000 | 55.171 | 17.634 | 3.13× |
| 500 000 | 257.394 | 90.207 | 2.85× |
| 1 000 000 | 602.433 | 250.242 | 2.41× |

The threaded path uses coarse CRT parallelism for the six forward transforms and three inverse transforms, plus chunked pointwise multiplication. Single-prime Goldilocks fallback also uses size-gated layer parallelism.

### Prepared CRT NTT operands

`tests/performance/prepared_ntt_bench.cpp` measures repeated multiplication by one fixed large operand. The prepared API caches that operand's CRT spectra and reuses it across partner operands.

| build | limb shape | count | normal total | prepared total | setup | steady speedup | amortized speedup |
|---|---|---:|---:|---:|---:|---:|---:|
| threaded default | `100000x100000` | 5 | 95.100 ms | 87.037 ms | 6.921 ms | 1.09× | 1.01× |
| threaded default | `100000x100000` | 20 | 372.320 ms | 346.550 ms | 6.996 ms | 1.07× | 1.05× |
| serial `-DBIGMATH_USE_THREADS=0` | `100000x100000` | 5 | 291.137 ms | 198.299 ms | 16.590 ms | 1.47× | 1.35× |

The default threaded gain is modest because regular CRT multiplication already overlaps the six forward transforms across the thread pool. The prepared API is still useful for repeated workloads, especially serial builds or CPU-budget-sensitive callers.

### Public prepared/cached wrappers

The `ops/` layer now exposes the same reuse model through public wrappers.

| wrapper | repeated workload | normal total | prepared/cached total | speedup |
|---|---|---:|---:|---:|
| `PreparedMultiplication` | `100000` fixed limbs, `500000` partner limbs, 5 calls | 513.957 ms | 515.918 ms | 0.996× |
| `CachedDivision` | `100000`-limb divisor, `500000`-limb dividends, 5 calls | 2 419.124 ms | 1 374.231 ms | 1.76× |

`PreparedMultiplication` was flat on this shape because the repeated work was already dominated by the partner operand and the one-time preparation did not amortize over only five calls. `CachedDivision` showed the expected win on repeated skewed division, where the reciprocal setup is reused across the calls.

### Shape-focused multiplication dispatch

`tests/performance/multiplication_shape_bench.cpp` was added to measure direct Classic, Karatsuba, NTT, dispatcher, and an experimental blockwise-skew prototype by limb shape. It found one production dispatch miss: tiny high-skew operands should not enter Karatsuba.

| limb shape | previous dispatch | retuned dispatch | impact |
|---|---:|---:|---|
| `1280x64` (`20n/n`) | 0.1522 ms | 0.0848 ms | 1.8× faster |
| `3200x64` (`50n/n`) | 0.4469 ms | ~0.21-0.25 ms | ~1.8-2.1× faster |

The blockwise-skew prototype wins some `min=64` microbenchmarks but loses at larger smaller-operand sizes, so it remains benchmark-only. A Barrett modular multiply experiment for CRT NTT primes was slower than the existing `% P` lowering and was rejected.

### Toom-3 dispatch band (2026-05-27)

Focused band scan around the Karatsuba → NTT crossover found a narrow Toom-3 window:

| total limbs | per-operand | kara ms | toom3 ms | ntt ms |
|---:|---:|---:|---:|---:|
| 2 560 | 1 280 | 0.42 | 0.41 | 0.77 |
| 3 072 | 1 536 | 0.55 | 0.55 | 0.78 |
| 3 584 | 1 792 | 0.64 | 0.60 | 0.78 |
| 4 096 | 2 048 | 0.83 | 0.75 | 0.78 |
| **4 608** | **2 304** | **1.06** | **1.10** | **1.65** ← NTT-length boundary regression |
| 5 120 | 2 560 | 1.31 | 1.24 | 0.59 |

Toom-3 ties Karatsuba below total 3 584 and beats both Karatsuba and NTT in [3 584, 4 608]. At total 4 608 NTT pays a length-boundary penalty (next pow-2 NTT size is wasteful here); Toom-3 sidesteps it cleanly. Dispatch now uses Toom-3 for total ∈ [2 560, 5 120), bumping `NTT_MULTIPLICATION_THRESHOLD` from 4 096 → 5 120. The total-4 608 case improved from NTT 1.65 ms → Toom-3 1.10 ms (33% faster) end-to-end through the `Multiply()` entry point.

Toom-5 was checked in the same scan and never wins decisively — it ties Karatsuba in 64-256 limbs per-operand and degrades sharply above 256. It stays excluded from dispatch.

### Toom-3 skew gate (2026-05-30)

A follow-up focused scan found that the Toom-3 band above is a **balanced-product** win only. For 2:1-or-more skewed operands inside the same total-limb window, Toom-3 pays evaluation/interpolation overhead but does not get enough balanced subproblem work back. Dispatch now keeps Toom-3 only when `maxSize < 2 * minSize`; otherwise the pre-NTT window falls back to Karatsuba.

Focused limb benchmark, AppleClang `-O3 -march=native`, default `BIGMATH_LIMB_64=1`:

| limb shape | total limbs | Karatsuba ms | Toom-3 ms | NTT ms | old dispatch ms | new dispatch ms |
|---|---:|---:|---:|---:|---:|---:|
| `2560x128` | 2 688 | 0.279 | 0.607 | 0.776 | 0.603 | 0.280 |
| `2560x256` | 2 816 | 0.341 | 0.733 | 0.778 | 0.730 | 0.343 |
| `2048x1024` | 3 072 | 0.556 | 0.605 | 0.780 | 0.601 | 0.557 |
| `3072x1024` | 4 096 | 0.843 | 1.171 | 0.782 | 1.162 | 0.844 |
| `3072x1536` | 4 608 | 1.078 | 1.185 | 1.705 | 1.180 | 1.113 |

Balanced rows remain in the original Toom-3/NTT decision band. The scan still supports keeping `NTT_MULTIPLICATION_THRESHOLD=5120`: NTT wins at total 5120+, but regresses around total 4608 due to the transform-length boundary.

### Scratch-buffer reuse check (2026-05-30)

The hot-path scratch reuse pass keeps temporary buffers alive across calls in the CRT NTT MFA path and in Newton division. It was benchmarked with the existing dispatch and focused shape harnesses, but it did not move the production thresholds or produce a stable end-to-end band shift on this machine. The useful effect is lower allocation churn, not a headline crossover change.

Representative reruns after the change:

| shape | metric | ms |
|---|---|---:|
| `4096x8192` multiplication | dispatch | 0.606 |
| `8192x16384` multiplication | dispatch | 1.115 |
| `16000x32000` multiplication | dispatch | 2.204 |
| `32768x8192` division | Newton | 23.477 |
| `65536x16384` division | Newton | 50.712 |

These sit in the same band as the pre-change runs, so the doc-level takeaway is that scratch reuse is a cleanup win rather than a dispatch or asymptotic improvement.

---

## Division

Balanced (`a.size() == b.size()`) — quotient is 1-2 limbs, both libraries short-circuit. Numbers are sub-microsecond noise; reported for completeness.

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 1 000 × 1 000 | <0.001 | <0.001 | 10.90× |
| 5 000 × 5 000 | 0.001 | <0.001 | 6.41× |
| 10 000 × 10 000 | 0.002 | <0.001 | 6.55× |
| 50 000 × 50 000 | <0.001 | <0.001 | 0.43× |
| 100 000 × 100 000 | <0.001 | 0.001 | 0.23× |
| 500 000 × 500 000 | 0.001 | 0.005 | 0.18× |
| 1 000 000 × 1 000 000 | 0.001 | 0.010 | 0.10× |
| 5 000 000 × 5 000 000 | 1.367 | 0.189 | 7.23× |

Skewed (`a.size() >> b.size()`) — Newton/BZ band, real algorithmic work:

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 40 000 × 10 000 | 0.732 | 0.218 | 3.36× |
| 100 000 × 10 000 | 2.178 | 0.450 | 4.84× |
| 200 000 × 50 000 | 11.088 | 1.704 | 6.51× |
| 500 000 × 100 000 | 17.742 | 4.608 | 3.85× |
| 1 000 000 × 200 000 | 33.632 | 10.022 | 3.36× |
| 2 000 000 × 500 000 | 82.969 | 24.724 | 3.36× |
| 5 000 000 × 1 000 000 | 200.641 | 69.900 | 2.87× |
| 10 000 000 × 2 000 000 | 432.299 | 154.709 | 2.79× |
| 20 000 000 × 4 000 000 | 982.906 | 353.134 | 2.78× |
| 50 000 000 × 10 000 000 | 2 464.018 | 1 324.260 | 1.86× |
| 100 000 000 × 20 000 000 | 6 119.817 | 2 594.778 | 2.36× |
| **200 000 000 × 40 000 000** | **12 979.547** | **4 550.205** | **2.85×** ← MFA flowing through Newton |

**Observations:**

- **Skewed division ratio narrows from ~6.5× at 200k×50k peak to 1.86× at 50M×10M, then rises back to 2.85× at 200M×40M.** Newton inherits the multiplication wins in the large regime, including MFA once the internal products cross the useful MFA band. The residual gap is still division-structure overhead, not decimal I/O.
- 200k×50k stays the worst point: divisor sits below the Newton band (2596 limbs), goes through BZ which loses ~6.5× to GMP's `mpn_dcpi1_div_q` at this size.
- Residual 2.8-3.9× gap in the 1M-20M skewed band is the structural cost of Newton's chunked iteration vs GMP's single recursive divide with precomputed inverse.
- Balanced cases route through FastDivision short-circuits and aren't algorithmically meaningful at this size profile. The 5M×5M balanced case was regressing 27.03× before PR #56 fix (BZ misroute on degenerate quotient); now 7.20× via FastDivision short-circuit.

### Shape-focused division dispatch

`tests/performance/division_shape_bench.cpp` was added after the table above to benchmark meaningful limb-shape division directly. The pass found that the default Base2_64 dispatcher was excluding BZ and sending many near-balanced cases through FastDivision. Dispatch now enables BZ for Base2_64, raises normal Newton entry to 4096-limb divisors, and keeps a high-skew Newton band from 2048-limb divisors at `a >= 8b`.

Representative Base2_64 results:

| shape (limbs) | old dispatch ms | new dispatch ms | main winner |
|---|---:|---:|---|
| `4096/2048` | 5.99 | 3.04 | BZ |
| `6144/4096` | 11.15 | 3.61 | BZ |
| `10240/4096` | 32.73 | 10.24 | BZ |
| `12288/8192` | 44.41 | 7.94 | BZ |
| `20480/2048` (`10n/n`) | 19.15 | 19.19 | Newton |

---

## Decimal parse (string → BigInteger)

| size (digits) | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 1 000 | 0.002 | 0.002 | 1.50× |
| 10 000 | 0.115 | 0.038 | 3.01× |
| 50 000 | 1.375 | 0.387 | 3.43× |
| 100 000 | 3.263 | 1.047 | 3.10× |
| 500 000 | 22.262 | 9.270 | 2.40× |
| 1 000 000 | 49.037 | 21.071 | 2.33× |
| 2 000 000 | 106.786 | 48.047 | 2.22× |
| 5 000 000 | 269.753 | 150.835 | 1.79× |
| 10 000 000 | 585.600 | 350.776 | 1.67× |
| 20 000 000 | 1 270.039 | 814.355 | 1.56× |
| 50 000 000 | 5 168.943 | 2 689.193 | 1.92× |

**Observation:** ratio narrows through the 10M-20M sweet spot (3.2× at 100k → **1.56× at 20M**) where BigMath's NTT overtakes GMP's basecase. It widens back to 1.92× at 50M as GMP's SSA activates.

---

## ToString (BigInteger → string)

| size (digits) | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 1 000 | 0.007 | 0.004 | 1.82× |
| 10 000 | 0.274 | 0.078 | 3.50× |
| 50 000 | 4.044 | 0.861 | 4.70× |
| 100 000 | 19.803 | 2.357 | 8.40× |
| 200 000 | 40.760 | 6.417 | 6.35× |
| 500 000 | 106.486 | 20.667 | 5.15× |
| 1 000 000 | 224.755 | 50.174 | 4.48× |
| 2 000 000 | 481.526 | 119.342 | 4.03× |
| 5 000 000 | 1 104.781 | 386.146 | 2.86× |
| 10 000 000 | 2 416.172 | 913.741 | 2.64× |
| 20 000 000 | 5 432.991 | 2 155.828 | 2.52× |

**Observation:** narrowest gap at 1k (the linear leaf, where GM div2by1 already runs); peaks at 100k (D&C overhead + Newton recip setup not yet amortized); narrows again from 200k onward as D&C asymptotic + NTT inheriting from multiplication's overtake compound — **8.09× at 100k → 2.57× at 20M**.

### ToString focused warm benchmark

After the table above, `tests/performance/tostring_bench.cpp` was added to isolate BigMath conversion cost without GMP and to measure warm-cache optimization deltas. The 2026-05-26 ToString pass added bit-length digit estimation, thread-local cached D&C divider chains, 10¹⁹ linear formatting chunks, digit-pair ASCII emission, and a boundary `NewtonDivision::Divider::DivideAndRemainderInto` API.

| size (digits) | pre-pass ms | post-pass ms | speedup |
|---|---:|---:|---:|
| 1 000 | 0.0079 | 0.0062 | 1.27× |
| 10 000 | 1.1224 | 0.5097 | 2.20× |
| 50 000 | 9.4674 | 4.4633 | 2.12× |
| 100 000 | 19.8435 | 9.4409 | 2.10× |
| 200 000 | 41.1607 | 20.5682 | 2.00× |
| 500 000 | — | 61.1404 | — |
| 1 000 000 | — | 135.8015 | — |
| 2 000 000 | — | 299.5458 | — |

The dominant win is cached divider-chain reuse. The `DivideAndRemainderInto` API is currently neutral because it delegates to existing Newton internals; deeper scratch-buffer reuse would need to be pushed into `DivideChunk`, multiplication, and subtraction temporaries.

---

## BigDecimal

Arbitrary-precision signed decimal. GMP has no native BigDecimal; compared against `mpf_t` at matching precision (`bits = total_digits * 3.322 + 64`).

### Arithmetic

Balanced operands (e.g. `100.10` = 100 integer + 10 fractional digits):

| op | size (int.frac) | BigMath ms | GMP (mpf) ms | BM/GMP |
|---|---|---:|---:|---:|
| Add | 100.10 | 0.000 | 0.000 | 0.00× |
| Add | 1000.100 | 0.000 | 0.000 | 0.00× |
| Add | 5000.500 | 0.001 | 0.000 | 5.66× |
| Add | 20000.2000 | 0.003 | 0.001 | 5.00× |
|---|---|---|---|---|
| Mul | 100.10 | 0.000 | 0.000 | 3.95× |
| Mul | 1000.100 | 0.003 | 0.001 | 2.33× |
| Mul | 5000.500 | 0.038 | 0.018 | 2.12× |
| Mul | 20000.2000 | 0.345 | 0.091 | 3.78× |
|---|---|---|---|---|
| Div | 100.10 (10 dp) | 0.001 | 0.000 | 6.53× |
| Div | 1000.100 (100 dp) | 0.003 | 0.002 | 1.41× |
| Div | **5000.500 (500 dp)** | **0.023** | **0.026** | **0.89×** ← BigMath faster |
| Div | 20000.2000 (2000 dp) | 0.233 | 0.198 | 1.18× |

Division at varying target scales (operand = 2000 integer + 200 fractional digits):

| scale (dp) | BigMath ms | GMP (mpf) ms | BM/GMP |
|---|---:|---:|---:|
| **0** | **0.001** | **0.005** | **0.20×** ← BigMath faster |
| **10** | **0.003** | **0.005** | **0.61×** ← BigMath faster |
| **100** | **0.005** | **0.005** | **0.85×** ← BigMath faster |
| 1000 | 0.015 | 0.010 | 1.57× |
| 5000 | 0.061 | 0.034 | 1.79× |

### Parse and ToString

Parse (`string → BigDecimal`):

| size (digits) | BigMath ms | GMP (mpf) ms | BM/GMP |
|---|---:|---:|---:|
| 100 | 0.001 | 0.001 | 1.00× |
| 1 000 | 0.005 | 0.005 | 1.10× |
| 10 000 | 0.138 | 0.107 | 1.29× |
| 50 000 | 1.504 | 1.017 | 1.48× |

ToString (`BigDecimal → string`):

| size (digits) | BigMath ms | GMP (mpf) ms | BM/GMP |
|---|---:|---:|---:|
| 100 | 0.000 | 0.000 | 0.64× |
| 1 000 | 0.006 | 0.004 | 1.41× |
| 10 000 | 0.269 | 0.102 | 2.63× |
| 50 000 | 4.003 | 1.094 | 3.66× |

**Observations:**
- BigDecimal division beats GMP at 5000.500 to 500 dp (0.023 vs 0.026 ms).
- Scale-aware division avoids unnecessary high-precision work; up to **5× faster** than `mpf_div` at small target scales.
- Parse stays within 1.0-1.5× of GMP across 100-50k digits; ToString stays within 0.7-3.7×.

---

## Headline summary

- **Multiplication 5M-20M balanced:** BigMath is faster at 5M and 10M, then slips back to parity at 20M. The current 10M balanced row is **106.856 ms vs 209.888 ms**, or 0.51×.
- **Multiplication ≥50M balanced:** GMP still wins via SSA. At 200M the gap widens again to 2.16× on this snapshot.
- **Multiplication skewed:** parity around 500k×50k and 1M×100k, with a slight BigMath lead at 2M×200k. BigMath falls behind at 50M×5M, 100M×10M, and 200M×20M.
- **Division skewed:** 1.86×-6.51× behind GMP in the main band, with the new 200M×40M row at 2.85×. Worst remains 200k×50k (BZ band). The large-multiplication speedups still flow through Newton, but the structure overhead never disappears.
- **Division balanced 5M×5M:** PR #56 fix routes degenerate-quotient cases to FastDivision; 27.03× → 7.20×.
- **Parse:** **1.56× at 20M** (best), widens to 1.92× at 50M as GMP's SSA path dominates.
- **ToString:** 2.52× at 20M, narrowing from an 8.40× peak at 100k.
- **BigDecimal division:** beats GMP at small target scales (0.20-1.00×) and at the 500 dp near-parity point (0.89×).

For optimizations considered and rejected with measurement evidence, see the **Explored but rejected** sections of each subsystem doc. The 2026-05 optimization stack (LIMB_64 + CRT NTT + threading + M-G reciprocals + BZ Knuth fix + degenerate-quotient guard + radix-4/radix-8 fused butterflies + MFA) closed the GMP gap by 3-5× across every band up to the 20M crossover; past 20M, GMP's SSA still wins but the loss factor halved.

---

## Dispatch validity after MFA addition

PR #65 added MFA / Bailey 6-step routing inside `NTTMultiplicationCrt`. The default was retuned from `n ≥ 2^21` to `n ≥ 2^24` transform coefficients on 2026-05-27. For Base2_64 balanced multiplication that starts around **2M limbs per operand** (≈40M decimal digits), because each 64-bit limb produces two CRT coefficients and the convolution has about `4L` coefficients. The retune avoids the measured `2^21`-`2^23` regression band while preserving the `2^24+` wins. All other dispatch thresholds are unchanged in semantics. Verified post-MFA via `tests/performance/dispatch_tuner --full` (2026-05-30):

| Knob | Current | Tuner suggestion | Verdict |
|---|---:|---:|---|
| `CLASSIC_MULTIPLICATION_THRESHOLD` (total limbs) | 96 | 96 (= 48 per-operand × 2) | Match |
| `TOOM3_MULTIPLICATION_THRESHOLD` (total limbs) | 2560 | (added 2026-05-27 from focused band scan) | New |
| `NTT_MULTIPLICATION_THRESHOLD` (total limbs) | 5120 | (raised from 4096 by focused scan) | Updated |
| `TOOM3_SKEW_RATIO` | 2 | (added 2026-05-30 from skew-focused scan) | New |
| `CLASSIC_MIN_LIMB_THRESHOLD` | 0 | 0 | Match |
| `NTT_SQUARE_THRESHOLD` (per limb, LIMB_64) | 2048 | 2048 (NTT wins by < 5% but does win) | Keep — current is valid |
| `BZ_DIVISOR_THRESHOLD` | 512 | 768 | Keep — current 512 already correct (BZ activates at b > 512; the 1024/512 tie is FastDiv noise) |
| `NEWTON_MEDIUM_B` | 4096 | 8192 | Keep — current dispatch uses NEWTON_MEDIUM_B together with `a ≥ 3b`, all measured 8192-divisor cases at ratio ≥ 4 still win Newton |
| `NEWTON_SKEW_NUMERATOR/DENOMINATOR` | 3/1 | 4/1 | Keep — bench rows from 5M×1M onward confirm Newton wins at ratio 5 (BigMath path) |
| `NEWTON_HIGH_SKEW_*` | 2048 / 8/1 | (n/a, no rows) | Keep |

**Toom-3 added to dispatch (2026-05-27, skew-gated 2026-05-30):** focused band scan (above, **Toom-3 dispatch band**) found Toom-3 wins by 4-10% in total ∈ [3 584, 4 096] and avoids a 33% NTT-length boundary regression at total 4 608. Dispatch slots Toom-3 between Karatsuba and NTT for total ∈ [2 560, 5 120), but only for non-skewed operands (`max < 2*min`). A skew-focused scan found 2:1+ shapes in that same window are faster through Karatsuba, with wins such as `2560x256` improving from 0.730 ms → 0.343 ms through `Multiply()`. `NTT_MULTIPLICATION_THRESHOLD` remains 5 120. Toom-5 still has no productive band — ties Karatsuba in 64-256 per-operand limbs and degrades sharply above 256, so stays excluded.
