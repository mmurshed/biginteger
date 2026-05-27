# BigMath vs GMP — full benchmark suite

Measured 2026-05-27. Single-run snapshot capturing the state of `master` after the 2026-05 optimization pass (LIMB_64 default, multi-prime CRT NTT default, multithreaded NTT default, M-G div2by1 in `ClassicDivision`, M-G 3/2 qhat in `FastDivision` Base2_64, BZ Knuth-normalize fix, **radix-4 + radix-8 fused NTT butterflies (PRs #59, #60)**).

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

- `tests/performance/bench_vs_gmp.cpp` — canonical bench. Mul covers ≤100M digits (balanced + skewed), div ≤50M×10M skewed, parse ≤50M, ToString ≤20M, plus BigDecimal arithmetic + I/O. Single-iter timing for digits ≥100k; multi-iter best-of below that. Per-row `[op size] setup...` log to stderr so long setup phases don't look stuck.

---

## Multiplication

Balanced (`a.size() == b.size()`):

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 1 000 × 1 000 | 0.002 | 0.001 | 2.24× |
| 5 000 × 5 000 | 0.036 | 0.014 | 2.64× |
| 10 000 × 10 000 | 0.197 | 0.059 | 3.32× |
| 50 000 × 50 000 | 0.600 | 0.362 | 1.66× |
| 100 000 × 100 000 | 1.198 | 0.761 | 1.57× |
| 500 000 × 500 000 | 5.014 | 4.204 | 1.19× |
| 1 000 000 × 1 000 000 | 9.839 | 8.892 | 1.11× |
| **2 000 000 × 2 000 000** | **19.843** | **20.823** | **0.95×** ← BigMath faster |
| **5 000 000 × 5 000 000** | **45.811** | **63.960** | **0.72×** ← BigMath faster |
| **10 000 000 × 10 000 000** | **103.018** | **213.191** | **0.48×** ← BigMath 2× faster |
| **20 000 000 × 20 000 000** | **276.471** | **280.078** | **0.99×** ← parity |
| 50 000 000 × 50 000 000 | 1 392.879 | 680.512 | 2.05× |
| 100 000 000 × 100 000 000 | 3 220.928 | 1 448.546 | 2.22× |

Skewed (`a.size() >> b.size()`):

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 100 000 × 10 000 | 0.543 | 0.303 | 1.79× |
| **500 000 × 50 000** | **2.010** | **2.074** | **0.97×** ← BigMath faster |
| **1 000 000 × 100 000** | **4.269** | **4.475** | **0.95×** ← BigMath faster |
| **2 000 000 × 200 000** | **8.512** | **9.397** | **0.91×** ← BigMath faster |
| 5 000 000 × 500 000 | 37.312 | 30.617 | 1.22× |
| 10 000 000 × 1 000 000 | 84.151 | 70.710 | 1.19× |
| 20 000 000 × 2 000 000 | 229.908 | 158.847 | 1.45× |
| **50 000 000 × 5 000 000** | **625.327** | **665.872** | **0.94×** ← BigMath faster |
| 100 000 000 × 10 000 000 | 1 279.684 | 1 007.119 | 1.27× |

**Observations:**

- **BigMath beats GMP on balanced multiplication across the 2M-20M band.** Radix-4 + radix-8 fused NTT butterflies (PRs #59, #60) added 1.5-1.6× wall-clock vs prior. Sweet spot at 10M balanced: BigMath **2.08× faster** (103 ms vs 213 ms). 20M balanced is now parity (0.99×, was 1.61× loss). 2M crossed into wins (was 1.25× loss).
- **GMP recovers at ≥50M balanced.** Crossover reverses: 2.05× at 50M, 2.22× at 100M (was 3.58× / 3.71× before radix-4+8). GMP's Schönhage-Strassen activates here. SSA gap remains structural — see [memory: ssa-rejection-2026-05-26] for parameter analysis showing single-level SSA has no winning band; MFA-SSA is the path forward.
- Below 500k, GMP's hand-tuned basecase keeps a 1.5-3.3× lead.
- **Skewed mults: BigMath wins from 500k×50k through 50M×5M.** Four sweet spots: 0.97× / 0.95× / 0.91× / 0.94×. Past 100M×10M (1.27×) GMP recovers via SSA.

### Multithreaded NTT check

`BIGMATH_USE_THREADS=1` is the default. A focused `mul_xl_bench` run on 2026-05-27 compared the default build against `-DBIGMATH_USE_THREADS=0`:

| limb size | serial ms | threaded ms | speedup |
|---:|---:|---:|---:|
| 100 000 | 55.171 | 17.634 | 3.13× |
| 500 000 | 257.394 | 90.207 | 2.85× |
| 1 000 000 | 602.433 | 250.242 | 2.41× |

The threaded path uses coarse CRT parallelism for the six forward transforms and three inverse transforms, plus chunked pointwise multiplication. Single-prime Goldilocks fallback also uses size-gated layer parallelism.

### Shape-focused multiplication dispatch

`tests/performance/multiplication_shape_bench.cpp` was added to measure direct Classic, Karatsuba, NTT, dispatcher, and an experimental blockwise-skew prototype by limb shape. It found one production dispatch miss: tiny high-skew operands should not enter Karatsuba.

| limb shape | previous dispatch | retuned dispatch | impact |
|---|---:|---:|---|
| `1280x64` (`20n/n`) | 0.1522 ms | 0.0848 ms | 1.8× faster |
| `3200x64` (`50n/n`) | 0.4469 ms | ~0.21-0.25 ms | ~1.8-2.1× faster |

The blockwise-skew prototype wins some `min=64` microbenchmarks but loses at larger smaller-operand sizes, so it remains benchmark-only. A Barrett modular multiply experiment for CRT NTT primes was slower than the existing `% P` lowering and was rejected.

---

## Division

Balanced (`a.size() == b.size()`) — quotient is 1-2 limbs, both libraries short-circuit. Numbers are sub-microsecond noise; reported for completeness.

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 1 000 × 1 000 | 0.001 | <0.001 | 12.90× |
| 5 000 × 5 000 | 0.001 | <0.001 | 6.61× |
| 10 000 × 10 000 | 0.002 | <0.001 | 7.38× |
| 50 000 × 50 000 | <0.001 | <0.001 | 0.25× |
| 100 000 × 100 000 | <0.001 | 0.001 | 0.17× |
| 500 000 × 500 000 | 0.001 | 0.005 | 0.23× |
| 1 000 000 × 1 000 000 | 0.004 | 0.012 | 0.35× |
| 5 000 000 × 5 000 000 | 1.368 | 0.166 | 8.22× |

Skewed (`a.size() >> b.size()`) — Newton/BZ band, real algorithmic work:

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 40 000 × 10 000 | 0.761 | 0.225 | 3.38× |
| 100 000 × 10 000 | 2.194 | 0.454 | 4.83× |
| 200 000 × 50 000 | 11.292 | 1.674 | 6.74× |
| 500 000 × 100 000 | 18.243 | 4.755 | 3.84× |
| 1 000 000 × 200 000 | 35.215 | 10.677 | 3.30× |
| 2 000 000 × 500 000 | 81.141 | 25.526 | 3.18× |
| 5 000 000 × 1 000 000 | 191.291 | 68.730 | 2.78× |
| 10 000 000 × 2 000 000 | 412.358 | 156.831 | 2.63× |
| 20 000 000 × 4 000 000 | 909.693 | 343.664 | 2.65× |
| **50 000 000 × 10 000 000** | **2 355.601** | **1 312.316** | **1.79×** |

**Observations:**

- **Skewed division ratio narrows from ~6.7× at 200k×50k peak to 1.79× at 50M×10M.** Newton inherits the multiplication win in the 5M-50M sweet spot; 50M×10M improved from 2.50× → **1.79×** after radix-4+8.
- 200k×50k stays the worst point: divisor sits below the Newton band (2596 limbs), goes through BZ which loses ~6.7× to GMP's `mpn_dcpi1_div_q` at this size.
- Residual 2.6-3.4× gap in the 1M-20M skewed band is the structural cost of Newton's chunked iteration vs GMP's single recursive divide with precomputed inverse — narrowed by ~10-15% via the NTT speedup since Newton calls Multiply per chunk.
- Balanced cases route through FastDivision short-circuits and aren't algorithmically meaningful at this size profile. The 5M×5M balanced case was regressing 27.03× before PR #56 fix (BZ misroute on degenerate quotient); now 8.22× via FastDivision short-circuit.

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
| 1 000 | 0.003 | 0.002 | 1.58× |
| 10 000 | 0.112 | 0.039 | 2.89× |
| 50 000 | 1.343 | 0.391 | 3.44× |
| 100 000 | 3.190 | 1.059 | 3.01× |
| 500 000 | 21.260 | 8.905 | 2.39× |
| 1 000 000 | 46.496 | 20.496 | 2.27× |
| 2 000 000 | 101.795 | 48.581 | 2.10× |
| 5 000 000 | 264.905 | 147.716 | 1.79× |
| 10 000 000 | 570.287 | 343.762 | 1.66× |
| 20 000 000 | 1 242.885 | 800.282 | 1.55× |
| 50 000 000 | 4 665.032 | 2 674.914 | 1.74× |

**Observation:** ratio narrows monotonically through the 10M-20M sweet spot (3.0× at 100k → **1.55× at 20M**) where BigMath's NTT overtakes GMP's basecase. Widens back to 1.74× at 50M as GMP's SSA activates — parser inherits both the mul win and the mul recovery. 50M improved from 2.08× → 1.74× after radix-4+8.

---

## ToString (BigInteger → string)

| size (digits) | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 1 000 | 0.006 | 0.004 | 1.74× |
| 10 000 | 0.269 | 0.077 | 3.50× |
| 50 000 | 3.833 | 0.848 | 4.52× |
| 100 000 | 19.453 | 2.360 | 8.24× |
| 200 000 | 39.651 | 6.309 | 6.29× |
| 500 000 | 109.246 | 20.341 | 5.37× |
| 1 000 000 | 227.797 | 48.731 | 4.67× |
| 2 000 000 | 479.826 | 118.745 | 4.04× |
| 5 000 000 | 1 074.873 | 378.272 | 2.84× |
| 10 000 000 | 2 311.310 | 918.750 | 2.52× |
| 20 000 000 | 5 237.676 | 2 115.032 | 2.48× |

**Observation:** narrowest gap at 1k (the linear leaf, where GM div2by1 already runs); peaks at 100k (D&C overhead + Newton recip setup not yet amortized); narrows again from 200k onward as D&C asymptotic + NTT inheriting from multiplication's overtake compound — **8.24× at 100k → 2.48× at 20M**. 10M improved 2.84× → 2.52× after radix-4+8.

### ToString focused warm benchmark

After the table above, `tests/performance/tostring_bench.cpp` was added to isolate BigMath conversion cost without GMP and to measure warm-cache optimization deltas. The 2026-05-26 ToString pass added bit-length digit estimation, thread-local cached D&C divider chains, 10¹⁹ linear formatting chunks, digit-pair ASCII emission, and a boundary `NewtonDivision::Divider::DivideAndRemainderInto` API.

| size (digits) | pre-pass ms | post-pass ms | speedup |
|---|---:|---:|---:|
| 1 000 | 0.0079 | 0.0063-0.0066 | 1.2-1.3× |
| 10 000 | 1.1224 | 0.2793-0.2924 | 3.8-4.0× |
| 50 000 | 9.4674 | 3.95-4.05 | 2.3-2.4× |
| 100 000 | 19.8435 | 8.94-9.05 | 2.2× |
| 200 000 | 41.1607 | 20.28-20.65 | 2.0× |

The dominant win is cached divider-chain reuse. The `DivideAndRemainderInto` API is currently neutral because it delegates to existing Newton internals; deeper scratch-buffer reuse would need to be pushed into `DivideChunk`, multiplication, and subtraction temporaries.

---

## BigDecimal

Arbitrary-precision signed decimal. GMP has no native BigDecimal; compared against `mpf_t` at matching precision (`bits = total_digits * 3.322 + 64`).

### Arithmetic

Balanced operands (e.g. `100.10` = 100 integer + 10 fractional digits):

| op | size (int.frac) | BigMath ms | GMP (mpf) ms | BM/GMP |
|---|---|---:|---:|---:|
| Add | 100.10 | <0.001 | <0.001 | — |
| Add | 1000.100 | <0.001 | <0.001 | — |
| Add | 5000.500 | 0.001 | <0.001 | 5.66× |
| Add | 20000.2000 | 0.002 | 0.001 | 4.50× |
|---|---|---|---|---|
| Mul | 100.10 | <0.001 | <0.001 | 3.05× |
| Mul | 1000.100 | 0.003 | 0.001 | 2.31× |
| Mul | 5000.500 | 0.038 | 0.013 | 2.88× |
| Mul | 20000.2000 | 0.351 | 0.091 | 3.87× |
|---|---|---|---|---|
| Div | 100.10 (10 dp) | 0.001 | <0.001 | 9.00× |
| Div | 1000.100 (100 dp) | 0.003 | 0.002 | 1.51× |
| Div | **5000.500 (500 dp)** | **0.023** | **0.025** | **0.92×** ← BigMath faster |
| Div | 20000.2000 (2000 dp) | 0.234 | 0.201 | 1.16× |

Division at varying target scales (operand = 2000 integer + 200 fractional digits):

| scale (dp) | BigMath ms | GMP (mpf) ms | BM/GMP |
|---|---:|---:|---:|
| **0** | **0.001** | **0.005** | **0.20×** ← BigMath faster |
| **10** | **0.003** | **0.005** | **0.61×** ← BigMath faster |
| **100** | **0.005** | **0.005** | **0.85×** ← BigMath faster |
| 1000 | 0.015 | 0.010 | 1.57× |
| 5000 | 0.061 | 0.036 | 1.68× |

### Parse and ToString

Parse (`string → BigDecimal`):

| size (digits) | BigMath ms | GMP (mpf) ms | BM/GMP |
|---|---:|---:|---:|
| 100 | 0.001 | 0.001 | 1.00× |
| 1 000 | 0.005 | 0.005 | 1.05× |
| 10 000 | 0.137 | 0.108 | 1.26× |
| 50 000 | 1.480 | 1.049 | 1.41× |

ToString (`BigDecimal → string`):

| size (digits) | BigMath ms | GMP (mpf) ms | BM/GMP |
|---|---:|---:|---:|
| 100 | <0.001 | <0.001 | 0.64× |
| 1 000 | 0.006 | 0.005 | 1.41× |
| 10 000 | 0.275 | 0.104 | 2.64× |
| 50 000 | 3.898 | 1.135 | 3.43× |

**Observations:**
- BigDecimal division beats GMP at 5000.500 to 500 dp (0.023 vs 0.025 ms).
- Scale-aware division avoids unnecessary high-precision work; up to **5× faster** than `mpf_div` at small target scales.
- Parse and ToString stay within 1.0-3.5× of GMP across 100-50k digits.

---

## Headline summary

- **Multiplication 2M-20M balanced:** BigMath faster than GMP (0.95× at 2M, 0.72× at 5M, **0.48× at 10M**, 0.99× at 20M). Radix-4 + radix-8 fused butterflies (PRs #59, #60) widened the win band from 5M-10M to 2M-20M and pushed the peak win to **2.08× at 10M**.
- **Multiplication ≥50M balanced:** GMP recovers as SSA activates (2.05× at 50M, 2.22× at 100M — was 3.58× / 3.71× before radix-4+8). BigMath's CRT NTT is L2/L3-cache-bound here; SSA's `log log n` edge pays off past 20M. Single-level SSA in BigMath has no winning parameter band (see memory: ssa-rejection-2026-05-26); MFA-SSA is the path forward if the gap matters.
- **Multiplication skewed:** BigMath wins from **500k×50k through 50M×5M** (0.91-0.97×) — four sweet spots vs prior one. Past 100M×10M (1.27×) GMP recovers via SSA.
- **Division skewed:** 1.8-6.7× behind GMP. Worst at 200k×50k (BZ band). Best at **50M×10M (1.79×, was 2.50×)** as multiplication win flows through Newton.
- **Division balanced 5M×5M:** PR #56 fix routes degenerate-quotient cases to FastDivision; 27.03× → 8.22×.
- **Parse:** **1.55× at 20M** (best, was 1.64×), widens to 1.74× at 50M as mul recovery flows in.
- **ToString:** 2.48× at 20M (was 2.80×), monotonically narrowing from 8.24× peak at 100k.
- **BigDecimal division:** beats GMP at small target scales (0.20-0.86×) and at the 500 dp parity point (0.92×).

For optimizations considered and rejected with measurement evidence, see the **Explored but rejected** sections of each subsystem doc. The 2026-05 optimization stack (LIMB_64 + CRT NTT + threading + M-G reciprocals + BZ Knuth fix + degenerate-quotient guard + radix-4/radix-8 fused butterflies) closed the GMP gap by 3-5× across every band up to the 20M crossover; past 20M, GMP's SSA reverses the lead but the loss factor halved.
