# BigMath vs GMP — full benchmark suite

Measured 2026-05-26. Single-run snapshot capturing the state of `master` after the 2026-05 optimization pass (LIMB_64 default, multi-prime CRT NTT default, multithreaded NTT default, M-G div2by1 in `ClassicDivision`, M-G 3/2 qhat in `FastDivision` Base2_64, BZ Knuth-normalize fix).

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
| 1 000 × 1 000 | 0.002 | 0.001 | 1.82× |
| 5 000 × 5 000 | 0.031 | 0.012 | 2.69× |
| 10 000 × 10 000 | 0.093 | 0.028 | 3.30× |
| 50 000 × 50 000 | 0.776 | 0.441 | 1.76× |
| 100 000 × 100 000 | 1.189 | 0.725 | 1.64× |
| 500 000 × 500 000 | 6.781 | 4.248 | 1.60× |
| 1 000 000 × 1 000 000 | 12.803 | 8.961 | 1.43× |
| 2 000 000 × 2 000 000 | 26.322 | 21.027 | 1.25× |
| **5 000 000 × 5 000 000** | **58.266** | **63.340** | **0.92×** ← BigMath faster |
| **10 000 000 × 10 000 000** | **144.480** | **204.184** | **0.71×** ← BigMath faster |
| 20 000 000 × 20 000 000 | 432.353 | 268.309 | 1.61× ← GMP recovers |
| 50 000 000 × 50 000 000 | 2 423.370 | 676.038 | 3.58× |
| 100 000 000 × 100 000 000 | 5 313.011 | 1 433.202 | 3.71× |

Skewed (`a.size() >> b.size()`):

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 100 000 × 10 000 | 0.506 | 0.302 | 1.67× |
| 500 000 × 50 000 | 2.431 | 2.089 | 1.16× |
| 1 000 000 × 100 000 | 4.806 | 4.557 | 1.05× |
| **2 000 000 × 200 000** | **9.536** | **9.381** | **1.02×** ← parity |
| 5 000 000 × 500 000 | 45.262 | 30.901 | 1.46× |
| 10 000 000 × 1 000 000 | 111.951 | 71.620 | 1.56× |
| 20 000 000 × 2 000 000 | 389.544 | 159.945 | 2.44× |
| 50 000 000 × 5 000 000 | 940.123 | 678.768 | 1.39× |
| 100 000 000 × 10 000 000 | 2 224.512 | 1 000.753 | 2.22× |

**Observations:**

- **BigMath beats GMP on balanced multiplication in the 5M-10M band.** NTT's asymptotic edge overtakes GMP's hand-tuned Karatsuba/Toom assembly. At 10M balanced, BigMath is **29% faster** (144ms vs 204ms).
- **GMP recovers at ≥20M balanced.** Crossover reverses: 1.61× at 20M, 3.58× at 50M, 3.71× at 100M. GMP's Schönhage-Strassen activates here, and its FFT scales sub-linearly above the BigMath NTT band. BigMath's CRT NTT is bound by L2/L3 cache footprint (~250 MB working set at 100M operands) and the basecase butterfly inner loop. The 5M-10M sweet spot is where BigMath's CRT-parallelized NTT beats GMP's pre-SSA path.
- Below 1M, GMP's hand-tuned basecase keeps a 1.4-3.3× lead.
- Skewed mults stay within a 1.0-1.7× band up to 10M×1M; widens to 2.2-2.4× at 50M-100M / 10:1 ratio as the dividend hits the same post-10M NTT-loses-to-SSA regime.

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
| 1 000 × 1 000 | <0.001 | <0.001 | 5.02× |
| 5 000 × 5 000 | 0.001 | <0.001 | 8.03× |
| 10 000 × 10 000 | 0.002 | <0.001 | 6.56× |
| 50 000 × 50 000 | <0.001 | <0.001 | 0.28× |
| 100 000 × 100 000 | <0.001 | 0.001 | 0.17× |
| 500 000 × 500 000 | 0.003 | 0.005 | 0.64× |
| 1 000 000 × 1 000 000 | 0.004 | 0.010 | 0.39× |
| 5 000 000 × 5 000 000 | 1.220 | 0.168 | 7.25× |

Skewed (`a.size() >> b.size()`) — Newton/BZ band, real algorithmic work:

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 40 000 × 10 000 | 0.745 | 0.220 | 3.39× |
| 100 000 × 10 000 | 2.181 | 0.449 | 4.86× |
| 200 000 × 50 000 | 11.200 | 1.699 | 6.59× |
| 500 000 × 100 000 | 18.546 | 4.533 | 4.09× |
| 1 000 000 × 200 000 | 39.216 | 9.955 | 3.94× |
| 2 000 000 × 500 000 | 88.915 | 24.511 | 3.63× |
| 5 000 000 × 1 000 000 | 214.274 | 69.274 | 3.09× |
| 10 000 000 × 2 000 000 | 463.840 | 156.530 | 2.96× |
| 20 000 000 × 4 000 000 | 1 090.923 | 349.810 | 3.12× |
| 50 000 000 × 10 000 000 | 3 303.338 | 1 322.413 | 2.50× |

**Observations:**

- **Skewed division ratio narrows from ~6.6× at 200k×50k peak to ~2.5-3× at 10M-50M.** Newton inherits the multiplication win in the 5M-10M sweet spot. At 50M×10M the ratio is 2.50× — Newton's per-chunk multiplications fall in the BigMath-wins band.
- 200k×50k stays the worst point: divisor sits below the Newton band (2596 limbs), goes through BZ which loses ~6.6× to GMP's `mpn_dcpi1_div_q` at this size.
- Residual ~2.5-4× gap is the structural cost of Newton's chunked iteration vs GMP's single recursive divide with precomputed inverse.
- Balanced cases route through FastDivision short-circuits and aren't algorithmically meaningful at this size profile. The 5M×5M balanced case was regressing 27.03× before PR #56 fix (BZ misroute on degenerate quotient); now 7.25× via FastDivision short-circuit.

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
| 1 000 | 0.003 | 0.002 | 1.63× |
| 10 000 | 0.115 | 0.038 | 3.00× |
| 50 000 | 1.350 | 0.390 | 3.46× |
| 100 000 | 3.230 | 1.058 | 3.05× |
| 500 000 | 21.879 | 8.864 | 2.47× |
| 1 000 000 | 48.063 | 20.745 | 2.32× |
| 2 000 000 | 106.185 | 48.690 | 2.18× |
| 5 000 000 | 277.930 | 149.275 | 1.86× |
| 10 000 000 | 599.505 | 348.250 | 1.72× |
| 20 000 000 | 1 317.623 | 804.087 | 1.64× |
| 50 000 000 | 5 571.066 | 2 680.756 | 2.08× |

**Observation:** ratio narrows monotonically through the 5M-10M sweet spot (3.0× at 100k → 1.64× at 20M) where BigMath's NTT overtakes GMP's basecase. Widens back to 2.08× at 50M as GMP's SSA activates — parser inherits both the mul win and the mul recovery.

---

## ToString (BigInteger → string)

| size (digits) | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 1 000 | 0.006 | 0.003 | 1.82× |
| 10 000 | 0.274 | 0.079 | 3.50× |
| 50 000 | 3.874 | 0.854 | 4.54× |
| 100 000 | 19.637 | 2.387 | 8.23× |
| 200 000 | 40.453 | 6.152 | 6.58× |
| 500 000 | 115.198 | 20.407 | 5.64× |
| 1 000 000 | 242.821 | 49.591 | 4.90× |
| 2 000 000 | 513.827 | 119.879 | 4.29× |
| 5 000 000 | 1 161.977 | 381.892 | 3.04× |
| 10 000 000 | 2 564.660 | 904.627 | 2.84× |
| 20 000 000 | 5 951.212 | 2 129.169 | 2.80× |

**Observation:** narrowest gap at 1k (the linear leaf, where GM div2by1 already runs); peaks at 100k (D&C overhead + Newton recip setup not yet amortized); narrows again from 200k onward as D&C asymptotic + NTT inheriting from multiplication's overtake compound — 8.23× at 100k → 2.80× at 20M.

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

- **Multiplication 5M-10M balanced:** BigMath faster than GMP (0.92× at 5M, 0.71× at 10M).
- **Multiplication ≥20M balanced:** GMP recovers as SSA activates (1.61× at 20M → 3.71× at 100M). BigMath's CRT NTT is L2/L3-cache-bound here; SSA's `log log n` edge pays off past 10M operand sizes.
- **Multiplication skewed:** parity 1.0× at 2M×200k; 1.5-1.7× at 5M-10M skewed; 2.2-2.4× at 50M-100M (10:1 ratio).
- **Division skewed:** 2.5-6.6× behind GMP. Worst at 200k×50k (BZ band). Best at 50M×10M (2.5×) as multiplication win flows through Newton.
- **Division balanced 5M×5M:** PR #56 fix routes degenerate-quotient cases to FastDivision; 27.03× → 7.25×.
- **Parse:** 1.64× at 20M (best), widens to 2.08× at 50M as mul recovery flows in.
- **ToString:** 2.80× at 20M, monotonically narrowing from 8.23× peak at 100k.
- **BigDecimal division:** beats GMP at small target scales (0.20-0.85×) and at the 500 dp parity point (0.92×).

For optimizations considered and rejected with measurement evidence, see the **Explored but rejected** sections of each subsystem doc. The 2026-05 optimization stack (LIMB_64 + CRT NTT + threading + M-G reciprocals + BZ Knuth fix + degenerate-quotient guard) closed the GMP gap by 3-5× across every band up to the 10M crossover; past 10M, GMP's SSA reverses the lead.
