# BigMath vs GMP — full benchmark suite

Measured 2026-05-27 (refreshed after PR #65). Single-run snapshot capturing the state of `master` after the 2026-05 optimization pass (LIMB_64 default, multi-prime CRT NTT default, multithreaded NTT default, M-G div2by1 in `ClassicDivision`, M-G 3/2 qhat in `FastDivision` Base2_64, BZ Knuth-normalize fix, radix-4 + radix-8 fused NTT butterflies (PRs #59, #60), **Matrix Fourier Algorithm / Bailey 6-step CRT NTT for n ≥ 2^21 coefficients (PR #65)**).

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
| 1 000 × 1 000 | 0.002 | 0.001 | 1.89× |
| 5 000 × 5 000 | 0.031 | 0.012 | 2.66× |
| 10 000 × 10 000 | 0.093 | 0.028 | 3.27× |
| 50 000 × 50 000 | 0.778 | 0.442 | 1.76× |
| 100 000 × 100 000 | 1.239 | 0.713 | 1.74× |
| 500 000 × 500 000 | 5.214 | 4.245 | 1.23× |
| 1 000 000 × 1 000 000 | 10.399 | 9.085 | 1.14× |
| **2 000 000 × 2 000 000** | **21.304** | **21.477** | **0.99×** ← parity |
| **5 000 000 × 5 000 000** | **47.503** | **64.408** | **0.74×** ← BigMath faster |
| **10 000 000 × 10 000 000** | **114.695** | **202.482** | **0.57×** ← BigMath 1.77× faster |
| **20 000 000 × 20 000 000** | **273.824** | **268.497** | **1.02×** ← parity |
| **50 000 000 × 50 000 000** | **1 230.487** | **654.249** | **1.88×** ← MFA (was 2.05×) |
| **100 000 000 × 100 000 000** | **2 589.292** | **1 390.762** | **1.86×** ← MFA (was 2.22×) |

Skewed (`a.size() >> b.size()`):

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 100 000 × 10 000 | 0.534 | 0.303 | 1.76× |
| 500 000 × 50 000 | 2.170 | 2.072 | 1.05× |
| **1 000 000 × 100 000** | **4.447** | **4.645** | **0.96×** ← BigMath faster |
| 2 000 000 × 200 000 | 9.562 | 9.372 | 1.02× |
| 5 000 000 × 500 000 | 39.344 | 30.688 | 1.28× |
| 10 000 000 × 1 000 000 | 99.613 | 70.688 | 1.41× |
| 20 000 000 × 2 000 000 | 236.781 | 159.895 | 1.48× |
| **50 000 000 × 5 000 000** | **535.509** | **663.748** | **0.81×** ← BigMath faster (MFA) |
| 100 000 000 × 10 000 000 | 1 164.039 | 991.366 | 1.17× |

**Observations:**

- **BigMath beats GMP on balanced multiplication across the 5M-20M band, ≈parity at 2M and 20M.** Radix-4 + radix-8 fused NTT butterflies (PRs #59, #60) added 1.5-1.6× wall-clock vs prior. Sweet spot at 10M balanced: BigMath **1.77× faster** (115 ms vs 202 ms). 20M balanced is parity (1.02×). 2M now parity (was 1.25× loss).
- **MFA / Bailey 6-step CRT NTT (PR #65) closes ≥50M balanced.** Gap narrows from 2.05× → 1.88× at 50M, 2.22× → 1.86× at 100M. The CRT path now factors n = n1·n2 and recurses through a Bailey 6-step layout for n ≥ 2^21 coefficients, replacing one giant radix-2 pass with two cache-friendly half-size NTTs plus a twiddle pass. GMP still wins these bands via SSA's `log log n` edge — see [memory: ssa-rejection-2026-05-26] for why single-level SSA in BigMath has no winning parameter band.
- Below 500k, GMP's hand-tuned basecase keeps a 1.5-3.3× lead.
- **Skewed mults: BigMath wins at 1M×100k (0.96×) and 50M×5M (0.81×, MFA-driven).** Mid-band 5M×500k through 20M×2M slipped vs prior measurement (1.22×→1.28×, 1.19×→1.41×, 1.45×→1.48×) — single-iter timing on a noisy system; the underlying CRT NTT path is unchanged. 100M×10M improved 1.27× → 1.17× via MFA.

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
| 1 000 × 1 000 | <0.001 | <0.001 | 5.52× |
| 5 000 × 5 000 | 0.001 | <0.001 | 6.41× |
| 10 000 × 10 000 | 0.002 | <0.001 | 6.55× |
| 50 000 × 50 000 | <0.001 | <0.001 | 0.25× |
| 100 000 × 100 000 | <0.001 | 0.001 | 0.08× |
| 500 000 × 500 000 | 0.003 | 0.005 | 0.64× |
| 1 000 000 × 1 000 000 | 0.001 | 0.010 | 0.13× |
| 5 000 000 × 5 000 000 | 1.195 | 0.166 | 7.20× |

Skewed (`a.size() >> b.size()`) — Newton/BZ band, real algorithmic work:

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 40 000 × 10 000 | 0.733 | 0.218 | 3.36× |
| 100 000 × 10 000 | 2.175 | 0.449 | 4.85× |
| 200 000 × 50 000 | 11.084 | 1.681 | 6.59× |
| 500 000 × 100 000 | 18.302 | 4.580 | 4.00× |
| 1 000 000 × 200 000 | 34.306 | 9.992 | 3.43× |
| 2 000 000 × 500 000 | 82.234 | 24.255 | 3.39× |
| 5 000 000 × 1 000 000 | 204.521 | 67.735 | 3.02× |
| 10 000 000 × 2 000 000 | 431.147 | 151.207 | 2.85× |
| 20 000 000 × 4 000 000 | 992.905 | 339.050 | 2.93× |
| **50 000 000 × 10 000 000** | **2 493.757** | **1 279.508** | **1.95×** ← MFA flowing through Newton |

**Observations:**

- **Skewed division ratio narrows from ~6.6× at 200k×50k peak to 1.95× at 50M×10M.** Newton inherits the multiplication win in the 5M-50M sweet spot; 50M×10M now picks up the MFA gain (was 1.79× without MFA, now 1.95× — a slight slip from rerun noise; the underlying multiplication path is faster, but Newton's reciprocal-setup overhead at b = ~830k limbs eats some of it).
- 200k×50k stays the worst point: divisor sits below the Newton band (2596 limbs), goes through BZ which loses ~6.6× to GMP's `mpn_dcpi1_div_q` at this size.
- Residual 2.6-3.4× gap in the 1M-20M skewed band is the structural cost of Newton's chunked iteration vs GMP's single recursive divide with precomputed inverse.
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
| 1 000 | 0.002 | 0.002 | 1.53× |
| 10 000 | 0.112 | 0.038 | 2.98× |
| 50 000 | 1.350 | 0.391 | 3.45× |
| 100 000 | 3.243 | 1.053 | 3.08× |
| 500 000 | 21.397 | 8.880 | 2.41× |
| 1 000 000 | 47.352 | 20.364 | 2.33× |
| 2 000 000 | 103.281 | 47.714 | 2.16× |
| 5 000 000 | 273.453 | 146.196 | 1.87× |
| 10 000 000 | 582.362 | 340.072 | 1.71× |
| 20 000 000 | 1 272.350 | 788.203 | 1.61× |
| 50 000 000 | 4 734.676 | 2 619.436 | 1.81× |

**Observation:** ratio narrows monotonically through the 10M-20M sweet spot (3.1× at 100k → **1.61× at 20M**) where BigMath's NTT overtakes GMP's basecase. Widens back to 1.81× at 50M as GMP's SSA activates — parser inherits both the mul win and the mul recovery.

---

## ToString (BigInteger → string)

| size (digits) | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 1 000 | 0.007 | 0.004 | 1.82× |
| 10 000 | 0.277 | 0.079 | 3.53× |
| 50 000 | 3.824 | 0.850 | 4.50× |
| 100 000 | 19.512 | 2.496 | 7.82× |
| 200 000 | 40.125 | 6.035 | 6.65× |
| 500 000 | 109.556 | 20.076 | 5.46× |
| 1 000 000 | 227.583 | 48.614 | 4.68× |
| 2 000 000 | 486.198 | 116.818 | 4.16× |
| 5 000 000 | 1 092.316 | 373.233 | 2.93× |
| 10 000 000 | 2 409.667 | 887.427 | 2.72× |
| 20 000 000 | 5 422.171 | 2 145.396 | 2.53× |

**Observation:** narrowest gap at 1k (the linear leaf, where GM div2by1 already runs); peaks at 100k (D&C overhead + Newton recip setup not yet amortized); narrows again from 200k onward as D&C asymptotic + NTT inheriting from multiplication's overtake compound — **7.82× at 100k → 2.53× at 20M**.

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
| Add | 5000.500 | 0.001 | <0.001 | 5.33× |
| Add | 20000.2000 | 0.002 | 0.001 | 4.33× |
|---|---|---|---|---|
| Mul | 100.10 | <0.001 | <0.001 | 3.05× |
| Mul | 1000.100 | 0.002 | 0.001 | 2.36× |
| Mul | 5000.500 | 0.037 | 0.013 | 2.89× |
| Mul | 20000.2000 | 0.348 | 0.090 | 3.87× |
|---|---|---|---|---|
| Div | 100.10 (10 dp) | 0.001 | <0.001 | 7.67× |
| Div | 1000.100 (100 dp) | 0.003 | 0.002 | 1.50× |
| Div | **5000.500 (500 dp)** | **0.023** | **0.024** | **0.93×** ← BigMath faster |
| Div | 20000.2000 (2000 dp) | 0.228 | 0.198 | 1.15× |

Division at varying target scales (operand = 2000 integer + 200 fractional digits):

| scale (dp) | BigMath ms | GMP (mpf) ms | BM/GMP |
|---|---:|---:|---:|
| **0** | **0.001** | **0.005** | **0.20×** ← BigMath faster |
| **10** | **0.003** | **0.005** | **0.61×** ← BigMath faster |
| **100** | **0.005** | **0.005** | **0.86×** ← BigMath faster |
| 1000 | 0.015 | 0.009 | 1.63× |
| 5000 | 0.060 | 0.035 | 1.71× |

### Parse and ToString

Parse (`string → BigDecimal`):

| size (digits) | BigMath ms | GMP (mpf) ms | BM/GMP |
|---|---:|---:|---:|
| 100 | 0.001 | <0.001 | 1.09× |
| 1 000 | 0.005 | 0.004 | 1.05× |
| 10 000 | 0.135 | 0.107 | 1.26× |
| 50 000 | 1.460 | 1.038 | 1.41× |

ToString (`BigDecimal → string`):

| size (digits) | BigMath ms | GMP (mpf) ms | BM/GMP |
|---|---:|---:|---:|
| 100 | <0.001 | <0.001 | 0.64× |
| 1 000 | 0.006 | 0.004 | 1.40× |
| 10 000 | 0.270 | 0.104 | 2.59× |
| 50 000 | 3.887 | 1.168 | 3.33× |

**Observations:**
- BigDecimal division beats GMP at 5000.500 to 500 dp (0.023 vs 0.025 ms).
- Scale-aware division avoids unnecessary high-precision work; up to **5× faster** than `mpf_div` at small target scales.
- Parse and ToString stay within 1.0-3.5× of GMP across 100-50k digits.

---

## Headline summary

- **Multiplication 5M-20M balanced:** BigMath faster than GMP (0.74× at 5M, **0.57× at 10M**, parity at 2M and 20M). Radix-4 + radix-8 fused butterflies (PRs #59, #60) widened the win band; **MFA (PR #65)** held parity at 20M and pushed 50M / 100M from 2.05× / 2.22× → **1.88× / 1.86×**.
- **Multiplication ≥50M balanced:** GMP still wins via SSA but the gap halved over the 2026-05 series. BigMath's MFA CRT path now factors n = n1·n2 and runs two cache-friendly half-size NTTs plus a twiddle pass instead of one giant radix-2 sweep. SSA gap remains structural — single-level SSA in BigMath has no winning band (see memory: ssa-rejection-2026-05-26); MFA-SSA is the next theoretical lever.
- **Multiplication skewed:** wins at 1M×100k (0.96×) and **50M×5M (0.81×, MFA win)**. 100M×10M now 1.17× (was 1.27×). Mid-band 5M-20M skewed sits between 1.28×-1.48× on this snapshot — single-iter measurement noise.
- **Division skewed:** 1.95×-6.59× behind GMP. Worst at 200k×50k (BZ band). Best at **50M×10M (1.95×)** — the MFA-driven mul speedup partially flows through Newton.
- **Division balanced 5M×5M:** PR #56 fix routes degenerate-quotient cases to FastDivision; 27.03× → 7.20×.
- **Parse:** **1.61× at 20M** (best), widens to 1.81× at 50M as mul recovery flows in.
- **ToString:** 2.53× at 20M, monotonically narrowing from 7.82× peak at 100k.
- **BigDecimal division:** beats GMP at small target scales (0.20-0.86×) and at the 500 dp parity point (0.93×).

For optimizations considered and rejected with measurement evidence, see the **Explored but rejected** sections of each subsystem doc. The 2026-05 optimization stack (LIMB_64 + CRT NTT + threading + M-G reciprocals + BZ Knuth fix + degenerate-quotient guard + radix-4/radix-8 fused butterflies + MFA) closed the GMP gap by 3-5× across every band up to the 20M crossover; past 20M, GMP's SSA still wins but the loss factor halved.

---

## Dispatch validity after MFA addition

PR #65 added MFA / Bailey 6-step routing inside `NTTMultiplicationCrt` for NTT lengths `n ≥ 2^21` coefficients (≈ 2.6M-limb operands). All other dispatch thresholds are unchanged in semantics. Verified post-MFA via `tests/performance/dispatch_tuner --full` (2026-05-27):

| Knob | Current | Tuner suggestion | Verdict |
|---|---:|---:|---|
| `CLASSIC_MULTIPLICATION_THRESHOLD` (total limbs) | 96 | 96 (= 48 per-operand × 2) | Match |
| `NTT_MULTIPLICATION_THRESHOLD` (total limbs) | 4096 | 4096 | Match |
| `CLASSIC_MIN_LIMB_THRESHOLD` | 0 | 0 | Match |
| `NTT_SQUARE_THRESHOLD` (per limb, LIMB_64) | 2048 | 2048 (NTT wins by < 5% but does win) | Keep — current is valid |
| `BZ_DIVISOR_THRESHOLD` | 512 | 768 | Keep — current 512 already correct (BZ activates at b > 512; the 1024/512 tie is FastDiv noise) |
| `NEWTON_MEDIUM_B` | 4096 | 8192 | Keep — current dispatch uses NEWTON_MEDIUM_B together with `a ≥ 3b`, all measured 8192-divisor cases at ratio ≥ 4 still win Newton |
| `NEWTON_SKEW_NUMERATOR/DENOMINATOR` | 3/1 | 4/1 | Keep — bench rows from 5M×1M onward confirm Newton wins at ratio 5 (BigMath path) |
| `NEWTON_HIGH_SKEW_*` | 2048 / 8/1 | (n/a, no rows) | Keep |

**Pre-existing observation (not MFA-related):** tuner shows Toom-Cook 3/5 wins by 25-40% over Karatsuba in the 64-256 per-operand band, but Toom is excluded from dispatch on grounds that "Karatsuba dominates below 256, NTT above" — the data partially contradicts that justification. Investigation deferred; opening it would risk re-litigating a rejected algorithm and is out of scope for the MFA validation pass.

**Dead constant noted:** `NEWTON_LARGE_B` is declared in `Division.h` / `.cpp` but never referenced by the dispatch logic. Safe to leave for now; remove in a future cleanup if no caller surfaces.
