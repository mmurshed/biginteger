# BigMath vs GMP — full benchmark suite

Canonical run **2026-11-21** for division, decimal parse, ToString and
BigDecimal division (post #141–#149: the division pass — raw-buffer basecase,
in-place Burnikel–Ziegler, Newton reciprocal rework, dispatch re-sweep,
quotient-only division, lock-free pool pickup). The multiplication tables are
from the **2026-11-20** run (post #129–#139, the multiplication pass) and
were not re-measured. Single `bench_vs_gmp` pass, with warm steady-state
re-measurement for the chain-/plan-heavy rows where single-iteration cold
timing is not representative (see Methodology). Tables carry the previous
canonical ratios for comparison.

![BigMath vs GMP ratio by operand size](docs/images/bigmath_vs_gmp.png)

*Figure uses the canonical tables below: multiplication 2026-11-20;
division and decimal I/O 2026-11-21. Multiplication at ≥50M digits uses
warm timings; ToString cold and measured warm timings are shown separately.*

**Headlines:**

- **Balanced multiplication beats GMP from 10k digits up** (was 500k) —
  0.21–0.63× from 50k to 20M digits, peaking at 4.8× faster (10M×10M digits:
  43 vs 207 ms), and 0.63–0.72× warm through 200M digits.
- **Skewed multiplication wins every measured shape** (0.27–0.51×), including
  100k×10k, which was 1.21×.
- **Skewed division beats GMP from 200k-digit dividends** (was 500k; 20M at
  v13.0): 0.19–0.67×, e.g. 50M ÷ 10M digits in 244 ms against GMP's 1.31 s.
  The 12 skewed division rows are 1.14–1.58× faster than the 2026-11-20 run.
- **Decimal I/O**: parse beats GMP from 100k digits (0.23–0.54×); ToString
  from 200k digits cold (0.38–0.74×; was 1M) and from 100k warm
  (0.30–0.57×).
- Below ~10k digits multiplication is 1.4–1.5× behind GMP instead of
  2.3–3.2×. The earlier "portable-C++ wall" was a loop idiom, not a language
  limit — see [docs/MULTIPLICATION.md](docs/MULTIPLICATION.md). What still
  trails: division below ~200k-digit dividends (1.4–2.0×, was 2.0–2.9×),
  small decimal I/O (1.4–2.3×), and 1–4-limb products where the result
  allocation dominates.

For per-subsystem deep dives — algorithms, dispatch, optimization history,
rejected approaches — see [docs/MULTIPLICATION.md](docs/MULTIPLICATION.md),
[docs/DIVISION.md](docs/DIVISION.md),
[docs/STRING_CONVERSION.md](docs/STRING_CONVERSION.md),
[docs/BASE.md](docs/BASE.md). Release history: [CHANGELOG.md](CHANGELOG.md).
This file is the raw numbers, all subsystems in one place.

---

The follow-up review records matched Q/Q and QR/QR timings at 100k ÷ 10k digits,
including the time spent in each BZ kernel. [Review and measurements](docs/PERFORMANCE_REVIEW_PR129.md#matched-division-and-the-remaining-gap)

## Methodology

- **Hardware:** Apple M1 Max (10-core, 8 performance + 2 efficiency), 64 GB RAM.
- **OS:** macOS (Darwin 27.2.0).
- **Build:** CMake Release (`-O3 -march=native`), Apple Clang 21, C++20.
- **Reference:** GMP 6.3.0 (Homebrew), timed in-process on the same operands.
  GMP is single-threaded; BigMath uses its 8-thread pool on NTT-sized work.
- **Default stack:** `BIGMATH_LIMB_64=1` (64-bit limbs), `BIGMATH_NTT_CRT=1`
  (3-prime CRT NTT), `BIGMATH_USE_THREADS=1` (8-thread pool),
  `BIGMATH_NEON=1` on aarch64, MFA gate `2^20` coefficients.
- **Timing:** the integer harness takes the minimum of 7 iterations through
  2k digits, 5 through 20k, 3 through 100k and 1 above 100k. ToString uses
  one iteration from 100k digits. Current output includes the repetition count.
  The machine was not idle for this run (load average 3–5 from desktop
  processes), so treat single-iteration rows as ±10%.
- **Warm steady state:** rows marked *warm* discard the first call (NTT plan
  build, divider-chain build, first-touch page faults) and take the best of
  the following runs. This is the authoritative comparison for ≥50M-digit
  multiplication and ≥100k-digit ToString: BigMath amortizes per-process
  caches that GMP doesn't have, so cold-vs-cold pays them once per process
  and never again.
- **Random inputs:** mt19937_64 seeded per case, full digit/limb range.

---

## Multiplication

Balanced (`a.size() == b.size()`):

| digits | BigMath ms | GMP ms | BM/GMP | v13.0 BM/GMP |
|---|---:|---:|---:|---:|
| 1 000 × 1 000 | 0.001 | 0.001 | 1.39× | 2.33× |
| 5 000 × 5 000 | 0.018 | 0.012 | 1.48× | 2.56× |
| **10 000 × 10 000** | **0.027** | **0.030** | **0.90×** | **3.19×** |
| **50 000 × 50 000** | **0.182** | **0.291** | **0.63×** | **1.24×** |
| **100 000 × 100 000** | **0.253** | **0.597** | **0.42×** | **1.29×** |
| **500 000 × 500 000** | **2.620** | **4.398** | **0.60×** | **0.89×** |
| **1 000 000 × 1 000 000** | **4.828** | **9.568** | **0.50×** | **0.85×** |
| **2 000 000 × 2 000 000** | **12.535** | **22.636** | **0.55×** | **0.79×** |
| **5 000 000 × 5 000 000** | **19.746** | **63.223** | **0.31×** | **0.52×** |
| **10 000 000 × 10 000 000** | **43.395** | **207.4** | **0.21×** | **0.31×** |
| **20 000 000 × 20 000 000** | **92.591** | **275.6** | **0.34×** | **0.47×** |
| **50 000 000 × 50 000 000** | **430.9 (warm)** | **684.9** | **0.63×** | **0.89× (warm)** | ← cold single-iter 621.4 ms
| **100 000 000 × 100 000 000** | **958.1 (warm)** | **1 407** | **0.68×** | **0.91× (warm)** | ← cold single-iter 1 655 ms
| **200 000 000 × 200 000 000** | **2 121 (warm)** | **2 957** | **0.72×** | **0.97× (warm)** | ← cold single-iter 4 420 ms

Skewed (`a.size() >> b.size()`):

| digits | BigMath ms | GMP ms | BM/GMP | v13.0 BM/GMP |
|---|---:|---:|---:|---:|
| **100 000 × 10 000** | **0.139** | **0.300** | **0.46×** | **1.21×** |
| **500 000 × 50 000** | **0.569** | **2.105** | **0.27×** | **0.68×** |
| **1 000 000 × 100 000** | **1.370** | **4.727** | **0.29×** | **0.61×** |
| **2 000 000 × 200 000** | **3.394** | **10.029** | **0.34×** | **0.60×** |
| **5 000 000 × 500 000** | **8.902** | **31.633** | **0.28×** | **0.73×** |
| **10 000 000 × 1 000 000** | **21.659** | **73.174** | **0.30×** | **0.66×** |
| **20 000 000 × 2 000 000** | **73.793** | **175.9** | **0.42×** | **0.58×** |
| **50 000 000 × 5 000 000** | **204.1** | **676.8** | **0.30×** | **0.40×** |
| **100 000 000 × 10 000 000** | **448.8** | **999.4** | **0.45×** | **0.61×** |
| **200 000 000 × 20 000 000** | **989.2** | **1 930** | **0.51×** | **0.71×** |

Observations:

- The crossover against GMP moved from ~500k digits to ~10k. Below it the
  flag-chain limb kernels put Karatsuba within 1.4–1.5× of GMP's assembly
  basecase plus Toom tiers; at 10k digits the product already runs on the
  CRT NTT (balanced warm entry at 320 limbs per operand, was 640).
- 50k–20M digits: NEON tail kernels and the leaner Shoup butterfly on the
  transforms, the six-unit inverse, the parallel finalize pass and the
  thread-pool rewrite; just above each power-of-two transform length the
  truncated `H + L` scheme removes most of the old staircase.
- 50M+ cold rows include the one-time MFA plan build and first-touch page
  faults for multi-GB buffers, which is why the warm numbers are
  authoritative there. Same-session warm comparison against the previous
  `main`: 50M 602 → 431 ms, 100M 1288 → 958 ms, 200M 2782 → 2121 ms.

---

## Division

Skewed (`a.size() >> b.size()`), `bench_vs_gmp`, 2026-11-21:

| digits | BigMath ms | GMP ms | BM/GMP | 2026-11-20 | v13.0 |
|---|---:|---:|---:|---:|---:|
| 40 000 × 10 000 | 0.302 | 0.222 | 1.36× | 1.99× | 2.64× |
| 100 000 × 10 000 | 0.900 | 0.447 | 2.01× | 2.91× | 3.81× |
| **200 000 × 50 000** | **1.160** | **1.721** | **0.67×** | 1.02× | 2.31× |
| **500 000 × 100 000** | **2.532** | **4.533** | **0.56×** | 0.78× | 1.59× |
| **1 000 000 × 200 000** | **5.017** | **10.244** | **0.49×** | 0.61× | 1.44× |
| **2 000 000 × 500 000** | **8.901** | **25.087** | **0.35×** | 0.51× | 1.18× |
| **5 000 000 × 1 000 000** | **24.965** | **69.610** | **0.36×** | 0.48× | 1.00× |
| **10 000 000 × 2 000 000** | **55.673** | **157.1** | **0.35×** | 0.56× | 1.18× |
| **20 000 000 × 4 000 000** | **105.5** | **351.3** | **0.30×** | 0.40× | 0.82× |
| **50 000 000 × 10 000 000** | **244.1** | **1 307** | **0.19×** | 0.26× | 0.45× |
| **100 000 000 × 20 000 000** | **715.2** | **2 693** | **0.27×** | 0.32× | 0.49× |
| **200 000 000 × 40 000 000** | **1 450** | **4 472** | **0.32×** | 0.47× | 0.69× |

**What this table times.** The November 21 run called `a / b` (quotient only) on the
BigMath side and `mpz_tdiv_qr` (quotient and remainder) on the GMP side. The current harness reports matched Q/Q and QR/QR
rows with first-call and warm times for each operation on the same input. At these ratios the quotient is longer than the divisor, so
BigMath's quotient-only shortcut (#147) does not apply and both sides do a
full division; GMP's own quotient-only `mpz_tdiv_q` is 17–25% faster than
its `tdiv_qr` at ratio 2 in separate measurements, so against that entry
point the ratios above would be correspondingly less favourable. The balanced
equal-size rows are not comparable at all — the quotient is 1–2 limbs,
BigMath now computes it from the top limbs in 0.1–4 µs, and GMP is still
asked for the remainder — so that table stays omitted.

Like for like — quotient and remainder on both sides (`DivideAndRemainder`
vs `mpz_tdiv_qr`), dividend 5× the divisor, best of several runs; GMP column
measured 2026-11-20 on the same machine:

| dividend digits | before the pass ms | after ms | GMP ms | after / GMP |
|---:|---:|---:|---:|---:|
| 40 000 | 0.883 | 0.280 | 0.200 | 1.40× |
| **200 000** | 1.674 | **1.178** | **1.533** | **0.77×** |
| **1 000 000** | 5.330 | **4.444** | **9.789** | **0.45×** |
| **5 000 000** | 39.262 | **25.182** | **65.791** | **0.38×** |
| **20 000 000** | 141.0 | **104.4** | **344.2** | **0.30×** |

Observations:

- **Beats GMP from 200k-digit dividends.** The basecase, Burnikel–Ziegler
  and the Newton reciprocal were each rebuilt in this pass, and the
  dispatcher re-swept on a (divisor length, ratio) grid; details and
  per-algorithm numbers in [docs/DIVISION.md](docs/DIVISION.md).
- **The remaining gap is divisors of a few hundred limbs at high ratio**
  (100k ÷ 10k digits: 2.0×). Both Burnikel–Ziegler and Newton land at
  ~1.8× of GMP there; GMP's divide-and-conquer sits on an assembly basecase
  and Toom multiplies that have no counterpart here.
- At limb level the schoolbook basecase is 1.1–1.6× of GMP's (was 2.2–3.5×)
  and the single-limb kernel 1.6× (4.9 vs 3.1 ns per limb).

---

## Decimal parse (string → BigInteger)

| digits | BigMath ms | GMP ms | BM/GMP | v13.0 BM/GMP |
|---|---:|---:|---:|---:|
| 1 000 | 0.002 | 0.002 | 1.49× | 1.55× |
| 10 000 | 0.064 | 0.038 | 1.71× | 2.26× |
| 50 000 | 0.539 | 0.394 | 1.37× | 2.77× |
| **100 000** | **0.578** | **1.074** | **0.54×** | **1.62×** |
| **500 000** | **3.346** | **9.154** | **0.37×** | **0.69×** |
| **1 000 000** | **7.097** | **20.773** | **0.34×** | **0.60×** |
| **2 000 000** | **15.444** | **48.025** | **0.32×** | **0.56×** |
| **5 000 000** | **36.536** | **148.8** | **0.25×** | **0.44×** |
| **10 000 000** | **83.586** | **347.3** | **0.24×** | **0.40×** |
| **20 000 000** | **182.7** | **803.9** | **0.23×** | **0.38×** |
| **50 000 000** | **743.6** | **2 646** | **0.28×** | **0.53×** |

Parse beats GMP from 100k digits (the parallel leaf fan-out threshold; it was
500k). Below ~50k GMP's basecase edge shows.

---

## ToString (BigInteger → string)

Cold = single run including the divider-chain build (a per-process cache GMP
has no equivalent of); warm = steady state, measured at three sizes in this
run.

| digits | BigMath cold ms | BigMath warm ms | GMP ms | cold BM/GMP | 2026-11-20 cold | v13.0 BM/GMP |
|---|---:|---:|---:|---:|---:|---:|
| 1 000 | 0.007 | — | 0.004 | 1.80× | 1.83× | 1.83× |
| 10 000 | 0.191 | — | 0.082 | 2.33× | 2.41× | 3.37× |
| 50 000 | 1.614 | — | 0.897 | 1.80× | 1.76× | 2.81× |
| 100 000 | 3.244 | **1.358 (0.57×)** | 2.371 | 1.37× | 2.15× | 1.46× (warm) |
| **200 000** | **4.435** | — | **6.028** | **0.74×** | 1.33× | 2.27× (cold) |
| **500 000** | **10.586** | — | **20.429** | **0.52×** | 1.03× | 0.62× (warm) |
| **1 000 000** | **20.949** | **14.776 (0.30×)** | **49.560** | **0.42×** | 0.65× | 0.55× (warm) |
| **2 000 000** | **45.565** | **—** | **118.7** | **0.38×** | 0.58× | 0.49× (warm) |
| **5 000 000** | **168.5** | **121.8 (0.32×)** | **382.1** | **0.44×** | 0.65× | 0.58× (warm) |
| **10 000 000** | **362.0** | **—** | **900.8** | **0.40×** | 0.53× | 0.54× (warm) |
| **20 000 000** | **817.7** | **—** | **2 114** | **0.39×** | 0.47× | 0.49× (warm) |

ToString beats GMP cold from 200k digits (was 1M) and warm from 100k digits.
The cold gain is the division pass: the divider-chain build and the
divide-and-conquer splits are divisions. Warm times moved 3–6%.

---

## BigDecimal

GMP has no native BigDecimal; compared against `mpf_t` at matching precision
(`bits = total_digits × 3.322 + 64`). Sizes are `integer.fraction` digits.

| op | size | BigMath ms | GMP (mpf) ms | BM/GMP |
|---|---|---:|---:|---:|
| add | 100.10 | 0.000 | 0.000 | 3.05× |
| add | 1000.100 | 0.000 | 0.000 | 0.00× |
| add | 5000.500 | 0.001 | 0.000 | 4.33× |
| add | 20000.2000 | 0.002 | 0.000 | 4.00× |
| mul | 100.10 | 0.000 | 0.000 | 2.02× |
| mul | 1000.100 | 0.002 | 0.001 | 1.39× |
| mul | 5000.500 | 0.020 | 0.013 | 1.58× |
| mul | 20000.2000 | 0.060 | 0.095 | 0.63× |
| div | 100.10 / 10 dp | 0.000 | 0.000 | 3.66× |
| div | 1000.100 / 100 dp | 0.001 | 0.002 | 0.55× |
| div | 5000.500 / 500 dp | 0.010 | 0.024 | 0.41× |
| div | 20000.2000 / 2000 dp | 0.068 | 0.195 | 0.35× |
| parse | 100 digits | 0.000 | 0.001 | 0.83× |
| parse | 1000 digits | 0.005 | 0.004 | 1.03× |
| parse | 10000 digits | 0.087 | 0.105 | 0.83× |
| parse | 50000 digits | 0.706 | 1.004 | 0.70× |
| tostr | 100 digits | 0.000 | 0.000 | 0.64× |
| tostr | 1000 digits | 0.006 | 0.004 | 1.40× |
| tostr | 10000 digits | 0.196 | 0.101 | 1.94× |
| tostr | 50000 digits | 1.436 | 1.082 | 1.33× |

Scale-aware division, 2000.200-digit operand at varying target scales:

| target scale | BigMath ms | GMP (mpf) ms | BM/GMP |
|---|---:|---:|---:|
| 0 dp | 0.001 | 0.005 | 0.14× |
| 10 dp | 0.001 | 0.005 | 0.16× |
| 100 dp | 0.002 | 0.005 | 0.30× |
| 1000 dp | 0.007 | 0.009 | 0.82× |
| 5000 dp | 0.030 | 0.034 | 0.89× |

The division rows and the scale table are from the 2026-11-21 run (division
is 1.7–1.9× faster than on 2026-11-20); the other rows are from 2026-11-20.
Addition is allocation-bound at these sizes (≤2 µs absolute; the timer
resolution makes those ratios coarse).

---

## Internal algorithm crossovers

`multperf_simple 5` (BigMath-internal, equal-digit random operands, mean of 5):

| digits | Classical ms | Karatsuba ms | Karatsuba× | NTT ms | NTT× |
|---:|---:|---:|---:|---:|---:|
| 10 000 | 0.158 | 0.055 | 2.9× | 0.093 | 1.7× |
| 100 000 | 14.943 | 2.133 | 7.0× | 0.477 | 31.3× |
| 500 000 | 378.445 | 28.359 | 13.3× | 1.606 | 235.7× |

(v13.0: Classical 0.342 / 29.3 / 691 ms, Karatsuba 0.119 / 3.48 / 45.6 ms,
NTT 0.208 / 0.865 / 2.97 ms — all three algorithms are about 2× faster.)
The default 64-bit NTT path requires at least 640 total limbs, 64 limbs in
the shorter operand and an area of at least 96,000 limb pairs. With a parked
pool the area floor is 262,144. Balanced entry is 320 limbs per operand
when warm and 512 when cold.

---

## Optimization history (condensed)

The margins above are the product of the 2026-05/06 campaign and the two
2026-11 passes. Full analysis,
profiles, and rejection evidence live in the subsystem docs; per-release
summary in [CHANGELOG.md](CHANGELOG.md). The major steps, in order:

| when | change | headline effect |
|---|---|---|
| 2026-05 | 64-bit limbs (`BIGMATH_LIMB_64`), 3-prime CRT NTT, thread pool, M-G reciprocals, radix-4/8 fused butterflies, MFA layout | closed the GMP gap 3–5× across every band up to ~20M digits |
| #82–#95 | division/decimal session: bit-shift normalization, wrap-around Newton family (cyclic remainders, top-limbs quotient estimate, invertappr reciprocal), Newton band frontier, quotient-sized division, decimal D&C retunes | skewed div 5M×1M 2.87× → ~1.0×; killed the BZ non-pow2 5–60× blowups (97.7× at 2^18+1) |
| #96–#99 | **NEON Shoup butterflies** + interleaved twiddle tables + tail layers; dispatch retune they unlocked (NTT entry 5120 → 1280 limbs, Toom-3 window retired, CRT always-on) | butterfly kernel ~4×; mul 100k digits −40% |
| #103 | **fix: MFA inverse cross-twiddle ordering** — every ≥2^24-coefficient product had been silently wrong for 12 days; Newton's fixup fallback masked it as "stuck" runs. `mfa_roundtrip` added to ctest; all earlier 50M+ rows were invalid timings | 200M×40M div: hung → 6.3 s |
| #107–#114 | fused-MFA pass fusion (pointwise rides operand B's last forward stage), row-chunked `ParallelDo(6)` stages, **MFA gate 2^24 → 2^20**, cyclic products on the fused pipeline up to the 2^26 CRT ceiling, F3 on-the-fly packing | mul n=2^22 3×; div 100M÷20M 1.37× → 0.55×; 200M÷40M 6.5 s → 3.1 s |
| #116–#117 | **BZ odd-size bottom-padding fix** (odd divisor limb counts fell back to quadratic Knuth D — 1.3–10× on ~half of real shapes), basecase 512 → 128, full Newton frontier re-sweep | 12289-limb divisor ratio 1.5: 95 → 8.1 ms |
| #118–#119 | parallel D&C subtree fan-out for ToString and parse | both directions flipped to a 1.4–2.6× win ≥500k digits |
| #120–#124 | codebase audit: FastDivision MG top-equal qhat fix (silent wrong quotients on structured operands), remainder-sign fix, CRT-squaring dispatch (2.5–5.8× on Square, parse-cold −14…25%), dedups, header hygiene | correctness + Square/parse-cold wins |
| #129–#139 | **multiplication pass** (2026-11-20): flag-chain limb kernels (`addmul_1` 2.5–3.3 → 1.7 cycles/limb), NEON tail kernels + leaner Shoup butterfly, thread-pool rewrite (claim-based, spin-before-park; also fixes concurrent callers), NEON + parallel Garner finalize, six-unit inverse, truncated `H + L` transform length, area-rule dispatch; **fix: Karatsuba workspace heap overflow on skewed operands** | mul 1.5–4.5× at every size, GMP crossover 500k → 10k digits; div/parse/ToString 1.3–2.8× for free |
| #141–#149 | **division pass** (2026-11-21): raw-buffer Knuth D on a flag-chain `submul_1` (basecase 2.2–3.5× → 1.1–1.6× of GMP), in-place divide-and-conquer Burnikel–Ziegler (1.5–2.4×), Newton reciprocal on a guarded precision schedule with no refinement pass, dispatcher re-swept to one rule per regime, quotient-only division from the operand tops, single-limb kernel, schoolbook leading block in Newton; thread pool: pool-temperature-aware small-NTT entry and lock-free job pickup | division beats GMP from 200k-digit dividends (was 500k) and is 1.14–1.58× faster across the 12 published skewed cases; ToString cold beats GMP from 200k digits (was 1M) |

Two history figures are kept with their original data:

![MFA transpose fusion — multiply speed and fusion speedup vs operand size](docs/images/mfa_fusion_speedup.png)

*MFA transpose fusion (#78): fusing the transpose into the adjacent row-FFT
cut transform traffic from 8n to 4n bytes — +7–11% across the MFA band,
bit-exact.*

![Near-balanced division dispatch — Newton balanced band vs Burnikel-Ziegler](docs/images/division_balanced_speedup.png)

*Newton balanced band (#79): BZ's recursive halving lands intermediate NTT
multiplies just over power-of-2 transform boundaries on non-pow2 divisor
sizes (5–60× blowup, 97.7× at 2^18+1); Newton pads once and stays flat. The
band's floor and the residual ratio < 4/3 slice were later rebuilt by #88/#89
(quotient-sized division) and re-swept after #116/#117.*

Process lessons recorded along the way:

- Profile-sample arithmetic overcounts parallel-overlapped work — CPU-time
  savings ≠ wall-clock savings on the threaded stack (#84's rejection).
- The exact Newton tower drifts up to ~B^6 ulps; residue-window sizing needs
  ~B^16 headroom (a B^4 window caused an 8× ToString regression via silent
  FastDivision fallbacks).
- Warm-state-only measurement for chain-/plan-heavy ops; single-iter cold
  rows on those shapes scatter ±20–60% with process state.
- A documented "wall" deserves a kernel-level measurement before it is
  accepted: the sub-NTT gap to GMP was attributed to assembly for months, and
  cycles-per-limb on `addmul_1` showed it was the carry idiom (#136).
- Wall-clock is not CPU time: the finalize pass looked like 3–5% in summed
  profiles and was 17–33% of a threaded multiply's wall-clock (#134).
- There is no sanitizer job in CI. The Karatsuba workspace overflow (#137)
  shipped for months and was only caught when a later change made it crash.
- Random-input cross-checks miss structured-operand bug classes entirely —
  the seeded adversarial-limb fuzz (`test_adversarial_fuzz.cpp`, #125) now
  gates every ctest run.
