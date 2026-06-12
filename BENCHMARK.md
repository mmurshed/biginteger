# BigMath vs GMP — full benchmark suite

Canonical run **2026-06-12 evening at v13.0** (post #120–#126: codebase audit
fixes, CRT-NTT squaring, dedup refactors, header hygiene). Single
`bench_vs_gmp` pass on a quiet machine, with warm steady-state re-measurement
for the chain-/plan-heavy rows where single-iteration cold timing is not
representative (see Methodology).

![BigMath vs GMP ratio by operand size](docs/images/bigmath_vs_gmp.png)

![Session before/after](docs/images/session_2026_06_11.png)

**Headlines:**

- **Balanced multiplication beats GMP at every size from 500k digits up** —
  peaking at 3.2× faster (10M×10M digits: 64 vs 209 ms) and staying at or
  below GMP through 200M digits warm (0.89–0.97×).
- **Skewed multiplication wins every measured shape from 500k×50k up**
  (0.40–0.73×).
- **Skewed division beats GMP from 20M-digit dividends (0.45–0.82×)**, parity
  at 5M (200M÷40M digits: 3.1 s vs GMP's 4.5 s).
- **Decimal I/O beats GMP from 500k digits in both directions**: parse
  0.38–0.69×, warm ToString 0.49–0.62×.
- Sub-NTT sizes (≲ 13k digits) remain 2–3× behind GMP's hand-tuned basecase —
  the documented portable-C++ wall.

For per-subsystem deep dives — algorithms, dispatch, optimization history,
rejected approaches — see [docs/MULTIPLICATION.md](docs/MULTIPLICATION.md),
[docs/DIVISION.md](docs/DIVISION.md),
[docs/STRING_CONVERSION.md](docs/STRING_CONVERSION.md),
[docs/BASE.md](docs/BASE.md). Release history: [CHANGELOG.md](CHANGELOG.md).
This file is the raw numbers, all subsystems in one place.

---

## Methodology

- **Hardware:** Apple M1 Max (10-core, 8 performance + 2 efficiency), 64 GB RAM.
- **OS:** macOS (Darwin 25.5.0).
- **Build:** CMake Release (`-O3 -march=native`), Apple Clang, C++20.
- **Reference:** GMP 6.3.0 (Homebrew), timed in-process on the same operands.
- **Default stack:** `BIGMATH_LIMB_64=1` (64-bit limbs), `BIGMATH_NTT_CRT=1`
  (3-prime CRT NTT), `BIGMATH_USE_THREADS=1` (8-thread pool),
  `BIGMATH_NEON=1` on aarch64, MFA gate `2^20` coefficients.
- **Timing:** min over 1–1000+ iterations depending on size
  (`tests/performance/bench_vs_gmp.cpp`); rows ≥100k digits are single-iter.
- **Warm steady state:** rows marked *warm* discard the first call (NTT plan
  build, divider-chain build, first-touch page faults) and take best-of-3.
  This is the authoritative comparison for ≥50M-digit multiplication and
  ≥100k-digit ToString: BigMath amortizes per-process caches that GMP doesn't
  have, so cold-vs-cold pays them once per process and never again.
- **Random inputs:** mt19937_64 seeded per case, full digit/limb range.

---

## Multiplication

Balanced (`a.size() == b.size()`):

| digits | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 1 000 × 1 000 | 0.002 | 0.001 | 2.33× |
| 5 000 × 5 000 | 0.033 | 0.013 | 2.56× |
| 10 000 × 10 000 | 0.100 | 0.031 | 3.19× |
| 50 000 × 50 000 | 0.363 | 0.292 | 1.24× |
| 100 000 × 100 000 | 0.797 | 0.617 | 1.29× |
| **500 000 × 500 000** | **4.081** | **4.595** | **0.89×** ← BigMath faster |
| **1 000 000 × 1 000 000** | **7.867** | **9.298** | **0.85×** |
| **2 000 000 × 2 000 000** | **16.606** | **21.030** | **0.79×** |
| **5 000 000 × 5 000 000** | **32.854** | **63.502** | **0.52×** ← 1.9× faster |
| **10 000 000 × 10 000 000** | **64.334** | **208.649** | **0.31×** ← 3.2× faster |
| **20 000 000 × 20 000 000** | **128.586** | **275.484** | **0.47×** ← 2.1× faster |
| **50 000 000 × 50 000 000** | **606.2 (warm)** | **684.5** | **0.89×** ← cold single-iter 1.08× |
| **100 000 000 × 100 000 000** | **1 302.5 (warm)** | **1 439.1** | **0.91×** ← cold 1.12× |
| **200 000 000 × 200 000 000** | **2 808.5 (warm)** | **2 907.3** | **0.97×** ← cold-state scatter up to 1.6× |

Skewed (`a.size() >> b.size()`):

| digits | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 100 000 × 10 000 | 0.365 | 0.303 | 1.21× |
| **500 000 × 50 000** | **1.463** | **2.149** | **0.68×** ← BigMath faster |
| **1 000 000 × 100 000** | **2.784** | **4.576** | **0.61×** |
| **2 000 000 × 200 000** | **5.694** | **9.414** | **0.60×** |
| **5 000 000 × 500 000** | **23.489** | **32.139** | **0.73×** |
| **10 000 000 × 1 000 000** | **46.937** | **70.970** | **0.66×** |
| **20 000 000 × 2 000 000** | **93.993** | **160.704** | **0.58×** |
| **50 000 000 × 5 000 000** | **271.616** | **677.097** | **0.40×** ← 2.5× faster |
| **100 000 000 × 10 000 000** | **620.034** | **1 019.765** | **0.61×** |
| **200 000 000 × 20 000 000** | **1 376.823** | **1 949.539** | **0.71×** |

Observations:

- The 500k–20M-digit balanced band is carried by the NEON Shoup butterflies +
  radix-8 fused chain + threaded 3-prime CRT; the fused-MFA pipeline
  (gate 2^20 coefficients) extends the win through 200M digits warm.
- Sub-50k digits stay on Karatsuba, where GMP's hand-tuned assembly basecase
  keeps a 2.3–3.2× lead; 50k–100k digits sit at 1.2–1.3×.
- 50M+ cold rows include the one-time MFA plan build and first-touch page
  faults for multi-GB buffers; run-to-run cold scatter at 200M digits reaches
  ±60%, which is why the warm numbers are authoritative there.

---

## Division

Skewed (`a.size() >> b.size()`) — the Newton/BZ/quotient-sized bands, real
algorithmic work:

| digits | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 40 000 × 10 000 | 0.576 | 0.218 | 2.64× |
| 100 000 × 10 000 | 1.754 | 0.460 | 3.81× |
| 200 000 × 50 000 | 4.054 | 1.755 | 2.31× |
| 500 000 × 100 000 | 7.671 | 4.810 | 1.59× |
| 1 000 000 × 200 000 | 14.465 | 10.056 | 1.44× |
| 2 000 000 × 500 000 | 28.920 | 24.466 | 1.18× |
| **5 000 000 × 1 000 000** | **69.446** | **69.387** | **1.00×** ← parity |
| 10 000 000 × 2 000 000 | 180.657 | 153.320 | 1.18× |
| **20 000 000 × 4 000 000** | **285.714** | **348.436** | **0.82×** ← BigMath faster |
| **50 000 000 × 10 000 000** | **589.797** | **1 317.839** | **0.45×** ← 2.2× faster |
| **100 000 000 × 20 000 000** | **1 263.489** | **2 581.606** | **0.49×** ← 2× faster |
| **200 000 000 × 40 000 000** | **3 119.348** | **4 543.900** | **0.69×** ← BigMath faster |

Balanced equal-size rows are degenerate — the quotient is 0–2 limbs and both
libraries short-circuit; absolute times are sub-15 µs through 1M digits and
the ratios are allocation noise, so the table is omitted.

Observations:

- **Beats GMP from 20M-digit dividends, parity at 5M.** Newton division
  inherits every multiplication win (its reciprocal iterations and wrapped
  remainders run cyclic mod-`B^L−1` products at half transform length on the
  fused MFA pipeline).
- The residual sub-200k-digit gap (2.3–3.8×) is the ≤520-limb-divisor
  basecase band below all NTT-era algorithms — the same portable-C++ wall as
  small multiplication, documented as accepted in
  [docs/DIVISION.md](docs/DIVISION.md).
- These canonical shapes land on even-ish divisor limb counts. Odd divisor
  sizes — ~half of real inputs — were 1.3–10× worse until the BZ odd-size
  padding fix (#116); they now track this table.

In-tree algorithm cross-check (`divperf_simple 5`, no GMP, limb shapes):
FastDivision vs BurnikelZiegler at `1024×512` 0.355/0.275 ms, `4096×2048`
7.99/2.20 ms (BZ 3.6×), `8192×512` 5.25/3.08 ms, `16384×512` 10.36/6.18 ms.

---

## Decimal parse (string → BigInteger)

| digits | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 1 000 | 0.002 | 0.002 | 1.55× |
| 10 000 | 0.087 | 0.038 | 2.26× |
| 50 000 | 1.081 | 0.390 | 2.77× |
| 100 000 | 1.739 | 1.072 | 1.62× |
| **500 000** | **6.234** | **8.989** | **0.69×** ← BigMath faster |
| **1 000 000** | **12.544** | **20.949** | **0.60×** |
| **2 000 000** | **26.787** | **48.148** | **0.56×** |
| **5 000 000** | **66.231** | **150.320** | **0.44×** |
| **10 000 000** | **141.762** | **350.371** | **0.40×** ← 2.5× faster |
| **20 000 000** | **306.340** | **814.925** | **0.38×** ← 2.7× faster |
| **50 000 000** | **1 451.183** | **2 744.997** | **0.53×** |

Parse beats GMP from 500k digits (the parallel leaf fan-out threshold sits at
100k digits; below it the serial divide-and-conquer numbers stand, and below
~50k GMP's basecase edge shows).

---

## ToString (BigInteger → string)

Cold = single run including the divider-chain build (a per-process cache GMP
has no equivalent of); warm = steady state, authoritative from 100k digits up.

| digits | BigMath cold ms | BigMath warm ms | GMP ms | warm BM/GMP |
|---|---:|---:|---:|---:|
| 1 000 | 0.006 | — | 0.003 | 1.83× |
| 10 000 | 0.262 | — | 0.078 | 3.37× |
| 50 000 | 2.509 | — | 0.893 | 2.81× |
| 100 000 | 7.904 | 3.41 | 2.337 | 1.46× |
| 200 000 | 14.422 | — | 6.342 | 2.27× (cold) |
| **500 000** | **35.363** | **13.5** | **21.884** | **0.62×** ← BigMath faster |
| **1 000 000** | **64.697** | **27.7** | **50.271** | **0.55×** |
| **2 000 000** | **128.493** | **58.4** | **119.721** | **0.49×** |
| **5 000 000** | **504.249** | **224.5** | **385.7** | **0.58×** |
| **10 000 000** | **913.196** | **488.3** | **901.2** | **0.54×** |
| **20 000 000** | **1 901.726** | **1 035.3** | **2 130.4** | **0.49×** ← faster even cold (0.88×) |

Warm ToString beats GMP from 500k digits (0.49–0.62×) via the parallel D&C
subtree fan-out (#118) on top of the cached Newton divider chains. The cold
chain build is the deferred one-shot-conversion lever; it only matters the
first time a process formats a given magnitude.

---

## BigDecimal

GMP has no native BigDecimal; compared against `mpf_t` at matching precision
(`bits = total_digits × 3.322 + 64`). Sizes are `integer.fraction` digits.

| op | size | BigMath ms | GMP (mpf) ms | BM/GMP |
|---|---|---:|---:|---:|
| add | 5000.500 | 0.001 | 0.000 | 5.0× |
| add | 20000.2000 | 0.002 | 0.001 | 4.3× |
| mul | 1000.100 | 0.002 | 0.001 | 2.23× |
| mul | 5000.500 | 0.031 | 0.013 | 2.46× |
| mul | 20000.2000 | 0.240 | 0.090 | 2.66× |
| div | 1000.100 / 100 dp | 0.002 | 0.002 | 1.16× |
| **div** | **5000.500 / 500 dp** | **0.020** | **0.024** | **0.82×** ← BigMath faster |
| div | 20000.2000 / 2000 dp | 0.207 | 0.197 | 1.05× |
| parse | 10 000 digits | 0.106 | 0.108 | 0.98× |
| parse | 50 000 digits | 1.064 | 1.004 | 1.06× |
| tostr | 10 000 digits | 0.261 | 0.109 | 2.41× |
| tostr | 50 000 digits | 2.439 | 1.156 | 2.11× |

Scale-aware division avoids unnecessary high-precision work — at small target
scales it runs up to 5× faster than `mpf_div` (2000.200-digit operand: 0 dp
0.18×, 10 dp 0.38×, 100 dp 0.65×); at deep scales (≥1000 dp) `mpf` leads
1.5–1.7×. Addition is allocation-bound at these sizes (≤2 µs absolute).
BigDecimal parse tracks GMP at ~1.0×; ToString carries the BigInteger
formatter's sub-100k-digit ratios.

---

## Internal algorithm crossovers

`multperf_simple 5` (BigMath-internal, equal-digit random operands, mean of 5):

| digits | Classical ms | Karatsuba ms | Karatsuba× | NTT ms | NTT× |
|---:|---:|---:|---:|---:|---:|
| 10 000 | 0.342 | 0.119 | 2.9× | 0.208 | 1.6× |
| 100 000 | 29.349 | 3.475 | 8.4× | 0.865 | 33.9× |
| 500 000 | 691.111 | 45.561 | 15.2× | 2.965 | 233.1× |

Dispatch enters NTT at 1280 total limbs (~13k digits in the default 64-bit
limb build); by 500k digits the CRT NTT is 233× over schoolbook and 15×
over Karatsuba.

---

## Optimization history (condensed)

The margins above are the product of the 2026-05/06 campaign. Full analysis,
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
- Random-input cross-checks miss structured-operand bug classes entirely —
  the seeded adversarial-limb fuzz (`test_adversarial_fuzz.cpp`, #125) now
  gates every ctest run.
