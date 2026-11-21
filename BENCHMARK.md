# BigMath vs GMP — full benchmark suite

Canonical run **2026-11-20** (post #129–#139: the multiplication pass — limb
kernels, NTT kernel/threading/truncated-length work, pool rewrite, dispatch
re-sweep). Single `bench_vs_gmp` pass, with warm steady-state re-measurement
for the chain-/plan-heavy rows where single-iteration cold timing is not
representative (see Methodology). Every table carries the previous canonical
ratio (2026-06-12, v13.0) for comparison.

![BigMath vs GMP ratio by operand size](docs/images/bigmath_vs_gmp.png)

*The figures were generated from the v13.0 run and have not been
regenerated; the tables below are current.*

**Headlines:**

- **Balanced multiplication beats GMP from 10k digits up** (was 500k) —
  0.21–0.63× from 50k to 20M digits, peaking at 4.8× faster (10M×10M digits:
  43 vs 207 ms), and 0.63–0.72× warm through 200M digits.
- **Skewed multiplication wins every measured shape** (0.27–0.51×), including
  100k×10k, which was 1.21×.
- **Skewed division beats GMP from 500k-digit dividends** (was 20M):
  0.26–0.78×; parity at 200k (1.02×).
- **Decimal I/O**: parse beats GMP from 100k digits (0.23–0.66×; was from
  500k); ToString from 1M digits cold (0.47–0.65×) and from 100k warm
  (0.31–0.63×).
- Below ~10k digits BigMath is now 1.4–1.5× behind GMP instead of 2.3–3.2×.
  The earlier "portable-C++ wall" was a loop idiom, not a language limit —
  see [docs/MULTIPLICATION.md](docs/MULTIPLICATION.md). What still trails:
  division below ~200k-digit dividends (2–2.9×), small decimal I/O
  (1.5–2.4×), and 1–4-limb products where the result allocation dominates.

For per-subsystem deep dives — algorithms, dispatch, optimization history,
rejected approaches — see [docs/MULTIPLICATION.md](docs/MULTIPLICATION.md),
[docs/DIVISION.md](docs/DIVISION.md),
[docs/STRING_CONVERSION.md](docs/STRING_CONVERSION.md),
[docs/BASE.md](docs/BASE.md). Release history: [CHANGELOG.md](CHANGELOG.md).
This file is the raw numbers, all subsystems in one place.

---

## Methodology

- **Hardware:** Apple M1 Max (10-core, 8 performance + 2 efficiency), 64 GB RAM.
- **OS:** macOS (Darwin 27.2.0).
- **Build:** CMake Release (`-O3 -march=native`), Apple Clang 21, C++20.
- **Reference:** GMP 6.3.0 (Homebrew), timed in-process on the same operands.
  GMP is single-threaded; BigMath uses its 8-thread pool on NTT-sized work.
- **Default stack:** `BIGMATH_LIMB_64=1` (64-bit limbs), `BIGMATH_NTT_CRT=1`
  (3-prime CRT NTT), `BIGMATH_USE_THREADS=1` (8-thread pool),
  `BIGMATH_NEON=1` on aarch64, MFA gate `2^20` coefficients.
- **Timing:** min over 1–1000+ iterations depending on size
  (`tests/performance/bench_vs_gmp.cpp`); rows ≥100k digits are single-iter.
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
  CRT NTT (entry at ~310 limbs per operand, was 640).
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

Skewed (`a.size() >> b.size()`) — the Newton/BZ/quotient-sized bands, real
algorithmic work:

| digits | BigMath ms | GMP ms | BM/GMP | v13.0 BM/GMP |
|---|---:|---:|---:|---:|
| 40 000 × 10 000 | 0.440 | 0.221 | 1.99× | 2.64× |
| 100 000 × 10 000 | 1.282 | 0.441 | 2.91× | 3.81× |
| 200 000 × 50 000 | 1.748 | 1.716 | 1.02× | 2.31× |
| **500 000 × 100 000** | **3.564** | **4.597** | **0.78×** | **1.59×** |
| **1 000 000 × 200 000** | **6.165** | **10.029** | **0.61×** | **1.44×** |
| **2 000 000 × 500 000** | **12.849** | **25.220** | **0.51×** | **1.18×** |
| **5 000 000 × 1 000 000** | **33.418** | **69.044** | **0.48×** | **1.00×** |
| **10 000 000 × 2 000 000** | **87.825** | **155.5** | **0.56×** | **1.18×** |
| **20 000 000 × 4 000 000** | **141.9** | **353.4** | **0.40×** | **0.82×** |
| **50 000 000 × 10 000 000** | **337.3** | **1 302** | **0.26×** | **0.45×** |
| **100 000 000 × 20 000 000** | **817.6** | **2 530** | **0.32×** | **0.49×** |
| **200 000 000 × 40 000 000** | **2 123** | **4 492** | **0.47×** | **0.69×** |

Balanced equal-size rows are degenerate — the quotient is 0–2 limbs and both
libraries short-circuit; absolute times are sub-15 µs through 1M digits and
the ratios are allocation noise, so the table is omitted.

Observations:

- **Beats GMP from 500k-digit dividends, parity at 200k.** No division code
  changed in this run: Newton and Burnikel–Ziegler inherit every
  multiplication win.
- The residual sub-200k-digit gap (2–2.9×) is the basecase band below the
  NTT-era algorithms. `FastDivision`'s `submul` loops still use the
  limb-carry idiom that the multiplication kernels dropped; the same
  restructuring applies and has not been done.
- The Newton/BZ dispatch frontier was tuned against the old multiply and has
  not been re-swept, so some of these shapes are likely not on their best
  algorithm.

---

## Decimal parse (string → BigInteger)

| digits | BigMath ms | GMP ms | BM/GMP | v13.0 BM/GMP |
|---|---:|---:|---:|---:|
| 1 000 | 0.002 | 0.002 | 1.51× | 1.55× |
| 10 000 | 0.064 | 0.038 | 1.70× | 2.26× |
| 50 000 | 0.620 | 0.392 | 1.58× | 2.77× |
| **100 000** | **0.701** | **1.062** | **0.66×** | **1.62×** |
| **500 000** | **3.243** | **9.159** | **0.35×** | **0.69×** |
| **1 000 000** | **7.238** | **21.990** | **0.33×** | **0.60×** |
| **2 000 000** | **18.899** | **47.509** | **0.40×** | **0.56×** |
| **5 000 000** | **36.469** | **146.9** | **0.25×** | **0.44×** |
| **10 000 000** | **80.995** | **346.6** | **0.23×** | **0.40×** |
| **20 000 000** | **183.9** | **815.9** | **0.23×** | **0.38×** |
| **50 000 000** | **707.6** | **2 655** | **0.27×** | **0.53×** |

Parse beats GMP from 100k digits (the parallel leaf fan-out threshold; it was
500k). Below ~50k GMP's basecase edge shows.

---

## ToString (BigInteger → string)

Cold = single run including the divider-chain build (a per-process cache GMP
has no equivalent of); warm = steady state, measured at three sizes in this
run.

| digits | BigMath cold ms | BigMath warm ms | GMP ms | cold BM/GMP | v13.0 BM/GMP |
|---|---:|---:|---:|---:|---:|
| 1 000 | 0.006 | — | 0.004 | 1.83× | 1.83× |
| 10 000 | 0.191 | — | 0.079 | 2.41× | 3.37× |
| 50 000 | 1.490 | — | 0.847 | 1.76× | 2.81× |
| **100 000** | **5.248** | **1.437 (0.63×)** | **2.441** | **2.15×** | **1.46× (warm)** |
| 200 000 | 8.403 | — | 6.306 | 1.33× | 2.27× (cold) |
| 500 000 | 21.500 | — | 20.894 | 1.03× | 0.62× (warm) |
| **1 000 000** | **33.193** | **15.360 (0.31×)** | **51.355** | **0.65×** | **0.55× (warm)** |
| **2 000 000** | **70.187** | **—** | **121.7** | **0.58×** | **0.49× (warm)** |
| **5 000 000** | **251.8** | **125.4 (0.33×)** | **387.2** | **0.65×** | **0.58× (warm)** |
| **10 000 000** | **478.5** | **—** | **909.9** | **0.53×** | **0.54× (warm)** |
| **20 000 000** | **1 010** | **—** | **2 151** | **0.47×** | **0.49× (warm)** |

ToString now beats GMP cold from 1M digits (it needed the warm state before)
and warm from 100k digits.

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
| div | 100.10 / 10 dp | 0.001 | 0.000 | 6.00× |
| div | 1000.100 / 100 dp | 0.002 | 0.002 | 1.00× |
| div | 5000.500 / 500 dp | 0.017 | 0.024 | 0.70× |
| div | 20000.2000 / 2000 dp | 0.128 | 0.195 | 0.66× |
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
| 0 dp | 0.001 | 0.005 | 0.16× |
| 10 dp | 0.002 | 0.005 | 0.36× |
| 100 dp | 0.003 | 0.005 | 0.56× |
| 1000 dp | 0.012 | 0.009 | 1.36× |
| 5000 dp | 0.053 | 0.034 | 1.54× |

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
Dispatch enters the NTT when `la·lb ≥ 96000` with the shorter operand at
least 64 limbs — balanced, ~310 limbs per operand (~6k digits).

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
| #129–#139 | **multiplication pass** (2026-11-20): flag-chain limb kernels (`addmul_1` 2.5–3.3 → 1.7 cycles/limb), NEON tail kernels + leaner Shoup butterfly, thread-pool rewrite (claim-based, spin-before-park; also fixes concurrent callers), NEON + parallel Garner finalize, six-unit inverse, truncated `H + L` transform length, area-rule dispatch; **fix: Karatsuba workspace heap overflow on skewed operands** | mul 1.5–4.5× at every size, GMP crossover 500k → 10k digits; div/parse/ToString 1.3–2.8× for free |

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
