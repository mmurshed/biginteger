# BigMath vs GMP — full benchmark suite

Measured 2026-05-30 (fresh local rerun). Single-run snapshot capturing the state of `master` after the 2026-05 optimization pass (LIMB_64 default, multi-prime CRT NTT default, multithreaded NTT default, M-G div2by1 in `ClassicDivision`, M-G 3/2 qhat in `FastDivision` Base2_64, BZ Knuth-normalize fix, radix-4 + radix-8 fused NTT butterflies (PRs #59, #60), **Matrix Fourier Algorithm / Bailey 6-step CRT NTT retuned to n ≥ 2^24 coefficients**).

**All canonical tables re-measured 2026-06-12 evening (post PR #107–#111: MFA pass fusion, row-chunked stages, gate 2^24→2^20, cyclic MFA + 2^26 cap, on-the-fly packing).** Single-iter rows on chain-heavy ops (ToString ≥5M digits, balanced mul 500k–2M digits) vary ±20–60% run-to-run with process state — where it matters, rows are annotated with warm steady-state numbers, and the warm-state sections below are authoritative. Earlier note: **all ≤20M-digit tables previously re-measured 2026-06-12 morning on a quiet machine** after the two-day optimization run (PRs #82–#99): the division/decimal session (#82–#95 — bit-shift normalization, wraparound-Newton family, dispatch bands, quotient-sized division, cyclic gate) and the **NEON Shoup NTT** (#96–#99 — vectorized radix-8 butterflies + interleaved twiddle tables + tail layers, plus the dispatch retune they unlocked: NTT entry 5120 → 1280 limbs, Toom-3 retired, CRT always-on). 50M+ rows are annotated where not re-measured.

![BigMath vs GMP ratio by operand size](docs/images/bigmath_vs_gmp.png)

![Session before/after](docs/images/session_2026_06_11.png)

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
| 1 000 × 1 000 | 0.002 | 0.001 | 2.17× |
| 5 000 × 5 000 | 0.030 | 0.012 | 2.55× |
| 10 000 × 10 000 | 0.095 | 0.030 | 3.19× |
| 50 000 × 50 000 | 0.338 | 0.281 | 1.20× |
| 100 000 × 100 000 | 0.658 | 0.612 | 1.08× |
| **500 000 × 500 000** | **4.155** | **4.484** | **0.93×** ← single-iter cold; warm best 2.3 ms (0.55×) |
| **1 000 000 × 1 000 000** | **7.815** | **9.594** | **0.81×** ← BigMath faster |
| **2 000 000 × 2 000 000** | **16.580** | **20.897** | **0.79×** ← BigMath faster |
| **5 000 000 × 5 000 000** | **32.615** | **64.423** | **0.51×** ← BigMath 2× faster |
| **10 000 000 × 10 000 000** | **65.349** | **206.022** | **0.32×** ← BigMath 3.2× faster (was 0.49×) |
| **20 000 000 × 20 000 000** | **128.820** | **275.617** | **0.47×** ← BigMath 2.1× faster (was 1.05×) |
| **50 000 000 × 50 000 000** | **741.2** | **672.3** | **1.10×** ← cold; 0.90× warm (was 1.28×) |
| **100 000 000 × 100 000 000** | **1 574.3** | **1 420.7** | **1.11×** ← cold; 0.94× warm (was 1.30×) |
| **200 000 000 × 200 000 000** | **3 452.7** | **3 307.1** | **1.04×** ← cold; 0.97× warm (was 1.33×) |

Skewed (`a.size() >> b.size()`):

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 100 000 × 10 000 | 0.356 | 0.304 | 1.17× |
| **500 000 × 50 000** | **1.370** | **2.110** | **0.65×** ← BigMath faster |
| **1 000 000 × 100 000** | **2.459** | **4.571** | **0.54×** ← BigMath 1.9× faster |
| **2 000 000 × 200 000** | **5.741** | **9.509** | **0.60×** ← BigMath faster |
| **5 000 000 × 500 000** | **24.471** | **31.971** | **0.77×** ← BigMath faster |
| **10 000 000 × 1 000 000** | **47.700** | **72.288** | **0.66×** ← BigMath faster (was 1.29×) |
| **20 000 000 × 2 000 000** | **94.096** | **162.222** | **0.58×** ← BigMath 1.7× faster (was 1.51×) |
| **50 000 000 × 5 000 000** | **266.8** | **671.5** | **0.40×** ← BigMath 2.5× faster |
| **100 000 000 × 10 000 000** | **618.4** | **1 014.1** | **0.61×** ← BigMath faster |
| **200 000 000 × 20 000 000** | **1 323.4** | **1 960.3** | **0.68×** ← BigMath faster |

**Observations (updated 2026-06-12, post-NEON):**

- **BigMath beats GMP on balanced multiplication at every size from 500k digits up** — peaking at 3.2× (10M×10M: 65 vs 206 ms) and staying ≤1.0× warm through 200M digits. The NEON Shoup butterflies (PR #96), interleaved twiddle tables (#98), and tail layers (#99) compound with the radix-8 fused chain and multithreaded CRT; the fused-MFA work (#107–#111) extended the win through the former 20M–200M-digit gap.
- **Skewed multiplication beats GMP at every measured shape from 500k×50k up (0.40–0.77×)** — the former 10M×1M / 20M×2M losses (1.29×/1.51×) flipped to 0.66×/0.58× with the MFA gate retune.
- Sub-NTT sizes (≤ ~13k digits) stay on Karatsuba where GMP's hand-tuned basecase keeps a 2.2–3.3× lead; 50k–100k digits are near parity.
- **50M+ rows re-measured 2026-06-12 after the MFA correctness fix (PR #103)**: the MFA inverse cross-twiddle had been mis-ordered since 2026-05-31 (PR #72), so every pre-#103 measurement in this band was a wrong-result timing. Post-fix: balanced 50M–200M at 1.28–1.33×, skewed 50M–200M **beating GMP at 0.47–0.83×**. **Superseded in this band by PR #107 — see the next section: balanced 50M–200M now 0.90–0.97×.**

### MFA pass fusion + row-chunked stages (PR #107, 2026-06-12)

PR #107 fused the pointwise multiply into operand B's last forward MFA stage (`FusedForwardBMul` — fb's spectrum never reaches DRAM) and split the forward and inverse MFA stages into `ParallelDo(6)` (prime, half-row-range) work units; the old inverse batched whole per-prime transforms in `ParallelDo(3)`, idling cores for a third of the multiply. Net ≈1.20× in the MFA band — **BigMath now meets or beats GMP on balanced multiplication at every size from 500k digits up, and on every measured skew shape ≥500k digits.**

Warm-state methodology (NOT comparable to the cold single-iter rows above): raw-limb operands, first call discarded (plan build + first-touch page faults), best-of-3 warm × 3 interleaved rounds vs pre-#107 build, GMP timed identically in-process. M1 Max, load < 2.7 during measurement. `≈digits = limbs × 19.27`:

| shape | limbs | pre-#107 ms | post ms | GMP ms | BM/GMP was → now |
|---|---|---:|---:|---:|---|
| balanced ≈50M digits | 2.6M × 2.6M | 738.1 | **610.5** | 675.1 | 1.09× → **0.90×** |
| balanced ≈100M digits | 5.2M × 5.2M | 1 569.5 | **1 314.1** | 1 397.5 | 1.12× → **0.94×** |
| balanced ≈200M digits | 10.4M × 10.4M | 3 363.3 | **2 766.8** | 2 846.6 | 1.17× → **0.97×** |
| 10:1 skew ≈50M×5M | 2.6M × 260k | 311.5 | **260.9** | 666.1 | 0.47× → **0.39×** |
| 10:1 skew ≈100M×10M | 5.2M × 520k | 746.7 | **618.7** | 996.5 | 0.75× → **0.62×** |
| 10:1 skew ≈200M×20M | 10.4M × 1.04M | 1 581.9 | **1 319.3** | 1 893.7 | 0.84× → **0.70×** |
| div ≈100M÷20M digits | 5.2M ÷ 1.04M | 3 091.0 | 3 151.9 | 2 305.4 | 1.34× → 1.37× (no change) |
| div ≈200M÷40M digits | 10.4M ÷ 2.08M | 6 109.7 | 6 319.3 | 4 415.1 | 1.36× → 1.43× (no change) |

Two structural findings from this work:

- **The band is no longer purely memory-bandwidth-bound post-NEON.** A 2-concurrent-process probe shows one multiply draws ~64–69% of the M1 Max's DRAM bandwidth; pure pass-cutting that sacrificed work-unit concurrency regressed 1.47×. Parallelism is the current lever (mem_pass_fusion.md carries the revised analysis).
- **Division did not inherit the win** (~1% movement): Newton's inner products route through the cyclic `MultiplyMod2km1` path (capped at `n ≤ 2^22`, below the MFA gate by design) and transforms at `n = 2^22–2^23` — below the 2^24 MFA gate. Addressed by the gate retune below.

### MFA gate retune 2^24 → 2^20 (2026-06-12, post-#107)

A `sample` profile of 200M÷40M-digit division showed ~88% of active compute in plain (whole-transform) radix-8 NTT layers with half of all thread-time in `psynch_cvwait` — Newton's transforms sit at `n = 2^22–2^23`, below the MFA gate, on the non-MFA path whose `ParallelDo(6)`/`ParallelDo(3)` whole-transform units idle cores. The 2^24 break-even was measured before PR #107's row-chunked stages existed; re-sweeping the gate (warm-state, quiet machine, best-of-2×2–3 interleaved rounds) flipped it decisively:

| op | size (limbs) | n | gate 2^24 | gate 2^20 | BM/GMP was → now |
|---|---|---|---:|---:|---|
| mul | 1M × 1M | 2^22 | 299.1 | **100.7** | 1.12× → **0.38×** |
| mul | 2M × 2M | 2^23 | 621.5 | **278.9** | 1.15× → **0.52×** |
| mul | 500k × 500k | 2^21 | 102.5 | **50.4** | 1.06× → **0.51×** |
| mul | 250k × 250k | 2^20 | 25.9 | 25.5 | 0.59× (wash) |
| div | 1M ÷ 200k | — | 358.6 | **288.5** | 1.26× → **1.01×** |
| div | 2.6M ÷ 520k | — | 904.0 | **626.8** | 1.03× → **0.71×** |
| div | 5.2M ÷ 1.04M | — | 3 166.0 | **1 700.9** | 1.37× → **0.74×** |
| div | 10.4M ÷ 2.08M | — | 6 492.4 | **4 426.9** | 1.47× → **1.00×** |

**Division now beats GMP at ≈50–100M-digit dividends (0.71–0.74×) and reaches parity at 200M÷40M digits.** Lower gates (2^18) measured ≈ wash vs 2^20 (≤4% div); 2^20 keeps the gate out of the latency-sensitive sub-ms band. The remaining 200M÷40M gap is the cyclic `MultiplyMod2km1` transforms (always plain, `n ≤ 2^22` cap) — addressed in the next section.

Correctness: raw-limb GMP probe OK at every newly-MFA size (n=2^20–2^23 balanced + 8:1 skew at n=2^19 high-skew gate); full unit/roundtrip/mult/div correctness batteries green.

### Cyclic-product MFA routing + cap raise (2026-06-12)

`MultiplyMod2km1` — the wrap-around product Newton division uses for its reciprocal iterations and wrapped remainders — always ran plain whole-transform `ParallelDo(6/3)` units and was capped at `n ≤ 2^22`. The cap predated MFA support; the MFA permutation in fact cancels between forward and pointwise (both operands share it) and the inverse returns natural order, so the cyclic path now runs the same fused pipeline as the linear multiply (helpers `MfaFusedForwardPointwise` / `MfaFusedInverse`, extracted from `Multiply`) at `n ≥` the MFA gate, and the cap rose to the CRT ceiling `2^26` (coefficient-sum headroom `N·2^64 = 2^90 < p1·p2·p3 ≈ 2^90.5`). Newton's two `nCyc` gates rose to match — large iterations again use cyclic transforms at HALF the linear product's length.

Warm-state, best-of-3 × 3 interleaved rounds vs post-gate-retune main, quiet machine:

| div (limbs) | ≈digits | before ms | after ms | BM/GMP was → now |
|---|---|---:|---:|---|
| 1M ÷ 200k | 20M÷4M | 288 | 290 | ~1.0× (wash) |
| 2.6M ÷ 520k | 50M÷10M | 617 | 603 | 0.69× → **0.67×** |
| 5.2M ÷ 1.04M | 100M÷20M | 1 668 | **1 301** | 0.71× → **0.55×** |
| 10.4M ÷ 2.08M | 200M÷40M | 4 194 | **3 074** | 0.95× → **0.68×** |

**Division now beats GMP at every measured size from 20M-digit dividends up.** Cumulative across PR #107 + the gate retune + this change, 200M÷40M-digit division went 6.5s → 3.1s (1.47× → 0.68× vs GMP, a 2.1× swing).

Correctness: probe extended with cyclic rows (vs GMP fold-mod `2^{64L}−1`) at n=2^19 (plain), 2^20/2^22 (MFA), 2^23/2^24 (beyond the old cap, incl. uneven operands) and full-pipeline division quotient checks at both big tiers — all green; transient-break validated (wrong inverse planes → the four MFA-routed cyclic rows FAIL, the plain row stays OK); full unit/roundtrip/mult/div batteries green.

### F3 on-the-fly packing (PR #111, 2026-06-12)

In the fused MFA window, stage A's gather packs each element straight from the operand's limbs (`PackedOperandView` / `FusedForwardAPack`) — the pre-packed coefficient planes, their padding zero-fills, and the `fb` planes disappear from the fused path entirely, and the `% P` work moves into the parallel gather. Warm-state: 2–5% across the band (mul 1M limbs −4.4% → 0.36× vs GMP; div 100M÷20M digits −4.6% → 0.53×; 200M÷40M −3.8% → 0.66×).

**Bandwidth saturation tracking** (2-concurrent-process probe at 10.4M-limb mul, per-process slowdown / aggregate throughput): pre-#107 +28% / 1.56× (~64% of DRAM bandwidth used by one multiply) → post-#109 +38% / 1.45× (~69%) → post-#111 +45% / 1.38× (~72%). Still unsaturated — parallelism-first remains the right bias; pure byte-cutting (F2 finalize fusion) stays deprioritized until this approaches ~90%.

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

### MFA transpose fusion (PR #78, 2026-05-31)

The MFA forward/inverse transforms ran each axis as a full out-of-place transpose immediately followed by a full row-FFT sweep over the just-written scratch buffer — two extra round-trips through DRAM per axis. At MFA sizes that scratch is 64–256 MB per prime, and the large-MFA kernel is memory-bandwidth bound on this machine, so those round-trips dominate. Fusing the transpose into the adjacent row-FFT — gather a tile of rows into an L2-resident buffer, FFT (and cross-twiddle) there, write the result once — cuts transform traffic from `8n` to `4n` bytes per forward/inverse. The per-row FFT math is unchanged, so output is **bit-exact**: verified against the unfused path (FNV checksum of product limbs) at 2M–8M limbs, plus `mult_correctness` (18/0). Gated by `BIGMATH_NTT_MFA_FUSE` (default on), tile rows by `BIGMATH_NTT_MFA_FUSE_TILE` (default 16).

![MFA transpose fusion — multiply speed and fusion speedup vs operand size](docs/images/mfa_fusion_speedup.png)

Full sweep, balanced Base2_64, M1 Max, best-of-3 interleaved (fusion default vs `-DBIGMATH_NTT_MFA_FUSE=0`). `≈digits = limbs × 19.27`:

| limbs/operand | ≈digits | path | baseline ms | fusion ms | speedup |
|---:|---:|---|---:|---:|---:|
| 3 000 | 58K | non-MFA | 0.51 | 0.54 | −5.0% |
| 10 000 | 193K | non-MFA | 1.89 | 1.92 | −1.5% |
| 50 000 | 963K | non-MFA | 8.69 | 9.08 | −4.2% |
| 200 000 | 3.9M | non-MFA | 40.02 | 40.06 | −0.1% |
| 700 000 | 13.5M | non-MFA | 250.91 | 241.81 | +3.8% |
| 1 300 000 | 25M | non-MFA | 648.92 | 649.04 | −0.0% |
| 2 000 000 | 38.5M | MFA (gate) | 689.23 | 678.93 | +1.5% |
| **2 600 000** | **50M** | **MFA** | **1 121.25** | **1 011.62** | **+10.8%** |
| **4 000 000** | **77M** | **MFA** | **1 159.97** | **1 046.13** | **+10.9%** |
| **5 200 000** | **100M** | **MFA** | **2 379.93** | **2 163.77** | **+10.0%** |
| **6 000 000** | **116M** | **MFA** | **2 393.66** | **2 181.19** | **+9.7%** |
| **8 000 000** | **154M** | **MFA** | **2 473.24** | **2 231.99** | **+10.8%** |
| **10 000 000** | **193M** | **MFA** | **4 978.10** | **4 636.81** | **+7.4%** |

The win switches on exactly at the MFA gate (~2–2.6M limbs / ~40–50M digits) and holds **+7–11%** through 193M digits. Below the gate the fused code path is not reached, so those rows run byte-identical machine code — the ±5% scatter is sub-millisecond measurement noise, not a regression.

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
| 1 000 × 1 000 | <0.001 | <0.001 | — |
| 5 000 × 5 000 | 0.001 | <0.001 | — |
| 10 000 × 10 000 | <0.001 | <0.001 | — |
| 50 000 × 50 000 | 0.001 | <0.001 | 2.33× |
| 100 000 × 100 000 | 0.001 | 0.001 | 2.31× |
| 500 000 × 500 000 | 0.007 | 0.004 | 1.59× |
| 1 000 000 × 1 000 000 | 0.014 | 0.009 | 1.62× |
| 5 000 000 × 5 000 000 | 1.355 | 0.162 | 8.39× |

Skewed (`a.size() >> b.size()`) — Newton/BZ band, real algorithmic work:

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 40 000 × 10 000 | 0.741 | 0.221 | 3.35× |
| 100 000 × 10 000 | 2.161 | 0.446 | 4.84× |
| **200 000 × 50 000** | **3.825** | **1.691** | **2.26×** ← was 6.65× at session start |
| **500 000 × 100 000** | **7.647** | **4.629** | **1.65×** ← was 3.85× |
| **1 000 000 × 200 000** | **14.607** | **10.018** | **1.46×** ← was 3.36× |
| **2 000 000 × 500 000** | **28.387** | **24.623** | **1.15×** ← was 3.36× |
| **5 000 000 × 1 000 000** | **70.481** | **69.836** | **1.01×** ← parity, was 2.87× |
| **10 000 000 × 2 000 000** | **173.530** | **152.696** | **1.14×** ← was 2.79× |
| **20 000 000 × 4 000 000** | **282.254** | **346.173** | **0.82×** ← BigMath faster, was 2.78× |
| **50 000 000 × 10 000 000** | **587.0** | **1 306.9** | **0.45×** ← BigMath 2.2× faster (was 1.05×) |
| **100 000 000 × 20 000 000** | **1 255.4** | **2 550.1** | **0.49×** ← BigMath 2× faster (was 1.23×) |
| **200 000 000 × 40 000 000** | **3 010.4** | **4 491.6** | **0.67×** ← BigMath faster (was 1.37×; 2.85× at session start) |

**Observations (updated 2026-06-12, post-NEON):**

- **Skewed division beats GMP from 20M-digit dividends up (0.45–0.82×) and is at parity from 5M.** At the session start the ≥5M rows sat at 2.8–2.9×; the morning run brought them to parity, and the MFA work (#109/#110: gate retune + cyclic fused routing) flipped 50M–200M to a 1.5–2.2× BigMath win. The wraparound-Newton family (PRs #85–#87) halved the band first (serial-transform-chain cuts: cyclic mod-B^L−1 products, top-limbs quotient estimates, invertappr reciprocal); the NEON butterflies then accelerated every NTT inside Newton.
- The former 200k×50k worst point (6.65×) sits at 2.23× after the band retunes (#92: `NEWTON_MEDIUM` 5/2 @ 2560 + 8/5 band @ 6144 — fractional ratios kill one-limb knife-edge cliffs) plus NEON. The residual worst point is now 100k×10k (4.92×, 520-limb divisor below all NTT-era bands).
- Balanced equal-size rows are degenerate (quotient 0–1 limbs, both libraries short-circuit; sub-15µs absolute through 1M digits). FastDivision's scalar-normalization fix (PR #82, bit-shift normalize, 3.5× on that path) shows up in driver benchmarks rather than these noise-level rows.
- The 200M×40M row's history is a story: the 2026-05-30 value (12.98 s, 2.85×) and every later attempt ran the **broken MFA band** (PR #103) — re-measure runs "hung" because Newton's fixup cap detected the corrupt products and fell back to quadratic FastDivision. Post-fix it completes normally at **6.26 s, 1.37×**.

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

### Near-balanced division — Newton balanced band (2026-05-31)

Profiling the dispatch path on near-balanced (ratio ≈ 2) division at large
divisor sizes surfaced two compounding problems, both fixed in
`perf/newton-balanced-division`:

1. **Newton single-block pathology.** The single-block path (`na ≤ 2n+1`) ran a
   `2n+1`-limb chunk through the truncated `(chunk·R) >> 2n` quotient estimate.
   The truncation error scales with `chunk/B^(2n)`, which reaches ~`B` once the
   chunk exceeds `2n` limbs (the `+1` limb routinely appears from the Knuth
   normalize shift on `a ≈ 2n`). The estimate then underestimated `Q` by more
   than the 8-step fixup cap and bailed to quadratic `FastDivision`. At
   `nb = 50000` this was a 23× spike (5200 ms vs ~360 ms at neighbouring sizes).
   Fix: route `na > 2n` through the blockwise path so every chunk stays `≤ 2n`.

2. **Burnikel-Ziegler non-power-of-2 blowup.** BZ's recursive `2n/n` halving
   lands its intermediate NTT multiplies just over power-of-2 length boundaries
   for non-power-of-2 divisor sizes, where the transform length doubles. The
   constant factor compounds across the recursion depth into a **5–60× slowdown**
   versus Newton, worst right above a power of two (`n = 2^k + 1`). Newton pads
   once to the working size and stays flat. Fix: a Newton *balanced band* —
   `b ≥ NEWTON_BALANCED_B (98304)` and `a ≥ 2b` — takes the near-balanced large
   case that BZ previously owned. (Exact-power-of-2 divisor sizes, BZ's best
   case where it ties Newton, regress ~4 %; they are rare in practice.)

Dispatch wall-clock, BigMath Base2_64, M1 Max, before (`main` `d38f3c8`, → BZ)
vs after (Newton balanced band):

| shape (limbs) | divisor | before ms | after ms | speedup |
|---|---|---:|---:|---:|
| `200000/100000` | non-pow2 | 496.6 | 249.9 | **1.99×** |
| `500000/250000` | non-pow2 | 5 177.7 | 588.3 | **8.80×** |
| `1000000/500000` | non-pow2 | 10 554.7 | 1 485.7 | **7.10×** |
| `2000000/1000000` | non-pow2 | 21 006.9 | 3 747.6 | **5.61×** |
| `524288/262144` | pow2 (BZ best) | 656.4 | 684.0 | 0.96× |
| `524290/262145` | 2¹⁸+1 (BZ worst) | 92 060.5 | 941.9 | **97.7×** |
| `3000000/1000000` | ratio 3 (control) | 4 687.7 | 4 635.7 | 1.01× |

![Near-balanced division dispatch — Newton balanced band vs Burnikel-Ziegler](docs/images/division_balanced_speedup.png)

Harness: `tests/performance/division_balanced_bench.cpp`. Plot regenerated by
`docs/images/make_division_balanced_plot.py`. Cross-checked against BZ
limb-for-limb across `na ∈ {2n-1, 2n, 2n+1, 2n+2, 3n}` × `nb` near the boundary
(120 cases) and the canonical `div_correctness` suite — all match.

---

## 2026-06-11/12 optimization run (PRs #82-#99)

Two-day profiling-driven pass: division and decimal I/O first, then NEON across the CRT NTT.
Full analysis and rejection notes in [docs/DIVISION.md](docs/DIVISION.md),
[docs/STRING_CONVERSION.md](docs/STRING_CONVERSION.md), and
[docs/MULTIPLICATION.md](docs/MULTIPLICATION.md); this is the headline summary. The refreshed
tables above are the post-run state.

![Session before/after](docs/images/session_2026_06_11.png)

| PR | change | headline win |
|---|---|---|
| #82 | FastDivision bit-shift normalization (replaces scalar-d normalize; kills per-limb `__udivmodti4` remainder denorm) | balanced div 1M digits 0.89 → 0.25 ms (3.5×) |
| #83 | decimal D&C thresholds: parse 8192 → 2048, tostr 2048 → 1024 (May sweep inverted by the radix-8/MFA/Newton changes) | parse −12…−53% (4k–1M digits), tostr −19…−33% (1.5k–10k) |
| #84 | doc-only: prepared-transform reuse in `DivideChunk` REJECTED — cut CPU 11% but wall-clock flat (the 6 forward transforms already run concurrently; only serial-chain cuts convert to wall-clock on the threaded stack) | — |
| #85 | cyclic wrap-around remainder: `Q·b_norm` replaced by its residue mod `B^L−1` at half transform length | skewed div −10–11%, tostr −7% |
| #86 | top-limbs quotient estimate (GMP `mu_divappr` style): multiply only top n+1 chunk limbs against R | skewed div −20% at transform-boundary sizes, tostr −11% |
| #87 | invertappr-style wrapped Newton iteration in `ApproxReciprocal` (exact E from cyclic residue + top-slice correction) | one-shot skewed div −24–25%, `Divider` setup −29% |
| #88 | Newton balanced band ratio 2/1 → 4/3 (post-#85–87 crossover re-measured) | nb=131073 limbs ratio 1.5: 10.7 s → 157 ms |
| #89 | **quotient-sized division** (`QuotientSizedDivision.h`): divides operand tops for short quotients; cost scales with quotient, not divisor; balanced-band floor 98304 → 24576 limbs | 2^k+1 pathology ratios 1.05–1.25: 1.07–5.35 s → 28–71 ms (38–75×) |
| #92 | medium band 3/1 @ 4096 → **5/2 @ 2560** + new **8/5 band @ 6144** (fractional ratios kill the one-limb knife-edges where digit-derived operands fall off the band onto BZ's non-pow2 blowups) | 200k×50k digits 6.65× → 3.76×; cliff shapes 300k×100k 24.6× → 3.07×, 400k×200k 20.7× → 3.57× |
| #94 | thin-quotient extension of the quotient-sized band (`b ≥ 8192`, `Δ ≤ b/8`) | 2^14+1-limb divisor ratio 1.1: 66.4 → 8.3 ms |
| #95 | separate gate for cyclic NTT products (`BIGMATH_CYCLIC_NTT_THRESHOLD = 1280` vs the 5120 full-product threshold — cyclic runs at half transform length) | cold tostr 100k ratio −30% |
| #96 | **NEON radix-8 butterflies via Shoup multiplication** (30-bit CRT primes in 32-bit lanes; `Plan` carries `floor(w·2^32/P)` companions; bit-exact) | butterfly kernel 4×; mul 100k −40% |
| #97 | dispatch retune the NEON win unlocked: NTT entry 5120 → **1280** limbs, Toom-3 window retired, CRT gate 5000 → 256 (always-CRT) | mul 30k digits 2.3×; tostr 1M −25% |
| #98 | **interleaved (twiddle, Shoup) tables** — one cache line per pair; removes the n > 2^20 memory-bound regression and the size gate | mul 10M −13%, 1M −30% vs separate tables |
| #99 | NEON tail layers (radix-4/2 collapse to one broadcast Shoup pair + `vld4q`) and inverse 1/n scaling | mul 10M −10% |
| #101 | Newton band floors post-NEON: high-skew 2048 → 768, medium 2560 → 1024, 8/5 band 6144 → 4096 | 1536-limb divisor ratio 10: 12.0 → 4.9 ms |
| #102 | ToString chain-top rounding (1/16-octave grid; cache no longer keyed on exact digit count) + Karatsuba leaf threshold 48 → 32; 2-row leaf unrolling rejected (M1 OoO already hides the carry chain) | mixed-size ~100k workload −34% |
| #103 | **fix: MFA inverse cross-twiddle ordering** — mis-fused since PR #72 (2026-05-31); every ≥2^24-coeff product was silently wrong for 12 days; Newton's fixup fallback masked it as "stuck" runs. `mfa_roundtrip` now in ctest. | 200M×40M div: hung → 6.26 s (1.37×) |

Cumulative on the standard grid across the run: balanced mul 5M×5M 48.9 → **26.2 ms (0.41× vs
GMP — 2.4× faster than GMP)**; skewed div 5M×1M 200.6 → **69.3 ms (2.87× → 1.03×)** and 10M×2M
at **0.94× — BigMath faster than GMP**; tostr 1M 224.8 → **101.5 ms**; tostr 20M 2.52× → **1.16×**;
parse 20M 1.50× → **1.09×**.

Division dispatch now has no known residual shape holes: ratio ∈ (1, 4/3) at b ≥ 24576 limbs routes to
quotient-sized division (which also beats BZ at BZ's exact-power-of-2 best case), 4/3 ≤ ratio bands route
to Newton, small/degenerate shapes keep FastDivision.

Two recorded process lessons: profile-sample arithmetic overcounts parallel-overlapped work (CPU-time
savings ≠ wall-clock savings — see PR #84's rejection), and the exact Newton tower drifts up to ~B^6 ulps,
so any residue-window sizing needs ~B^16 headroom (a B^4 window caused an 8× ToString regression via
silent FastDivision fallbacks before being caught with a `10^500000`-divisor repro).

---

## In-tree simple harnesses (k-averaged, no GMP)

`tests/multperf_simple.cpp` and `tests/divperf_simple.cpp` are the standalone
in-tree harnesses (no GMP dependency). Each compares BigMath's own algorithms
against each other and validates the results limb-for-limb. Both take an
optional run-count `k` (default 3); every size is timed `k` times and the
**arithmetic mean** is reported:

```
multperf_simple 5      # avg of 5 runs per size
divperf_simple 5
```

Measured 2026-05-31, M1 Max, `-O2`, default stack (Base2_32 internal repr,
`k = 5`).

### Multiplication — `multperf_simple 5`

Two random decimal numbers of equal digit length. Speedup is vs Classical.

| digits | Classical ms | Karatsuba ms | Karatsuba× | NTT ms | NTT× |
|---:|---:|---:|---:|---:|---:|
| 100 | 0.000 | 0.000 | 0.5× | 0.007 | — |
| 10 000 | 0.274 | 0.094 | 2.9× | 0.391 | 0.7× |
| 100 000 | 28.917 | 4.147 | 7.0× | 1.143 | 25.3× |
| 500 000 | 692.964 | 45.330 | 15.3× | 4.585 | 151.1× |

At 100 digits both classical and Karatsuba are sub-microsecond (rounds to
0.000 ms); NTT loses to schoolbook at that size, as expected. Karatsuba and NTT
cross over between 10k and 100k digits — by 500k digits NTT is **151× faster**
than schoolbook. All three products matched exactly.

### Division — `divperf_simple 5`

`FastDivision` (Knuth D variant) vs `BurnikelZieglerDivision`. Speedup is
`fast / bz` (>1 means BZ faster). Results cross-checked limb-for-limb.

| shape (limbs) | FastDivision ms | BurnikelZiegler ms | fast/bz |
|---|---:|---:|---:|
| 1024 × 512 | 0.845 | 0.864 | 0.98× |
| 4096 × 2048 | 7.833 | 3.218 | 2.43× |
| 8192 × 512 | 5.154 | 5.075 | 1.02× |
| 16384 × 512 | 10.306 | 10.311 | 1.00× |

BZ's recursive halving pays off on the near-balanced `4096 × 2048` shape
(2.43× over FastDivision). On the skewed `n × 512` shapes the two are at parity
— the divisor is small enough that FastDivision's quadratic inner loop is short.

---

## Decimal parse (string → BigInteger)

| size (digits) | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 1 000 | 0.003 | 0.002 | 1.60× |
| 10 000 | 0.087 | 0.038 | 2.26× |
| 50 000 | 1.081 | 0.393 | 2.75× |
| 100 000 | 2.431 | 1.052 | 2.31× |
| 500 000 | 15.054 | 8.907 | 1.69× |
| 1 000 000 | 32.478 | 20.621 | 1.58× |
| 2 000 000 | 69.853 | 47.904 | 1.46× |
| 5 000 000 | 177.393 | 149.434 | 1.19× |
| **10 000 000** | **380.508** | **345.231** | **1.10×** |
| **20 000 000** | **813.649** | **803.663** | **1.01×** ← parity |
| **50 000 000** | **2 766.0** | **2 659.0** | **1.04×** ← near parity (was 1.38×) |

**Observation (updated 2026-06-12 evening, post #107–#111):** parse is at **1.01–1.19× of GMP from 5M digits up** — the 20M and 50M rows reached parity once the MFA gate retune (#109) put parse's internal multiplies on the fused path. The `DecimalDcThreshold` retune (8192 → 2048, PR #83) shaved 10–26% through the 10k–2M band.

---

## ToString (BigInteger → string)

| size (digits) | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 1 000 | 0.006 | 0.003 | 1.84× |
| 10 000 | 0.262 | 0.078 | 3.38× |
| 50 000 | 2.395 | 0.856 | 2.80× |
| **100 000** | **9.657** | **2.406** | **4.01×** ← was 8.40× at session start |
| **200 000** | **20.001** | **6.110** | **3.27×** ← was 6.35× |
| **500 000** | **53.410** | **20.777** | **2.57×** ← was 5.15× |
| **1 000 000** | **105.587** | **48.963** | **2.16×** ← was 4.48× |
| **2 000 000** | **220.377** | **119.129** | **1.85×** ← was 4.03× |
| **5 000 000** | **830.148** | **398.580** | **2.08×** ← cold-chain; 510 ms warm (1.28×) |
| **10 000 000** | **1 536.425** | **908.349** | **1.69×** ← cold-chain; 1 062 ms warm (1.17×) |
| **20 000 000** | **3 156.844** | **2 149.992** | **1.47×** ← cold-chain (was 1.16× in a chain-warm process) |

**Observation (updated 2026-06-12 evening):** ToString is divider-chain-bound, so it compounds every division and multiplication win of the run: 1M digits went 224.8 → **106 ms** across the two days (4.48× → **2.16×** vs GMP) and 10M digits is **1 062 ms warm (1.17×, near parity)**. The ≥5M single-run rows include the cold chain build and vary ±60% with process state (direct old-vs-new A/B at 5M/10M digits: new code equal or faster in both cold and warm) — the warm numbers in the row annotations are the steady state. The cyclic-gate split (PR #95, `BIGMATH_CYCLIC_NTT_THRESHOLD = 1280`) plus the NTT-entry retune moved the chain reciprocal towers off Karatsuba. Gap still peaks around 100k cold (4.31×, chain build not amortized). Methodology matches the original table: warm best-of-N below 100k, single cold run (chain build included) at 100k+.

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
- **Division near-balanced (ratio ≈ 2, large `b`):** PR #79 added a Newton balanced band (`b ≥ 98304, a ≥ 2b`) that dodges BZ's 5–60× non-power-of-2 blowup — 8.8× at 500k/250k, 97.7× at the 2¹⁸+1 worst case. Division is now NTT-bound for ratio ≥ 2 and ratio ≥ 3. **Residual:** ratio ∈ (1, 2) at large `b` still on BZ (~2.7× slower than Newton at ratio 1.5).
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
| `NEWTON_BALANCED_B` / ratio | 98304 / 2/1 | — | Added PR #79 — near-balanced band; below it BZ wins, above it BZ blows up 5–60× on non-pow2 sizes |
| `NEWTON_HIGH_SKEW_*` | 2048 / 8/1 | (n/a, no rows) | Keep |

**Toom-3 added to dispatch (2026-05-27, skew-gated 2026-05-30):** focused band scan (above, **Toom-3 dispatch band**) found Toom-3 wins by 4-10% in total ∈ [3 584, 4 096] and avoids a 33% NTT-length boundary regression at total 4 608. Dispatch slots Toom-3 between Karatsuba and NTT for total ∈ [2 560, 5 120), but only for non-skewed operands (`max < 2*min`). A skew-focused scan found 2:1+ shapes in that same window are faster through Karatsuba, with wins such as `2560x256` improving from 0.730 ms → 0.343 ms through `Multiply()`. `NTT_MULTIPLICATION_THRESHOLD` remains 5 120. Toom-5 still has no productive band — ties Karatsuba in 64-256 per-operand limbs and degrades sharply above 256, so stays excluded.
