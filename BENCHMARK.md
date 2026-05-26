# BigMath vs GMP — full benchmark suite

Measured 2026-05-26. Single-run snapshot capturing the state of `master` after the 2026-05 optimization pass (LIMB_64 default, multi-prime CRT NTT default, multithreaded NTT default, M-G div2by1 in `ClassicDivision`, M-G 3/2 qhat in `FastDivision` Base2_64, BZ Knuth-normalize fix).

For per-subsystem deep dives — algorithms, dispatch, optimization history, rejected approaches — see:

- [docs/MULTIPLICATION.md](docs/MULTIPLICATION.md)
- [docs/DIVISION.md](docs/DIVISION.md)
- [docs/STRING_CONVERSION.md](docs/STRING_CONVERSION.md)
- [docs/BASE.md](docs/BASE.md)

This file is the **raw numbers**, all four subsystems in one place, kept for headline-level reference.

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

- `tests/performance/bench_vs_gmp.cpp` — canonical bench, covers sizes ≤1M digits, all four operations.
- Ad-hoc extension covering 1M-10M digit operands for multiplication, division, and parse; 200k-2M for ToString. Source not committed; reproducible from the same template against operand size lists.

---

## Multiplication

Balanced (`a.size() == b.size()`):

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 1 000 × 1 000 | 0.002 | 0.001 | 2.13× |
| 5 000 × 5 000 | 0.031 | 0.012 | 2.60× |
| 10 000 × 10 000 | 0.093 | 0.038 | 2.44× |
| 50 000 × 50 000 | 0.610 | 0.329 | 1.86× |
| 100 000 × 100 000 | 1.066 | 0.633 | 1.69× |
| 500 000 × 500 000 | 6.043 | 4.193 | 1.44× |
| 1 000 000 × 1 000 000 | 9.731 | 8.625 | 1.13× |
| 2 000 000 × 2 000 000 | 27.249 | 20.641 | 1.32× |
| **5 000 000 × 5 000 000** | **57.976** | **64.950** | **0.89×** ← BigMath faster |
| **10 000 000 × 10 000 000** | **140.222** | **214.225** | **0.65×** ← BigMath faster |

Skewed (`a.size() >> b.size()`):

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 100 000 × 10 000 | 0.537 | 0.301 | 1.78× |
| 500 000 × 50 000 | 2.256 | 2.107 | 1.07× |
| 1 000 000 × 100 000 | 5.030 | 4.473 | 1.12× |
| **2 000 000 × 200 000** | **9.478** | **9.675** | **0.98×** ← parity |
| 5 000 000 × 500 000 | 43.380 | 30.462 | 1.42× |
| 10 000 000 × 1 000 000 | 115.111 | 70.971 | 1.62× |

**Observations:**

- **BigMath beats GMP on balanced multiplication at ≥5M digits.** NTT's asymptotic edge overtakes GMP's hand-tuned Karatsuba/Toom assembly. At 10M balanced, BigMath is **35% faster** (140ms vs 214ms).
- Below 1M, GMP's hand-tuned basecase keeps a 1.4-2.6× lead.
- Skewed mults stay within a 1.0-1.6× band across the size range. The 2M×200k case is at parity; larger skewed sizes drag the smaller operand into BigMath's NTT-bound regime asymmetrically while GMP keeps the small operand in its well-tuned mid-band.

---

## Division

Balanced (`a.size() == b.size()`) — quotient is 1-2 limbs, both libraries short-circuit. Numbers are sub-microsecond noise; reported for completeness.

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 1 000 × 1 000 | 0.001 | 0.000 | 11.90× |
| 5 000 × 5 000 | 0.001 | 0.000 | 8.73× |
| 10 000 × 10 000 | 0.003 | 0.000 | 7.00× |
| 50 000 × 50 000 | 0.000 | 0.000 | 0.25× |
| 100 000 × 100 000 | 0.000 | 0.001 | 0.17× |
| 500 000 × 500 000 | 0.003 | 0.005 | 0.55× |
| 1 000 000 × 1 000 000 | 0.217 | 0.033 | 6.66× |
| 5 000 000 × 5 000 000 | 1.184 | 0.166 | 7.15× |

Skewed (`a.size() >> b.size()`) — Newton/BZ band, real algorithmic work:

| size | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 40 000 × 10 000 | 1.023 | 0.225 | 4.55× |
| 100 000 × 10 000 | 3.068 | 0.467 | 6.57× |
| 200 000 × 50 000 | 9.077 | 1.706 | 5.32× |
| 500 000 × 100 000 | 18.642 | 4.601 | 4.05× |
| 1 000 000 × 200 000 | 36.372 | 9.997 | 3.64× |
| 2 000 000 × 500 000 | 88.980 | 23.796 | 3.74× |
| 5 000 000 × 1 000 000 | 221.005 | 66.484 | 3.32× |
| 10 000 000 × 2 000 000 | 460.387 | 151.167 | 3.05× |

**Observations:**

- **Skewed division ratio narrows from ~6× at 100k-divisor band to ~3× at 10M/2M.** Newton inherits the multiplication win that overtakes GMP at 5M.
- Residual ~3× gap is the structural cost of Newton's chunked iteration vs GMP's `mpn_dcpi1_div_q` (a single recursive divide rather than reciprocal-then-chunk).
- Balanced cases route through trivial short-circuits and aren't algorithmically meaningful at this size profile.

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
| 1 000 | 0.003 | 0.002 | 1.80× |
| 10 000 | 0.114 | 0.037 | 3.06× |
| 50 000 | 1.395 | 0.387 | 3.60× |
| 100 000 | 3.281 | 1.045 | 3.14× |
| 500 000 | 21.903 | 9.138 | 2.40× |
| 1 000 000 | 48.160 | 20.412 | 2.36× |
| 2 000 000 | 106.767 | 48.054 | 2.22× |
| 5 000 000 | 277.142 | 147.445 | 1.88× |
| 10 000 000 | 598.642 | 346.629 | 1.73× |

**Observation:** ratio narrows monotonically with size (3.1× at 100k → 1.7× at 10M). The D&C parser's combine step is NTT-bound; as BigMath's NTT overtakes GMP at 5M, the parser inherits the win.

---

## ToString (BigInteger → string)

| size (digits) | BigMath ms | GMP ms | BM/GMP |
|---|---:|---:|---:|
| 1 000 | 0.009 | 0.004 | 2.33× |
| 10 000 | 0.875 | 0.083 | 10.61× |
| 50 000 | 9.436 | 0.845 | 11.17× |
| 100 000 | 21.516 | 2.471 | 8.71× |
| 200 000 | 42.937 | 6.286 | 6.83× |
| 500 000 | 118.444 | 21.054 | 5.63× |
| 1 000 000 | 247.620 | 50.602 | 4.89× |
| 2 000 000 | 532.934 | 119.440 | 4.46× |

**Observation:** narrowest gap at 1k (the linear leaf, where GM div2by1 already runs); peaks at 50k (D&C overhead + Newton recip setup not yet amortized); narrows again from 100k onward as D&C asymptotic + NTT inheriting from multiplication's overtake compound.

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

## Headline summary

- **Multiplication ≥5M digits:** BigMath faster than GMP (0.89× at 5M, 0.65× at 10M balanced).
- **Multiplication skewed 2M×200k:** parity (0.98×).
- **Division skewed:** ~3-6× behind GMP, narrowing slowly with size (NTT-bound — multiplication win flows through Newton's reciprocal iteration).
- **ToString:** 4.5× at 2M in the GMP comparison table, with focused warm BigMath-only results now showing a 2-4× improvement from cached divider chains and leaf formatting changes. The remaining gap is still concentrated in Newton's per-chunk divmod and NTT-backed multiplication inside division.
- **Parse:** 1.7× at 10M, monotonically narrowing.

For optimizations considered and rejected with measurement evidence, see the **Explored but rejected** sections of each subsystem doc. The 2026-05 optimization stack (LIMB_64 + CRT NTT + threading + M-G reciprocals + BZ Knuth fix) closed the GMP gap by 3-5× across every band; the remaining gap concentrates in the basecase NTT butterfly inner loop where GMP's hand-tuned ARM64 assembly maintains an edge.
