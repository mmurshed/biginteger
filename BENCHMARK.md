# BigMath vs GMP benchmark results

## Single-thread comparison

On November 21, 2026, BigMath and GMP 6.3.0 each used one thread on the Apple M1 Max (10 cores, 32 GiB). The full split benchmark ran the two libraries sequentially on the same 91 operand sets as the default-build run below. GMP was timed again. The Release build used Apple Clang 21, native arm64 instructions and the default 64-bit limbs. The library and its benchmark executable were both compiled with `-DBIGMATH_USE_THREADS=0`. [Measured revision](https://github.com/mmurshed/biginteger/tree/7fd44ea0c73330918ad7481da2f0bdc8247baac4)

The full CSV includes the first call, fastest subsequent call, warm repetition count and arithmetic check status for all 111 rows. Most large sizes have one warm observation per library. All 24 product checks and 40 division checks matched GMP’s reference hashes. Each division mode was checked against the corresponding GMP operation. [README performance table](README.md#performance-numbers) · [Single-thread CSV](docs/benchmarks/2026-11-21-m1-max-single-thread.csv) · [Single-thread figure](docs/images/bigmath_vs_gmp_single_thread.png)

To repeat the single-thread run, give the driver its own build and results directory. The figure generator uses the checked-in CSVs and requires Matplotlib.

```sh
CXXFLAGS=-DBIGMATH_USE_THREADS=0 BENCH_DIR="$PWD/.benchcache/single-thread" \
  HISTORY_MD="$PWD/.benchcache/single-thread/history.md" \
  bash tests/performance/benchsuite/run_benchmark.sh full --force-gmp
python3 docs/images/make_benchmark_plots.py --ratio-only
```

## Default build with up to eight BigMath threads

The November 21, 2026 run measured PR #152 on an Apple M1 Max, with balanced multiplication through 200 million digits per operand and skewed division through 200 million ÷ 40 million digits. Both libraries produced the requested outputs on the same inputs. All 24 product checks and 40 division checks matched the GMP reference hashes. [Full results CSV](docs/benchmarks/2026-11-21-pr152-m1-max.csv) · [Benchmarked revision](https://github.com/mmurshed/biginteger/tree/c9c3e56222a17f520f7725b65c03d33232c71416)

![Four plots of BigMath time divided by GMP time for multiplication, division, decimal parsing and formatting](docs/images/bigmath_vs_gmp.png)

The arithmetic panels show warm calls. The division panel covers dividend-to-divisor digit ratios from 4:1 to 10:1. The conversion panels show both the first and warm calls. The dashed horizontal line marks equal time.

### Methodology

- The host had 10 CPU cores and 32 GiB of memory and ran macOS 27.2 (Darwin 27.2.0). Apple Clang 21 compiled the C++20 Release build with `-O3 -march=native -DNDEBUG`. BigMath used its default 64-bit limbs, CRT NTT, NEON and up to eight threads. GMP 6.3.0 ran single-threaded.
- The full split benchmark generated 91 deterministic operand sets. It ran GMP and BigMath in separate processes, sequentially, reading identical decimal operands from disk. Each of the 20 division sets was timed in both output modes, giving 111 result rows. Parsing operands, disk I/O and result hashes are outside the arithmetic timers.
- **First** means the first timed call for that operation and input in an already running benchmark process. **Warm** is the minimum of the following repetitions on the same input. Earlier cases can populate caches. Most cases above 100k digits have a single warm observation. The CSV records the actual repetition count and both timings for every row.
- Ratios are BigMath time divided by GMP time. A ratio below 1 means BigMath was faster. A zero-duration observation is shown as “<1 tick”. Its ratio is omitted. Crossovers refer to sampled operand sizes on this host.
- Earlier published division tables used BigMath quotient-only against GMP quotient and remainder. Compare each new division series with its matching output mode. [Earlier tables](https://github.com/mmurshed/biginteger/blob/c9c3e56222a17f520f7725b65c03d33232c71416/BENCHMARK.md) · [PR #129–#150 review](docs/PERFORMANCE_REVIEW_PR129.md)

To repeat the default-build run with fresh GMP timings, run the full profile from the repository root. The driver saves the paired CSV in `$BENCH_DIR/results/`.

```sh
bash tests/performance/benchsuite/run_benchmark.sh full --force-gmp
```

A separate run of `./build/bench_vs_gmp` measured both libraries in one process on different seeded operands. It completed 151 measurements. Its single 200-million-digit balanced product took 3,721 ms in BigMath and 3,057 ms in GMP (1.22×). [In-process output](docs/benchmarks/2026-11-21-pr152-in-process.txt)

### Multiplication

The equal-length series measures two independently generated operands. The skewed series holds their digit-length ratio at 10:1.

#### Equal lengths

| Operand digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 × 1,000 | 0.07642 | 0.007625 | 10.02× | 0.001291 | 0.000875 | 1.48× | 7 |
| 5,000 × 5,000 | 0.022 | 0.01871 | 1.18× | 0.01696 | 0.01254 | 1.35× | 5 |
| 10,000 × 10,000 | 0.5888 | 0.03133 | 18.79× | 0.02283 | 0.02783 | 0.82× | 5 |
| 50,000 × 50,000 | 0.4726 | 0.2881 | 1.64× | 0.1107 | 0.2723 | 0.41× | 3 |
| 100,000 × 100,000 | 0.9058 | 0.6348 | 1.43× | 0.2152 | 0.5928 | 0.36× | 3 |
| 500,000 × 500,000 | 2.55 | 4.327 | 0.59× | 1.401 | 4.276 | 0.33× | 1 |
| 1,000,000 × 1,000,000 | 4.763 | 9.495 | 0.50× | 3.247 | 9.379 | 0.35× | 1 |
| 2,000,000 × 2,000,000 | 11.11 | 21.14 | 0.53× | 8.96 | 20.69 | 0.43× | 1 |
| 5,000,000 × 5,000,000 | 19.67 | 63.53 | 0.31× | 12.23 | 62.4 | 0.20× | 1 |
| 10,000,000 × 10,000,000 | 44.85 | 209.7 | 0.21× | 29.75 | 205.7 | 0.14× | 1 |
| 20,000,000 × 20,000,000 | 96.73 | 276.4 | 0.35× | 67.6 | 270.1 | 0.25× | 1 |
| 50,000,000 × 50,000,000 | 612.1 | 680.7 | 0.90× | 451.2 | 657.1 | 0.69× | 1 |
| 100,000,000 × 100,000,000 | 1,322 | 1,415 | 0.93× | 992 | 1,391 | 0.71× | 1 |
| 200,000,000 × 200,000,000 | 4,518 | 2,909 | 1.55× | 2,325 | 2,868 | 0.81× | 1 |

#### Skewed lengths

| Operand digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 100,000 × 10,000 | 0.2257 | 0.3299 | 0.68× | 0.1202 | 0.3011 | 0.40× | 3 |
| 500,000 × 50,000 | 0.6862 | 2.112 | 0.32× | 0.4202 | 2.043 | 0.21× | 1 |
| 1,000,000 × 100,000 | 1.303 | 4.488 | 0.29× | 1.279 | 4.437 | 0.29× | 1 |
| 2,000,000 × 200,000 | 3.107 | 9.346 | 0.33× | 3.191 | 9.384 | 0.34× | 1 |
| 5,000,000 × 500,000 | 10.11 | 31.88 | 0.32× | 7.651 | 31.79 | 0.24× | 1 |
| 10,000,000 × 1,000,000 | 22.7 | 70.75 | 0.32× | 19.28 | 70.55 | 0.27× | 1 |
| 20,000,000 × 2,000,000 | 70.97 | 157.8 | 0.45× | 69.07 | 157.6 | 0.44× | 1 |
| 50,000,000 × 5,000,000 | 204.7 | 680.4 | 0.30× | 203.3 | 670.3 | 0.30× | 1 |
| 100,000,000 × 10,000,000 | 471.3 | 1,022 | 0.46× | 450.2 | 1,009 | 0.45× | 1 |
| 200,000,000 × 20,000,000 | 977.2 | 1,915 | 0.51× | 1,013 | 1,915 | 0.53× | 1 |

### Division

Quotient-only uses BigMath `/` and GMP `mpz_tdiv_q`. Quotient and remainder uses BigMath `DivideAndRemainder` and GMP `mpz_tdiv_qr`. Both modes use the same dividend and divisor. Equal-digit operands can produce a quotient of zero or one, so those results depend on the sampled values.

#### Skewed dividends

| Dividend ÷ divisor digits | BigMath Q ms | GMP Q ms | Q ratio | BigMath QR ms | GMP QR ms | QR ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 40,000 ÷ 10,000 | 0.2941 | 0.1982 | 1.48× | 0.2863 | 0.2405 | 1.19× | 3 |
| 100,000 ÷ 10,000 | 0.8523 | 0.4327 | 1.97× | 0.8488 | 0.4658 | 1.82× | 3 |
| 200,000 ÷ 50,000 | 1.013 | 1.52 | 0.67× | 0.9871 | 1.642 | 0.60× | 1 |
| 500,000 ÷ 100,000 | 2.287 | 4.3 | 0.53× | 2.247 | 4.5 | 0.50× | 1 |
| 1,000,000 ÷ 200,000 | 4.794 | 9.531 | 0.50× | 4.769 | 10.14 | 0.47× | 1 |
| 2,000,000 ÷ 500,000 | 8.051 | 22.69 | 0.35× | 7.947 | 24.24 | 0.33× | 1 |
| 5,000,000 ÷ 1,000,000 | 23.99 | 64.15 | 0.37× | 25.14 | 68.85 | 0.37× | 1 |
| 10,000,000 ÷ 2,000,000 | 55.32 | 142.9 | 0.39× | 52.26 | 150.6 | 0.35× | 1 |
| 20,000,000 ÷ 4,000,000 | 104.7 | 321.1 | 0.33× | 107.2 | 341.4 | 0.31× | 1 |
| 50,000,000 ÷ 10,000,000 | 245.6 | 1,232 | 0.20× | 230.6 | 1,306 | 0.18× | 1 |
| 100,000,000 ÷ 20,000,000 | 532.9 | 2,316 | 0.23× | 538.6 | 2,515 | 0.21× | 1 |
| 200,000,000 ÷ 40,000,000 | 1,359 | 4,183 | 0.32× | 1,380 | 4,443 | 0.31× | 1 |

At 200,000 ÷ 50,000 digits, the in-process run measured warm ratios of 0.92× for quotient-only and 1.04× for quotient and remainder. Each suite took one warm measurement per library at that size. Both suites favored BigMath in both modes from the next sampled shape, 500,000 ÷ 100,000, through 200 million ÷ 40 million digits. [In-process output](docs/benchmarks/2026-11-21-pr152-in-process.txt)

#### Equal digit counts

| Dividend ÷ divisor digits | BigMath Q ms | GMP Q ms | Q ratio | BigMath QR ms | GMP QR ms | QR ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 ÷ 1,000 | 0.000083 | <1 tick | — | 0.000125 | 0.000083 | 1.51× | 7 |
| 5,000 ÷ 5,000 | 0.000083 | <1 tick | — | 0.000417 | 0.000208 | 2.00× | 5 |
| 10,000 ÷ 10,000 | <1 tick | <1 tick | — | 0.000208 | 0.000042 | 4.95× | 5 |
| 50,000 ÷ 50,000 | 0.000083 | 0.000041 | 2.02× | 0.003958 | 0.00175 | 2.26× | 3 |
| 100,000 ÷ 100,000 | 0.000084 | 0.000125 | 0.67× | 0.008042 | 0.0035 | 2.30× | 3 |
| 500,000 ÷ 500,000 | 0.000083 | 0.000167 | 0.50× | 0.009375 | 0.004292 | 2.18× | 1 |
| 1,000,000 ÷ 1,000,000 | 0.000125 | 0.000458 | 0.27× | 0.01604 | 0.008625 | 1.86× | 1 |
| 5,000,000 ÷ 5,000,000 | 0.000084 | 0.000208 | 0.40× | 0.08446 | 0.043 | 1.96× | 1 |

### Decimal parse (string → BigInteger)

| Digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.007208 | 0.003458 | 2.08× | 0.002375 | 0.001583 | 1.50× | 7 |
| 10,000 | 0.07496 | 0.05517 | 1.36× | 0.06383 | 0.03917 | 1.63× | 5 |
| 50,000 | 0.5564 | 0.4261 | 1.31× | 0.5193 | 0.3922 | 1.32× | 3 |
| 100,000 | 0.9724 | 1.149 | 0.85× | 0.5443 | 1.059 | 0.51× | 3 |
| 500,000 | 3.578 | 9.147 | 0.39× | 3.464 | 9.109 | 0.38× | 1 |
| 1,000,000 | 6.578 | 20.61 | 0.32× | 6.389 | 20.61 | 0.31× | 1 |
| 2,000,000 | 15.46 | 47.61 | 0.32× | 15.04 | 47.29 | 0.32× | 1 |
| 5,000,000 | 36.12 | 151.5 | 0.24× | 37.05 | 147.6 | 0.25× | 1 |
| 10,000,000 | 77.92 | 344.1 | 0.23× | 82.61 | 341.9 | 0.24× | 1 |
| 20,000,000 | 180.6 | 792.2 | 0.23× | 180 | 789.7 | 0.23× | 1 |
| 50,000,000 | 707.4 | 2,633 | 0.27× | 718.5 | 2,647 | 0.27× | 1 |

### ToString (BigInteger → string)

| Digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.007833 | 0.006208 | 1.26× | 0.006292 | 0.003375 | 1.86× | 7 |
| 10,000 | 0.3796 | 0.09775 | 3.88× | 0.1792 | 0.0765 | 2.34× | 5 |
| 50,000 | 1.818 | 0.9104 | 2.00× | 1.348 | 0.843 | 1.60× | 3 |
| 100,000 | 2.228 | 2.357 | 0.95× | 1.462 | 2.392 | 0.61× | 1 |
| 200,000 | 4.335 | 6.117 | 0.71× | 2.732 | 6.039 | 0.45× | 1 |
| 500,000 | 11 | 20.44 | 0.54× | 6.849 | 20.37 | 0.34× | 1 |
| 1,000,000 | 22.06 | 49.42 | 0.45× | 14.94 | 49.8 | 0.30× | 1 |
| 2,000,000 | 44.94 | 119.5 | 0.38× | 34.17 | 118.8 | 0.29× | 1 |
| 5,000,000 | 165 | 381.9 | 0.43× | 123.1 | 382.5 | 0.32× | 1 |
| 10,000,000 | 359 | 908.3 | 0.40× | 276.4 | 905.6 | 0.31× | 1 |
| 20,000,000 | 796.3 | 2,128 | 0.37× | 613.8 | 2,127 | 0.29× | 1 |

### BigDecimal

BigMath uses decimal fixed-point arithmetic and HALF_EVEN rounding for division. The comparison uses GMP binary floating-point (`mpf_t`), with precision derived from the digit count and requested scale. The split runner times these operations without an equality check. A size of `1,000.100` denotes 1,000 integer digits and 100 fractional digits.

| Operation | Size | BigMath warm ms | GMP warm ms | Ratio | Warm reps |
|---|---:|---:|---:|---:|---:|
| add | 100.10 digits | 0.000125 | <1 tick | — | 7 |
| add | 1000.100 digits | 0.000166 | 0.000041 | 4.05× | 7 |
| add | 5000.500 digits | 0.000375 | 0.000125 | 3.00× | 5 |
| add | 20000.2000 digits | 0.00125 | 0.0005 | 2.50× | 3 |
| mul | 100.10 digits | 0.000041 | 0.000041 | 1.00× | 7 |
| mul | 1000.100 digits | 0.0015 | 0.001083 | 1.39× | 7 |
| mul | 5000.500 digits | 0.0205 | 0.01275 | 1.61× | 5 |
| mul | 20000.2000 digits | 0.05021 | 0.09546 | 0.53× | 3 |
| div | 100.10 / 10 dp | 0.000458 | 0.000125 | 3.66× | 7 |
| div | 1000.100 / 100 dp | 0.001 | 0.001875 | 0.53× | 7 |
| div | 5000.500 / 500 dp | 0.00975 | 0.02554 | 0.38× | 5 |
| div | 20000.2000 / 2000 dp | 0.06404 | 0.1923 | 0.33× | 3 |
| div | 2000.200 / 0 dp | 0.000542 | 0.005209 | 0.10× | 5 |
| div | 2000.200 / 10 dp | 0.00075 | 0.005208 | 0.14× | 5 |
| div | 2000.200 / 100 dp | 0.001458 | 0.005459 | 0.27× | 5 |
| div | 2000.200 / 1000 dp | 0.007375 | 0.00925 | 0.80× | 5 |
| div | 2000.200 / 5000 dp | 0.03083 | 0.03358 | 0.92× | 5 |
| parse | 100 digits | 0.000375 | 0.0005 | 0.75× | 7 |
| parse | 1000 digits | 0.004666 | 0.004583 | 1.02× | 7 |
| parse | 10000 digits | 0.08542 | 0.1109 | 0.77× | 5 |
| parse | 50000 digits | 0.8206 | 1.002 | 0.82× | 3 |
| tostr | 100 digits | 0.000291 | 0.000458 | 0.64× | 7 |
| tostr | 1000 digits | 0.006583 | 0.004666 | 1.41× | 7 |
| tostr | 10000 digits | 0.1809 | 0.1014 | 1.78× | 5 |
| tostr | 50000 digits | 1.346 | 1.086 | 1.24× | 3 |

## Internal algorithm crossovers

The November 21 run of `multperf_simple 5` used equal-length random operands. Times are the mean of five calls, and speedups are relative to classical multiplication. All products matched at all four tested sizes.

| digits | Classical ms | Karatsuba ms | Karatsuba× | NTT ms | NTT× |
|---:|---:|---:|---:|---:|---:|
| 10,000 | 0.196 | 0.070 | 2.8× | 0.107 | 1.8× |
| 100,000 | 15.479 | 2.111 | 7.3× | 0.468 | 33.1× |
| 500,000 | 379.127 | 27.560 | 13.8× | 1.608 | 235.7× |

The default 64-bit NTT path requires at least 640 total limbs, 64 limbs in
the shorter operand and an area of at least 96,000 limb pairs. With a parked
pool the area floor is 262,144. Balanced entry is 320 limbs per operand
when warm and 512 when cold.

In `divperf_simple 5`, Burnikel–Ziegler and Knuth division returned matching quotients and remainders at all four tested shapes. At 4,096 ÷ 2,048 limbs they averaged 0.846 ms and 3.028 ms respectively. This harness draws 32 random bits per limb.

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
