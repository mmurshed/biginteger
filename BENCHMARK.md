# BigMath vs GMP benchmark results

## Single-thread comparison

On November 21, 2026, BigMath and GMP 6.3.0 each used one thread on the Apple M1 Max (10 cores, 32 GiB). The full split benchmark ran the two libraries sequentially on the same 91 operand sets as the default-build run below, and both comparisons use the same GMP pass. The Release build used Apple Clang 21, native arm64 instructions and the default 64-bit limbs. The library and its benchmark executable were both compiled with `-DBIGMATH_USE_THREADS=0`. The measured revision is `main` after #162, which divides the small levels of ToString's power-of-ten split with the regular division code instead of a cached Newton reciprocal (from 13,000 digits per level the reciprocal is still used). Against the previous run, warm ToString took 7–36% less time from 10,000 digits up (10,000 digits: 0.173 → 0.110 ms, 2.23× → 1.42× GMP's time) and now beats GMP at 1 million and 2 million digits; first calls fell from 3.11× to 1.61× GMP's time at 10,000 digits. BigDecimal formatting, which goes through the same code, took 36% less at 10,000 digits and 23% less at 50,000. 62 of the other 65 rows of 10 µs or longer moved by under 5%. [Measured revision](https://github.com/mmurshed/biginteger/tree/a3c6aaa) · [Previous single-thread run](docs/benchmarks/2026-11-21-ntt-trunc-m1-max-single-thread.csv)

The full CSV includes the first call, fastest subsequent call, warm repetition count and arithmetic check status for all 111 rows. Most large sizes have one warm observation per library. All 24 product checks and 40 division checks matched GMP’s reference hashes. Each division mode was checked against the corresponding GMP operation. [README performance table](README.md#performance-numbers) · [Single-thread CSV](docs/benchmarks/2026-11-21-tostring-m1-max-single-thread.csv) · [Single-thread figure](docs/images/bigmath_vs_gmp_single_thread.png)

To repeat the single-thread run, give the driver its own build and results directory. The figure generator uses the checked-in CSVs and requires Matplotlib.

```sh
CXXFLAGS=-DBIGMATH_USE_THREADS=0 BENCH_DIR="$PWD/.benchcache/single-thread" \
  HISTORY_MD="$PWD/.benchcache/single-thread/history.md" \
  bash tests/performance/benchsuite/run_benchmark.sh full --force-gmp
python3 docs/images/make_benchmark_plots.py --ratio-only
```

## Default build with up to eight BigMath threads

The November 21, 2026 run measured `main` after the ToString change (#162) on an Apple M1 Max, with balanced multiplication through 200 million digits per operand and skewed division through 200 million ÷ 40 million digits. Both libraries produced the requested outputs on the same inputs. All 24 product checks and 40 division checks matched the GMP reference hashes. Against the previous run, warm ToString took 5–37% less time from 10,000 to 10 million digits (10,000 digits: 2.28× → 1.45× GMP's time) and 4.5% less at 20 million, and the first call at 10,000 digits fell from 5.01× to 1.55× GMP's time. [Full results CSV](docs/benchmarks/2026-11-21-tostring-m1-max.csv) · [Measured revision](https://github.com/mmurshed/biginteger/tree/a3c6aaa) · [Previous run (large-NTT pass, #159–#160)](docs/benchmarks/2026-11-21-ntt-trunc-m1-max.csv) · [Run before that (Toom-3 tier, #157)](docs/benchmarks/2026-11-21-toom-m1-max.csv)

**Run-to-run variation.** GMP was not re-run for these CSVs; both comparisons reuse the GMP pass first published with the large-NTT run, which a same-day single-thread re-measurement matched within 3.6% on its 55 rows of 1 ms or longer (median +0.4%), except the 500,000-digit product (+12%). The published CSVs come from the second of two runs of both BigMath passes. The first ran straight after the build, and its small rows were slow again (single-threaded products of up to 100,000 digits took 1.8–2.5× as long as in the second run), as in the previous two benchmark sessions. Rows with a single warm observation still vary between runs, and the multithreaded ones most: 13 of the 65 non-ToString 8-thread rows of 10 µs or longer moved by 5% or more against the previous run, from −16% (2,000,000 × 200,000-digit product) to +24% (the 200-million-digit product, 1.72 → 2.13 s; that row has measured 1.72–2.32 s across the recent runs) and +128% (5,000,000-digit equal-length quotient and remainder, a 0.23 ms row). Single-thread rows moved much less; the three outside 5% were the 10,000,000 × 1,000,000-digit product (+6.0%) and two equal-length quotient-and-remainder rows of under 0.1 ms: 1,000,000 digits (+6.7%) and 5,000,000 digits (0.225 → 0.094 ms, back in the 0.071–0.094 ms range of the runs before the last).

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

A separate run of `./build/bench_vs_gmp` on PR #152 measured both libraries in one process on different seeded operands. It completed 151 measurements. Its single 200-million-digit balanced product took 3,721 ms in BigMath and 3,057 ms in GMP (1.22×). [In-process output](docs/benchmarks/2026-11-21-pr152-in-process.txt)

### Multiplication

The equal-length series measures two independently generated operands. The skewed series holds their digit-length ratio at 10:1.

#### Equal lengths

| Operand digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 × 1,000 | 0.01 | 0.1118 | 0.09× | 0.001291 | 0.000959 | 1.35× | 7 |
| 5,000 × 5,000 | 0.02221 | 0.03833 | 0.58× | 0.01642 | 0.01221 | 1.34× | 5 |
| 10,000 × 10,000 | 0.2412 | 0.03175 | 7.60× | 0.02521 | 0.02804 | 0.90× | 5 |
| 50,000 × 50,000 | 0.469 | 0.2815 | 1.67× | 0.1092 | 0.2734 | 0.40× | 3 |
| 100,000 × 100,000 | 0.8777 | 0.6323 | 1.39× | 0.2004 | 0.588 | 0.34× | 3 |
| 500,000 × 500,000 | 1.911 | 4.277 | 0.45× | 1.153 | 4.234 | 0.27× | 1 |
| 1,000,000 × 1,000,000 | 4.943 | 8.949 | 0.55× | 2.971 | 8.847 | 0.34× | 1 |
| 2,000,000 × 2,000,000 | 11.09 | 20.31 | 0.55× | 7.316 | 20.33 | 0.36× | 1 |
| 5,000,000 × 5,000,000 | 17.81 | 66.13 | 0.27× | 11.12 | 62.79 | 0.18× | 1 |
| 10,000,000 × 10,000,000 | 41.01 | 202.8 | 0.20× | 27.26 | 199.9 | 0.14× | 1 |
| 20,000,000 × 20,000,000 | 88.96 | 266.4 | 0.33× | 63.74 | 265.3 | 0.24× | 1 |
| 50,000,000 × 50,000,000 | 459.3 | 660.8 | 0.70× | 339.2 | 650.7 | 0.52× | 1 |
| 100,000,000 × 100,000,000 | 1,000 | 1,393 | 0.72× | 782.6 | 1,376 | 0.57× | 1 |
| 200,000,000 × 200,000,000 | 3,053 | 2,862 | 1.07× | 2,128 | 2,828 | 0.75× | 1 |

#### Skewed lengths

| Operand digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 100,000 × 10,000 | 0.1537 | 0.3218 | 0.48× | 0.126 | 0.3015 | 0.42× | 3 |
| 500,000 × 50,000 | 0.4763 | 2.068 | 0.23× | 0.4736 | 2.039 | 0.23× | 1 |
| 1,000,000 × 100,000 | 1.29 | 4.455 | 0.29× | 1.191 | 4.429 | 0.27× | 1 |
| 2,000,000 × 200,000 | 3.011 | 9.312 | 0.32× | 2.691 | 9.244 | 0.29× | 1 |
| 5,000,000 × 500,000 | 9.237 | 30.85 | 0.30× | 7.999 | 30.84 | 0.26× | 1 |
| 10,000,000 × 1,000,000 | 19.72 | 70.04 | 0.28× | 20.28 | 69.68 | 0.29× | 1 |
| 20,000,000 × 2,000,000 | 64.85 | 158.6 | 0.41× | 63.67 | 158.2 | 0.40× | 1 |
| 50,000,000 × 5,000,000 | 163.3 | 662.8 | 0.25× | 158.1 | 660.4 | 0.24× | 1 |
| 100,000,000 × 10,000,000 | 380.1 | 996.3 | 0.38× | 393.2 | 994.4 | 0.40× | 1 |
| 200,000,000 × 20,000,000 | 853.1 | 1,893 | 0.45× | 829.5 | 1,890 | 0.44× | 1 |

### Division

Quotient-only uses BigMath `/` and GMP `mpz_tdiv_q`. Quotient and remainder uses BigMath `DivideAndRemainder` and GMP `mpz_tdiv_qr`. Both modes use the same dividend and divisor. Equal-digit operands can produce a quotient of zero or one, so those results depend on the sampled values.

#### Skewed dividends

| Dividend ÷ divisor digits | BigMath Q ms | GMP Q ms | Q ratio | BigMath QR ms | GMP QR ms | QR ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 40,000 ÷ 10,000 | 0.2854 | 0.2001 | 1.43× | 0.2855 | 0.2249 | 1.27× | 3 |
| 100,000 ÷ 10,000 | 0.8476 | 0.4268 | 1.99× | 0.8468 | 0.4409 | 1.92× | 3 |
| 200,000 ÷ 50,000 | 0.9885 | 1.507 | 0.66× | 0.9855 | 1.626 | 0.61× | 1 |
| 500,000 ÷ 100,000 | 2.328 | 4.207 | 0.55× | 2.314 | 4.491 | 0.52× | 1 |
| 1,000,000 ÷ 200,000 | 4.246 | 9.476 | 0.45× | 4.33 | 10.08 | 0.43× | 1 |
| 2,000,000 ÷ 500,000 | 7.996 | 22.23 | 0.36× | 7.87 | 23.92 | 0.33× | 1 |
| 5,000,000 ÷ 1,000,000 | 27.05 | 63.21 | 0.43× | 24.25 | 67.66 | 0.36× | 1 |
| 10,000,000 ÷ 2,000,000 | 54.75 | 144.5 | 0.38× | 54.76 | 154.2 | 0.36× | 1 |
| 20,000,000 ÷ 4,000,000 | 105.1 | 319.7 | 0.33× | 102.3 | 342.2 | 0.30× | 1 |
| 50,000,000 ÷ 10,000,000 | 223.6 | 1,213 | 0.18× | 226.6 | 1,276 | 0.18× | 1 |
| 100,000,000 ÷ 20,000,000 | 510.7 | 2,288 | 0.22× | 506.9 | 2,490 | 0.20× | 1 |
| 200,000,000 ÷ 40,000,000 | 1,180 | 4,137 | 0.29× | 1,171 | 4,452 | 0.26× | 1 |

At 200,000 ÷ 50,000 digits, the in-process run measured warm ratios of 0.92× for quotient-only and 1.04× for quotient and remainder. Each suite took one warm measurement per library at that size. Both suites favored BigMath in both modes from the next sampled shape, 500,000 ÷ 100,000, through 200 million ÷ 40 million digits. [In-process output](docs/benchmarks/2026-11-21-pr152-in-process.txt)

#### Equal digit counts

| Dividend ÷ divisor digits | BigMath Q ms | GMP Q ms | Q ratio | BigMath QR ms | GMP QR ms | QR ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 ÷ 1,000 | 0.000083 | <1 tick | — | 0.000166 | 0.000042 | 3.95× | 7 |
| 5,000 ÷ 5,000 | 0.000083 | <1 tick | — | 0.000416 | 0.000167 | 2.49× | 5 |
| 10,000 ÷ 10,000 | <1 tick | <1 tick | — | 0.000208 | 0.000042 | 4.95× | 5 |
| 50,000 ÷ 50,000 | 0.000083 | <1 tick | — | 0.003917 | 0.001666 | 2.35× | 3 |
| 100,000 ÷ 100,000 | 0.000125 | 0.000125 | 1.00× | 0.008375 | 0.003583 | 2.34× | 3 |
| 500,000 ÷ 500,000 | 0.000084 | 0.000167 | 0.50× | 0.009125 | 0.004292 | 2.13× | 1 |
| 1,000,000 ÷ 1,000,000 | 0.000042 | 0.000167 | 0.25× | 0.01754 | 0.008583 | 2.04× | 1 |
| 5,000,000 ÷ 5,000,000 | 0.000125 | 0.00025 | 0.50× | 0.2255 | 0.04288 | 5.26× | 1 |

### Decimal parse (string → BigInteger)

| Digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.008458 | 0.004542 | 1.86× | 0.002416 | 0.001541 | 1.57× | 7 |
| 10,000 | 0.07979 | 0.05375 | 1.48× | 0.06383 | 0.03771 | 1.69× | 5 |
| 50,000 | 0.5838 | 0.4271 | 1.37× | 0.523 | 0.3946 | 1.33× | 3 |
| 100,000 | 0.6192 | 1.078 | 0.57× | 0.5464 | 1.061 | 0.52× | 3 |
| 500,000 | 2.994 | 9.019 | 0.33× | 2.957 | 8.93 | 0.33× | 1 |
| 1,000,000 | 6.606 | 20.99 | 0.31× | 6.607 | 20.78 | 0.32× | 1 |
| 2,000,000 | 14.98 | 49.72 | 0.30× | 15.12 | 47.15 | 0.32× | 1 |
| 5,000,000 | 36.14 | 146.3 | 0.25× | 36.31 | 146 | 0.25× | 1 |
| 10,000,000 | 85.04 | 341.6 | 0.25× | 78.92 | 343.4 | 0.23× | 1 |
| 20,000,000 | 175.3 | 805.9 | 0.22× | 177.8 | 802.5 | 0.22× | 1 |
| 50,000,000 | 663.4 | 2,627 | 0.25× | 666.2 | 2,650 | 0.25× | 1 |

### ToString (BigInteger → string)

| Digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.007875 | 0.0065 | 1.21× | 0.00625 | 0.003375 | 1.85× | 7 |
| 10,000 | 0.1442 | 0.09275 | 1.55× | 0.1118 | 0.07733 | 1.45× | 5 |
| 50,000 | 1.458 | 0.888 | 1.64× | 1.001 | 0.8432 | 1.19× | 3 |
| 100,000 | 1.878 | 2.285 | 0.82× | 1.146 | 2.267 | 0.51× | 1 |
| 200,000 | 3.812 | 6.026 | 0.63× | 2.357 | 5.961 | 0.40× | 1 |
| 500,000 | 8.934 | 20.1 | 0.44× | 5.649 | 20.08 | 0.28× | 1 |
| 1,000,000 | 18 | 49.15 | 0.37× | 12.42 | 49.38 | 0.25× | 1 |
| 2,000,000 | 40.79 | 117.6 | 0.35× | 28.18 | 117.4 | 0.24× | 1 |
| 5,000,000 | 151.3 | 376.8 | 0.40× | 113.2 | 376.6 | 0.30× | 1 |
| 10,000,000 | 336.1 | 893.7 | 0.38× | 250.4 | 901.9 | 0.28× | 1 |
| 20,000,000 | 771.8 | 2,092 | 0.37× | 580.8 | 2,094 | 0.28× | 1 |

### BigDecimal

BigMath uses decimal fixed-point arithmetic and HALF_EVEN rounding for division. The comparison uses GMP binary floating-point (`mpf_t`), with precision derived from the digit count and requested scale. The split runner times these operations without an equality check. A size of `1,000.100` denotes 1,000 integer digits and 100 fractional digits.

| Operation | Size | BigMath warm ms | GMP warm ms | Ratio | Warm reps |
|---|---:|---:|---:|---:|---:|
| add | 100.10 digits | 0.000083 | <1 tick | — | 7 |
| mul | 100.10 digits | 0.000042 | 0.000041 | 1.02× | 7 |
| add | 1000.100 digits | 0.000125 | <1 tick | — | 7 |
| mul | 1000.100 digits | 0.0015 | 0.001042 | 1.44× | 7 |
| add | 5000.500 digits | 0.000291 | 0.000125 | 2.33× | 5 |
| mul | 5000.500 digits | 0.01887 | 0.01321 | 1.43× | 5 |
| add | 20000.2000 digits | 0.001125 | 0.0005 | 2.25× | 3 |
| mul | 20000.2000 digits | 0.046 | 0.09292 | 0.50× | 3 |
| div | 100.10 / 10 dp | 0.00025 | 0.000083 | 3.01× | 7 |
| div | 1000.100 / 100 dp | 0.000875 | 0.001833 | 0.48× | 7 |
| div | 5000.500 / 500 dp | 0.009583 | 0.02558 | 0.37× | 5 |
| div | 20000.2000 / 2000 dp | 0.06225 | 0.1921 | 0.32× | 3 |
| div | 2000.200 / 0 dp | 0.0005 | 0.005 | 0.10× | 5 |
| div | 2000.200 / 10 dp | 0.000708 | 0.005083 | 0.14× | 5 |
| div | 2000.200 / 100 dp | 0.001334 | 0.005458 | 0.24× | 5 |
| div | 2000.200 / 1000 dp | 0.007375 | 0.009333 | 0.79× | 5 |
| div | 2000.200 / 5000 dp | 0.02996 | 0.03367 | 0.89× | 5 |
| parse | 100 digits | 0.000375 | 0.000459 | 0.82× | 7 |
| parse | 1000 digits | 0.004583 | 0.004542 | 1.01× | 7 |
| parse | 10000 digits | 0.08417 | 0.1045 | 0.81× | 5 |
| parse | 50000 digits | 0.6158 | 0.9947 | 0.62× | 3 |
| tostr | 100 digits | 0.000292 | 0.0005 | 0.58× | 7 |
| tostr | 1000 digits | 0.006625 | 0.004583 | 1.45× | 7 |
| tostr | 10000 digits | 0.1124 | 0.1012 | 1.11× | 5 |
| tostr | 50000 digits | 1.012 | 1.084 | 0.93× | 3 |

## Small operands

`tests/performance/small_ops_vs_gmp.cpp` times the BigInteger operators against GMP on operands of 1–16 limbs (up to about 300 decimal digits), where per-call overhead and allocation decide the cost. It ran single-threaded on the same machine against two builds: **before** is `main` at 3dafa62, where every value lives in a heap-allocated vector; **after** is the inline-storage revision ([43aa64f](https://github.com/mmurshed/biginteger/tree/43aa64f), #156), which stores values of up to 4 limbs inside the object and computes operands of up to 16 limbs in a stack buffer. Each figure is the fastest of seven rounds of 400,000 calls. BigMath operators return a new value each call; the GMP loop writes into an output it reuses, which is GMP's usual pattern. Division and remainder use a 2n-limb dividend and an n-limb divisor. [CSV](docs/benchmarks/2026-11-21-small-operands.csv)

![Six plots of BigMath time divided by GMP time on 1 to 16 limb operands, before and after inline storage](docs/images/bigmath_vs_gmp_small.png)

| Operation | Limbs | Before ns | After ns | GMP ns | Before ratio | After ratio |
|---|---:|---:|---:|---:|---:|---:|
| Multiply | 1 | 21.1 | 8.0 | 3.8 | 5.71× | 2.11× |
| Multiply | 2 | 31.0 | 10.6 | 9.1 | 3.44× | 1.17× |
| Multiply | 4 | 41.1 | 29.4 | 14.7 | 2.82× | 2.00× |
| Multiply | 8 | 73.0 | 60.9 | 39.5 | 1.87× | 1.54× |
| Multiply | 16 | 170.6 | 152.5 | 128.1 | 1.31× | 1.19× |
| Add | 1 | 22.5 | 9.6 | 3.7 | 6.25× | 2.62× |
| Add | 2 | 23.4 | 9.8 | 3.7 | 6.27× | 2.63× |
| Add | 4 | 24.9 | 10.6 | 4.0 | 6.06× | 2.65× |
| Add | 8 | 27.1 | 26.2 | 4.8 | 5.82× | 5.46× |
| Add | 16 | 29.1 | 30.3 | 6.2 | 4.70× | 4.89× |
| Subtract | 1 | 23.6 | 10.9 | 3.7 | 6.12× | 2.92× |
| Subtract | 2 | 25.0 | 11.6 | 3.9 | 6.02× | 2.96× |
| Subtract | 4 | 26.3 | 12.4 | 4.3 | 6.05× | 2.86× |
| Subtract | 8 | 27.1 | 25.9 | 5.2 | 5.36× | 4.95× |
| Subtract | 16 | 31.4 | 30.0 | 6.8 | 4.59× | 4.43× |
| Divide (2n ÷ n) | 1 | 34.6 | 16.3 | 20.1 | 2.90× | 0.81× |
| Divide (2n ÷ n) | 2 | 68.2 | 47.8 | 22.7 | 3.02× | 2.10× |
| Divide (2n ÷ n) | 4 | 94.0 | 89.8 | 63.1 | 1.49× | 1.42× |
| Divide (2n ÷ n) | 8 | 170.7 | 163.9 | 113.5 | 1.50× | 1.44× |
| Divide (2n ÷ n) | 16 | 368.2 | 367.3 | 225.2 | 1.61× | 1.63× |
| Remainder (2n ÷ n) | 1 | 49.0 | 14.1 | 12.5 | 3.94× | 1.13× |
| Remainder (2n ÷ n) | 2 | 88.2 | 48.1 | 24.3 | 3.67× | 1.98× |
| Remainder (2n ÷ n) | 4 | 115.4 | 77.8 | 70.5 | 1.64× | 1.10× |
| Remainder (2n ÷ n) | 8 | 194.3 | 166.8 | 134.9 | 1.44× | 1.24× |
| Remainder (2n ÷ n) | 16 | 397.6 | 398.6 | 300.9 | 1.31× | 1.32× |
| Copy | 1 | 13.2 | 0.3 | 2.5 | 5.32× | 0.12× |
| Copy | 2 | 13.5 | 0.3 | 2.5 | 5.43× | 0.12× |
| Copy | 4 | 15.0 | 0.3 | 2.8 | 5.36× | 0.11× |
| Copy | 8 | 14.5 | 14.3 | 2.8 | 5.19× | 5.13× |
| Copy | 16 | 15.8 | 15.7 | 3.8 | 4.22× | 4.16× |

Values of up to 4 limbs no longer allocate: a copy takes 0.3 ns instead of 13 ns and is faster than GMP's `mpz_set`, and a 1-limb product takes 8 ns instead of 21 ns. Above 4 limbs every result still needs one allocation, so addition stays about 5× behind GMP there; the products and quotients gain only where the stack path replaced the dispatcher's temporary vectors.

To repeat: `cmake --build build --target small_ops_vs_gmp && ./build/small_ops_vs_gmp after`.

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
| #159–#160 | **large-NTT pass** (2026-11-21): padded MFA transpose tiles (the tile rows sat one L1 way apart, 16 KB, and thrashed the 8-way cache), and in single-thread builds the truncated `H + L` transform extended to `H = 2^25` on cache-blocked transforms, one prime at a time | mul 50M–200M digits −19–24% threaded; single-threaded 50M digits 3.31× → 1.24× GMP's time, 200M 3.75× → 1.52× |
| #162 | **ToString small levels** (2026-11-21): the power-of-ten split divides levels below 13,000 digits (6,500 threaded) with the regular division code and builds no reciprocal for them | warm ToString 10k digits 2.23× → 1.42× GMP single-threaded; first call 3.11× → 1.61×; beats GMP at 1M–2M digits single-threaded |

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
