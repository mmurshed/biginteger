# BigMath vs GMP benchmark results

## Single-thread comparison

The November 21, 2026 full split benchmark measured `main` after PR #165 on an Apple M1 Max. The single-thread build compiled both BigMath and its benchmark driver with `-DBIGMATH_USE_THREADS=0`. With one thread per library, a 200-million-digit balanced product took 2.91 s in BigMath and 2.87 s in GMP. At 50 million and 100 million digits BigMath used 0.96× and 0.98× GMP’s time. Both BigMath configurations used this run’s single-thread GMP 6.3.0 reference. [Measured revision](https://github.com/mmurshed/biginteger/tree/8e19932) · [Previous single-thread run](docs/benchmarks/2026-11-21-tostring-m1-max-single-thread.csv)

The full CSV includes the first call, fastest subsequent call, warm repetition count and arithmetic check status for all 111 rows. Most large sizes have one warm observation per library. All 24 product checks and 40 division checks matched GMP’s reference hashes. Each division mode was checked against the corresponding GMP operation. [README performance table](README.md#performance-numbers) · [Single-thread CSV](docs/benchmarks/2026-11-21-factored-mfa-m1-max-single-thread.csv) · [Single-thread figure](docs/images/bigmath_vs_gmp_single_thread.png)

To repeat the single-thread run, give the driver its own build and results directory. The figure generator uses the checked-in CSVs and requires Matplotlib.

```sh
CXXFLAGS=-DBIGMATH_USE_THREADS=0 BENCH_DIR="$PWD/.benchcache/single-thread" \
  HISTORY_MD="$PWD/.benchcache/single-thread/history.md" \
  bash tests/performance/benchsuite/run_benchmark.sh full --force-gmp
python3 docs/images/make_benchmark_plots.py --ratio-only
```

## Default build with up to eight BigMath threads

With up to eight BigMath threads, equal-length products at 50, 100 and 200 million digits took 164.7, 345.3 and 719.5 ms. At 200 million digits the threaded build used 0.25× GMP’s time. All 64 arithmetic result hashes matched GMP in each build. [Full results CSV](docs/benchmarks/2026-11-21-factored-mfa-m1-max.csv) · [Measured revision](https://github.com/mmurshed/biginteger/tree/8e19932) · [Previous full run](docs/benchmarks/2026-11-21-tostring-m1-max.csv)

**Run-to-run variation.** Each full-profile product from 500,000 digits up has one warm observation. The focused before/after experiment used three warm calls on the same large operands. Its 200-million-digit threaded result was 677.5 ms, versus 719.5 ms in this full run. [Focused measurements](#large-multiplication-follow-up)

![Four plots of BigMath time divided by GMP time for multiplication, division, decimal parsing and formatting](docs/images/bigmath_vs_gmp.png)

The arithmetic panels show warm calls. The division panel covers dividend-to-divisor digit ratios from 4:1 to 10:1. The conversion panels show both the first and warm calls. The dashed horizontal line marks equal time.

### Methodology

- The host had 10 CPU cores and 32 GiB of memory and ran macOS 27.2 (Darwin 27.2.0). Apple Clang 21 compiled the C++20 Release build with `-O3 -march=native -DNDEBUG`. BigMath used its default 64-bit limbs, CRT NTT, NEON and up to eight threads. GMP 6.3.0 ran single-threaded.
- This run reused the full profile’s 91 deterministic operand sets. GMP ran first, followed by the default and single-thread BigMath builds. Each process read the same decimal operands from disk. Each of the 20 division sets was timed in both output modes, giving 111 result rows. Parsing operands, disk I/O and result hashes are outside the arithmetic timers.
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
| 1,000 × 1,000 | 0.03479 | 0.01063 | 3.27× | 0.001833 | 0.001333 | 1.38× | 7 |
| 5,000 × 5,000 | 0.04688 | 0.02262 | 2.07× | 0.024 | 0.01508 | 1.59× | 5 |
| 10,000 × 10,000 | 0.539 | 0.041 | 13.15× | 0.03067 | 0.03679 | 0.83× | 5 |
| 50,000 × 50,000 | 0.631 | 0.3656 | 1.73× | 0.1295 | 0.36 | 0.36× | 3 |
| 100,000 × 100,000 | 1.056 | 0.7795 | 1.35× | 0.2641 | 0.7124 | 0.37× | 3 |
| 500,000 × 500,000 | 2.321 | 4.434 | 0.52× | 1.236 | 4.379 | 0.28× | 1 |
| 1,000,000 × 1,000,000 | 4.794 | 9.239 | 0.52× | 3.031 | 8.997 | 0.34× | 1 |
| 2,000,000 × 2,000,000 | 10.61 | 20.84 | 0.51× | 7.425 | 20.56 | 0.36× | 1 |
| 5,000,000 × 5,000,000 | 11.08 | 63.39 | 0.17× | 8.139 | 62.65 | 0.13× | 1 |
| 10,000,000 × 10,000,000 | 23.97 | 204.6 | 0.12× | 17.74 | 201.7 | 0.09× | 1 |
| 20,000,000 × 20,000,000 | 49.96 | 269.8 | 0.19× | 36.16 | 266 | 0.14× | 1 |
| 50,000,000 × 50,000,000 | 218.8 | 663.2 | 0.33× | 164.7 | 655 | 0.25× | 1 |
| 100,000,000 × 100,000,000 | 447.4 | 1,408 | 0.32× | 345.3 | 1,386 | 0.25× | 1 |
| 200,000,000 × 200,000,000 | 961.3 | 2,907 | 0.33× | 719.5 | 2,869 | 0.25× | 1 |

#### Skewed lengths

| Operand digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 100,000 × 10,000 | 0.1292 | 0.3245 | 0.40× | 0.1248 | 0.2997 | 0.42× | 3 |
| 500,000 × 50,000 | 0.5576 | 2.131 | 0.26× | 0.4331 | 2.071 | 0.21× | 1 |
| 1,000,000 × 100,000 | 1.095 | 4.476 | 0.24× | 1.07 | 4.437 | 0.24× | 1 |
| 2,000,000 × 200,000 | 2.896 | 9.416 | 0.31× | 2.453 | 9.263 | 0.26× | 1 |
| 5,000,000 × 500,000 | 7.349 | 31.3 | 0.23× | 8.403 | 31.16 | 0.27× | 1 |
| 10,000,000 × 1,000,000 | 21.07 | 71.01 | 0.30× | 18.93 | 70.65 | 0.27× | 1 |
| 20,000,000 × 2,000,000 | 37.79 | 158.3 | 0.24× | 34.98 | 157.9 | 0.22× | 1 |
| 50,000,000 × 5,000,000 | 77.3 | 671.8 | 0.12× | 78.47 | 675.2 | 0.12× | 1 |
| 100,000,000 × 10,000,000 | 170.7 | 1,004 | 0.17× | 173.6 | 995.1 | 0.17× | 1 |
| 200,000,000 × 20,000,000 | 362.2 | 1,914 | 0.19× | 352.1 | 1,901 | 0.19× | 1 |

### Division

Quotient-only uses BigMath `/` and GMP `mpz_tdiv_q`. Quotient and remainder uses BigMath `DivideAndRemainder` and GMP `mpz_tdiv_qr`. Both modes use the same dividend and divisor. Equal-digit operands can produce a quotient of zero or one, so those results depend on the sampled values.

#### Skewed dividends

| Dividend ÷ divisor digits | BigMath Q ms | GMP Q ms | Q ratio | BigMath QR ms | GMP QR ms | QR ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 40,000 ÷ 10,000 | 0.286 | 0.1981 | 1.44× | 0.2856 | 0.2228 | 1.28× | 3 |
| 100,000 ÷ 10,000 | 0.849 | 0.425 | 2.00× | 0.8472 | 0.4388 | 1.93× | 3 |
| 200,000 ÷ 50,000 | 0.9902 | 1.508 | 0.66× | 1.337 | 1.639 | 0.82× | 1 |
| 500,000 ÷ 100,000 | 2.418 | 4.215 | 0.57× | 2.357 | 4.491 | 0.52× | 1 |
| 1,000,000 ÷ 200,000 | 4.236 | 9.501 | 0.45× | 4.923 | 10.1 | 0.49× | 1 |
| 2,000,000 ÷ 500,000 | 8.265 | 22.74 | 0.36× | 8.763 | 24.5 | 0.36× | 1 |
| 5,000,000 ÷ 1,000,000 | 25.75 | 65.53 | 0.39× | 25.55 | 70.02 | 0.36× | 1 |
| 10,000,000 ÷ 2,000,000 | 57.94 | 142.3 | 0.41× | 55.59 | 152 | 0.37× | 1 |
| 20,000,000 ÷ 4,000,000 | 88.69 | 319.4 | 0.28× | 86.14 | 340.7 | 0.25× | 1 |
| 50,000,000 ÷ 10,000,000 | 157.2 | 1,222 | 0.13× | 159.7 | 1,303 | 0.12× | 1 |
| 100,000,000 ÷ 20,000,000 | 328.1 | 2,298 | 0.14× | 326.2 | 2,503 | 0.13× | 1 |
| 200,000,000 ÷ 40,000,000 | 690 | 4,193 | 0.16× | 687.4 | 4,436 | 0.15× | 1 |

At 200,000 ÷ 50,000 digits the current full run measured 0.66× GMP’s time for the quotient and 0.82× for quotient and remainder. The earlier in-process run on PR #152 measured 0.92× and 1.04× on different operands. Both suites took one warm measurement per library at that size. [In-process output](docs/benchmarks/2026-11-21-pr152-in-process.txt)

#### Equal digit counts

| Dividend ÷ divisor digits | BigMath Q ms | GMP Q ms | Q ratio | BigMath QR ms | GMP QR ms | QR ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 ÷ 1,000 | 0.000083 | <1 tick | — | 0.000166 | 0.000083 | 2.00× | 7 |
| 5,000 ÷ 5,000 | 0.000083 | <1 tick | — | 0.000416 | 0.000208 | 2.00× | 5 |
| 10,000 ÷ 10,000 | 0.000041 | <1 tick | — | 0.000208 | 0.000083 | 2.51× | 5 |
| 50,000 ÷ 50,000 | 0.000083 | <1 tick | — | 0.003959 | 0.001666 | 2.38× | 3 |
| 100,000 ÷ 100,000 | 0.000083 | 0.000125 | 0.66× | 0.008125 | 0.003583 | 2.27× | 3 |
| 500,000 ÷ 500,000 | 0.000166 | 0.000208 | 0.80× | 0.008042 | 0.004292 | 1.87× | 1 |
| 1,000,000 ÷ 1,000,000 | 0.000083 | 0.000208 | 0.40× | 0.0175 | 0.008625 | 2.03× | 1 |
| 5,000,000 ÷ 5,000,000 | 0.0035 | 0.000208 | 16.83× | 0.1256 | 0.04317 | 2.91× | 1 |

### Decimal parse (string → BigInteger)

| Digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.005083 | 0.00825 | 0.62× | 0.0025 | 0.001584 | 1.58× | 7 |
| 10,000 | 0.07646 | 0.05325 | 1.44× | 0.06325 | 0.03929 | 1.61× | 5 |
| 50,000 | 0.5469 | 0.4226 | 1.29× | 0.5072 | 0.395 | 1.28× | 3 |
| 100,000 | 0.6892 | 1.177 | 0.59× | 0.5789 | 1.114 | 0.52× | 3 |
| 500,000 | 3.049 | 9.082 | 0.34× | 3.251 | 8.83 | 0.37× | 1 |
| 1,000,000 | 7.147 | 20.61 | 0.35× | 7.791 | 20.27 | 0.38× | 1 |
| 2,000,000 | 14.82 | 46.99 | 0.32× | 16.76 | 46.64 | 0.36× | 1 |
| 5,000,000 | 38.69 | 148.4 | 0.26× | 39.15 | 148.3 | 0.26× | 1 |
| 10,000,000 | 77.45 | 344.1 | 0.23× | 76.8 | 343.9 | 0.22× | 1 |
| 20,000,000 | 160 | 799.7 | 0.20× | 160.6 | 797.3 | 0.20× | 1 |
| 50,000,000 | 526.6 | 2,622 | 0.20× | 533.1 | 2,633 | 0.20× | 1 |

### ToString (BigInteger → string)

| Digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.007458 | 0.006167 | 1.21× | 0.00625 | 0.003375 | 1.85× | 7 |
| 10,000 | 0.1424 | 0.0965 | 1.48× | 0.112 | 0.07679 | 1.46× | 5 |
| 50,000 | 1.407 | 0.9015 | 1.56× | 1.021 | 0.8699 | 1.17× | 3 |
| 100,000 | 1.893 | 2.335 | 0.81× | 1.157 | 2.295 | 0.50× | 1 |
| 200,000 | 3.806 | 6.111 | 0.62× | 2.449 | 6.015 | 0.41× | 1 |
| 500,000 | 10.38 | 20.85 | 0.50× | 5.806 | 20.46 | 0.28× | 1 |
| 1,000,000 | 20.18 | 49.43 | 0.41× | 12.76 | 49.72 | 0.26× | 1 |
| 2,000,000 | 41.5 | 119.2 | 0.35× | 30.41 | 118.9 | 0.26× | 1 |
| 5,000,000 | 153 | 382.6 | 0.40× | 109.5 | 381.8 | 0.29× | 1 |
| 10,000,000 | 317.7 | 902.1 | 0.35× | 235.1 | 909.4 | 0.26× | 1 |
| 20,000,000 | 671.4 | 2,121 | 0.32× | 502 | 2,124 | 0.24× | 1 |

### BigDecimal

BigMath uses decimal fixed-point arithmetic and HALF_EVEN rounding for division. The comparison uses GMP binary floating-point (`mpf_t`), with precision derived from the digit count and requested scale. The split runner times these operations without an equality check. A size of `1,000.100` denotes 1,000 integer digits and 100 fractional digits.

| Operation | Size | BigMath warm ms | GMP warm ms | Ratio | Warm reps |
|---|---:|---:|---:|---:|---:|
| add | 100.10 digits | 0.000125 | <1 tick | — | 7 |
| mul | 100.10 digits | 0.000041 | 0.000041 | 1.00× | 7 |
| add | 1000.100 digits | 0.000125 | <1 tick | — | 7 |
| mul | 1000.100 digits | 0.001541 | 0.001 | 1.54× | 7 |
| add | 5000.500 digits | 0.000292 | 0.000125 | 2.34× | 5 |
| mul | 5000.500 digits | 0.01892 | 0.01283 | 1.47× | 5 |
| add | 20000.2000 digits | 0.001125 | 0.000459 | 2.45× | 3 |
| mul | 20000.2000 digits | 0.04617 | 0.09638 | 0.48× | 3 |
| div | 100.10 / 10 dp | 0.000333 | 0.000083 | 4.01× | 7 |
| div | 1000.100 / 100 dp | 0.000916 | 0.001875 | 0.49× | 7 |
| div | 5000.500 / 500 dp | 0.009541 | 0.02546 | 0.37× | 5 |
| div | 20000.2000 / 2000 dp | 0.06604 | 0.1927 | 0.34× | 3 |
| div | 2000.200 / 0 dp | 0.000459 | 0.005042 | 0.09× | 5 |
| div | 2000.200 / 10 dp | 0.000708 | 0.005208 | 0.14× | 5 |
| div | 2000.200 / 100 dp | 0.001333 | 0.005458 | 0.24× | 5 |
| div | 2000.200 / 1000 dp | 0.007416 | 0.009292 | 0.80× | 5 |
| div | 2000.200 / 5000 dp | 0.02979 | 0.03383 | 0.88× | 5 |
| parse | 100 digits | 0.000416 | 0.0005 | 0.83× | 7 |
| parse | 1000 digits | 0.0045 | 0.004583 | 0.98× | 7 |
| parse | 10000 digits | 0.08525 | 0.1105 | 0.77× | 5 |
| parse | 50000 digits | 0.6339 | 0.9925 | 0.64× | 3 |
| tostr | 100 digits | 0.000292 | 0.000459 | 0.64× | 7 |
| tostr | 1000 digits | 0.006334 | 0.004792 | 1.32× | 7 |
| tostr | 10000 digits | 0.1129 | 0.1052 | 1.07× | 5 |
| tostr | 50000 digits | 1.02 | 1.098 | 0.93× | 3 |

## Large multiplication follow-up

The default threaded build uses a power-of-two NTT in this range. Two equal operands of `m` 64-bit limbs need `4m − 1` coefficients. The transform doubles just above 10.10M, 20.20M and 40.40M decimal digits per operand. The 10M and 20M products fill 99% of their transforms. The 50M, 100M and 200M products fill 62%. A sweep around the first boundary on the M1 Max measured 27.3 ms at 10M and 65.1 ms at 11M on revision `1fa44fa`, while GMP took 200.7 and 200.6 ms on those same inputs (median of three warm products, seeded binary operands).

The matrix Fourier algorithm also paid more per coefficient at larger sizes. It accessed a large root table at a stride that increases with the row index. The revised NEON kernel splits each bit-reversed frequency into a block factor and a lane factor and multiplies four adjacent coefficients at a time. At 200M digits in the threaded build, the top-level forward and inverse tables for the three primes together shrink from 2,304 MiB to 108 MiB, calculated from their table lengths. These figures exclude the transform buffers and leaf plans. Scheduling transpose tiles across the pool lets all eight threads work in stages previously split into six tasks.

On a seeded 200M-digit operand pair, the original code took 1,677 ms. Factoring the cross-twiddles and storing the compact tables brought that to 846 ms with the original six-task schedule. Distributing the tiles across all eight threads brought it to 721 ms. These are warm medians (three, three and five repetitions respectively). Each product matched GMP limb-for-limb.

The split benchmark compared the public multiplication operator before and after the change on the same saved decimal operands. Each entry below is the fastest of three warm products in milliseconds. GMP 6.3.0 used one thread. All four BigMath builds matched the GMP result hashes on all nine operand pairs. [Balanced and 10:1 results](docs/benchmarks/2026-11-21-xl-multiplication-m1-max.csv)

| Digits per operand | 8 threads before | 8 threads after | 1 thread before | 1 thread after | GMP |
|---:|---:|---:|---:|---:|---:|
| 10M | 27.6 | 15.4 | 130.4 | 109.7 | 201.6 |
| 20M | 62.3 | 35.2 | 257.3 | 218.3 | 264.4 |
| 50M | 338.0 | 163.0 | 823.0 | 618.6 | 655.5 |
| 100M | 772.6 | 327.5 | 1,833.9 | 1,341.8 | 1,384.2 |
| 200M | 1,710.4 | 677.5 | 4,366.8 | 2,834.6 | 2,860.4 |

The 10M × 1M skewed case uses the smaller truncated transform. Its threaded time went from 16.1 to 17.6 ms in this sweep. Four isolated reruns per build, alternating their order, gave medians of 18.5 and 19.0 ms. The larger 50M × 5M and 200M × 20M cases used the revised matrix transform and improved by 1.96× and 2.36× respectively.

To repeat the comparison, select multiplication cases 8–13, 19, 21 and 23 from the full split-benchmark manifest, set their warm repetition counts to three, and use the same operand files for each build. Run `bench_run_gmp` once, then `bench_run_bigmath` from each Release build against that reference. The baseline is `1fa44fa`. The one-thread builds use `-DBIGMATH_USE_THREADS=0`.

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
| #165 | Factored NEON cross-twiddles and tile scheduling across the pool | On identical 200M-digit operands, threaded time fell from 1,710.4 to 677.5 ms and single-thread time from 4,366.8 to 2,834.6 ms |

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
