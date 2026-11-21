# BigMath vs GMP benchmark results

## Single-thread comparison

On November 21, 2026, BigMath and GMP 6.3.0 each used one thread on the Apple M1 Max (10 cores, 32 GiB). The full split benchmark ran the two libraries sequentially on the same 91 operand sets as the default-build run below, and both comparisons use the same GMP pass. The Release build used Apple Clang 21, native arm64 instructions and the default 64-bit limbs. The library and its benchmark executable were both compiled with `-DBIGMATH_USE_THREADS=0`. The measured revision adds inline storage for small values (at most 4 limbs), which none of these operands are, so these results are expected to match the previous run within noise. [Measured revision](https://github.com/mmurshed/biginteger/tree/43aa64f)

The full CSV includes the first call, fastest subsequent call, warm repetition count and arithmetic check status for all 111 rows. Most large sizes have one warm observation per library. All 24 product checks and 40 division checks matched GMP’s reference hashes. Each division mode was checked against the corresponding GMP operation. [README performance table](README.md#performance-numbers) · [Single-thread CSV](docs/benchmarks/2026-11-21-sbo-m1-max-single-thread.csv) · [Single-thread figure](docs/images/bigmath_vs_gmp_single_thread.png)

To repeat the single-thread run, give the driver its own build and results directory. The figure generator uses the checked-in CSVs and requires Matplotlib.

```sh
CXXFLAGS=-DBIGMATH_USE_THREADS=0 BENCH_DIR="$PWD/.benchcache/single-thread" \
  HISTORY_MD="$PWD/.benchcache/single-thread/history.md" \
  bash tests/performance/benchsuite/run_benchmark.sh full --force-gmp
python3 docs/images/make_benchmark_plots.py --ratio-only
```

## Default build with up to eight BigMath threads

The November 21, 2026 run measured the inline-storage revision on an Apple M1 Max, with balanced multiplication through 200 million digits per operand and skewed division through 200 million ÷ 40 million digits. Both libraries produced the requested outputs on the same inputs. All 24 product checks and 40 division checks matched the GMP reference hashes. [Full results CSV](docs/benchmarks/2026-11-21-sbo-m1-max.csv) · [Measured revision](https://github.com/mmurshed/biginteger/tree/43aa64f) · [Previous run (PR #152)](docs/benchmarks/2026-11-21-pr152-m1-max.csv)

**Run-to-run variation.** The BigMath pass was run three times on the same build. The published CSV is the third run, the only one taken against a GMP pass whose small-size rows matched GMP in a tight loop (an earlier GMP pass measured 1.3 µs for the 1,000-digit product where a loop measures 0.88 µs). Rows with a single warm observation vary between runs. The largest spread was multithreaded parsing: 0.57–1.02 ms at 100,000 digits, 2.96–4.55 ms at 500,000 and 6.4–9.9 ms at 1,000,000 across the three runs, against 0.54, 3.46 and 6.39 ms in the previous published run. Single-thread parsing did not vary, which points at thread-pool timing rather than the arithmetic.

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
| 1,000 × 1,000 | 0.008 | 0.1106 | 0.07× | 0.00125 | 0.000958 | 1.30× | 7 |
| 5,000 × 5,000 | 0.02308 | 0.03571 | 0.65× | 0.01746 | 0.01246 | 1.40× | 5 |
| 10,000 × 10,000 | 0.2233 | 0.03083 | 7.24× | 0.02446 | 0.02792 | 0.88× | 5 |
| 50,000 × 50,000 | 0.4817 | 0.2775 | 1.74× | 0.1307 | 0.2717 | 0.48× | 3 |
| 100,000 × 100,000 | 0.8899 | 0.6574 | 1.35× | 0.232 | 0.5921 | 0.39× | 3 |
| 500,000 × 500,000 | 2.303 | 4.391 | 0.52× | 1.288 | 4.261 | 0.30× | 1 |
| 1,000,000 × 1,000,000 | 4.842 | 8.986 | 0.54× | 2.787 | 8.977 | 0.31× | 1 |
| 2,000,000 × 2,000,000 | 10.83 | 20.43 | 0.53× | 7.238 | 20.62 | 0.35× | 1 |
| 5,000,000 × 5,000,000 | 18.8 | 62.68 | 0.30× | 12.16 | 62.29 | 0.20× | 1 |
| 10,000,000 × 10,000,000 | 42.45 | 207.4 | 0.20× | 28.85 | 201.8 | 0.14× | 1 |
| 20,000,000 × 20,000,000 | 99.7 | 268.8 | 0.37× | 72.37 | 267 | 0.27× | 1 |
| 50,000,000 × 50,000,000 | 559.8 | 666.9 | 0.84× | 447.8 | 659.9 | 0.68× | 1 |
| 100,000,000 × 100,000,000 | 1,216 | 1,399 | 0.87× | 973.4 | 1,390 | 0.70× | 1 |
| 200,000,000 × 200,000,000 | 4,134 | 2,886 | 1.43× | 2,176 | 2,865 | 0.76× | 1 |

#### Skewed lengths

| Operand digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 100,000 × 10,000 | 0.1532 | 0.3157 | 0.49× | 0.1158 | 0.2985 | 0.39× | 3 |
| 500,000 × 50,000 | 0.4852 | 2.065 | 0.23× | 0.476 | 2.085 | 0.23× | 1 |
| 1,000,000 × 100,000 | 1.348 | 4.473 | 0.30× | 1.191 | 4.479 | 0.27× | 1 |
| 2,000,000 × 200,000 | 3.234 | 9.374 | 0.34× | 2.939 | 9.34 | 0.31× | 1 |
| 5,000,000 × 500,000 | 9.689 | 30.86 | 0.31× | 8.424 | 30.79 | 0.27× | 1 |
| 10,000,000 × 1,000,000 | 21.48 | 70.32 | 0.31× | 18.85 | 69.58 | 0.27× | 1 |
| 20,000,000 × 2,000,000 | 71.78 | 158.9 | 0.45× | 67.56 | 159 | 0.42× | 1 |
| 50,000,000 × 5,000,000 | 200.7 | 665.1 | 0.30× | 199.1 | 661.7 | 0.30× | 1 |
| 100,000,000 × 10,000,000 | 466 | 999.4 | 0.47× | 451.5 | 996.4 | 0.45× | 1 |
| 200,000,000 × 20,000,000 | 1,030 | 1,910 | 0.54× | 1,006 | 1,911 | 0.53× | 1 |

### Division

Quotient-only uses BigMath `/` and GMP `mpz_tdiv_q`. Quotient and remainder uses BigMath `DivideAndRemainder` and GMP `mpz_tdiv_qr`. Both modes use the same dividend and divisor. Equal-digit operands can produce a quotient of zero or one, so those results depend on the sampled values.

#### Skewed dividends

| Dividend ÷ divisor digits | BigMath Q ms | GMP Q ms | Q ratio | BigMath QR ms | GMP QR ms | QR ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 40,000 ÷ 10,000 | 0.2858 | 0.2007 | 1.42× | 0.2909 | 0.2266 | 1.28× | 3 |
| 100,000 ÷ 10,000 | 0.8487 | 0.4285 | 1.98× | 0.8481 | 0.441 | 1.92× | 3 |
| 200,000 ÷ 50,000 | 1.064 | 1.512 | 0.70× | 1.146 | 1.684 | 0.68× | 1 |
| 500,000 ÷ 100,000 | 2.299 | 4.289 | 0.54× | 2.616 | 4.506 | 0.58× | 1 |
| 1,000,000 ÷ 200,000 | 4.283 | 9.546 | 0.45× | 4.372 | 10.13 | 0.43× | 1 |
| 2,000,000 ÷ 500,000 | 7.869 | 22.39 | 0.35× | 8.226 | 24.05 | 0.34× | 1 |
| 5,000,000 ÷ 1,000,000 | 25.52 | 63.51 | 0.40× | 27.4 | 67.95 | 0.40× | 1 |
| 10,000,000 ÷ 2,000,000 | 52.73 | 145 | 0.36× | 52.45 | 154.6 | 0.34× | 1 |
| 20,000,000 ÷ 4,000,000 | 103 | 320.6 | 0.32× | 101.2 | 344.5 | 0.29× | 1 |
| 50,000,000 ÷ 10,000,000 | 244.8 | 1,220 | 0.20× | 230.3 | 1,297 | 0.18× | 1 |
| 100,000,000 ÷ 20,000,000 | 523.5 | 2,309 | 0.23× | 517.4 | 2,497 | 0.21× | 1 |
| 200,000,000 ÷ 40,000,000 | 1,356 | 4,189 | 0.32× | 1,341 | 4,441 | 0.30× | 1 |

At 200,000 ÷ 50,000 digits, the in-process run measured warm ratios of 0.92× for quotient-only and 1.04× for quotient and remainder. Each suite took one warm measurement per library at that size. Both suites favored BigMath in both modes from the next sampled shape, 500,000 ÷ 100,000, through 200 million ÷ 40 million digits. [In-process output](docs/benchmarks/2026-11-21-pr152-in-process.txt)

#### Equal digit counts

| Dividend ÷ divisor digits | BigMath Q ms | GMP Q ms | Q ratio | BigMath QR ms | GMP QR ms | QR ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 ÷ 1,000 | 0.000083 | <1 tick | — | 0.000166 | 0.000041 | 4.05× | 7 |
| 5,000 ÷ 5,000 | 0.000084 | <1 tick | — | 0.000458 | 0.000208 | 2.20× | 5 |
| 10,000 ÷ 10,000 | 0.000041 | <1 tick | — | 0.000208 | 0.000083 | 2.51× | 5 |
| 50,000 ÷ 50,000 | 0.000125 | 0.000041 | 3.05× | 0.003917 | 0.001666 | 2.35× | 3 |
| 100,000 ÷ 100,000 | 0.000125 | 0.000125 | 1.00× | 0.008459 | 0.003542 | 2.39× | 3 |
| 500,000 ÷ 500,000 | 0.000084 | 0.000167 | 0.50× | 0.01025 | 0.004333 | 2.37× | 1 |
| 1,000,000 ÷ 1,000,000 | 0.000292 | 0.000208 | 1.40× | 0.01667 | 0.008625 | 1.93× | 1 |
| 5,000,000 ÷ 5,000,000 | 0.000209 | 0.000208 | 1.00× | 0.09942 | 0.04308 | 2.31× | 1 |

### Decimal parse (string → BigInteger)

| Digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.007958 | 0.003458 | 2.30× | 0.002459 | 0.001583 | 1.55× | 7 |
| 10,000 | 0.07779 | 0.04779 | 1.63× | 0.06517 | 0.03733 | 1.75× | 5 |
| 50,000 | 0.6757 | 0.4162 | 1.62× | 0.5443 | 0.3952 | 1.38× | 3 |
| 100,000 | 1.091 | 1.072 | 1.02× | 1.023 | 1.058 | 0.97× | 3 |
| 500,000 | 4.649 | 9.03 | 0.51× | 4.547 | 8.975 | 0.51× | 1 |
| 1,000,000 | 9.408 | 20.98 | 0.45× | 9.931 | 20.93 | 0.47× | 1 |
| 2,000,000 | 19.25 | 47.36 | 0.41× | 14.98 | 47.06 | 0.32× | 1 |
| 5,000,000 | 39.97 | 146.4 | 0.27× | 38.22 | 147.7 | 0.26× | 1 |
| 10,000,000 | 80.03 | 345.2 | 0.23× | 78.18 | 342.5 | 0.23× | 1 |
| 20,000,000 | 177.9 | 801.3 | 0.22× | 177.1 | 799.8 | 0.22× | 1 |
| 50,000,000 | 706.1 | 2,628 | 0.27× | 704.4 | 2,646 | 0.27× | 1 |

### ToString (BigInteger → string)

| Digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.00825 | 0.006667 | 1.24× | 0.0065 | 0.003416 | 1.90× | 7 |
| 10,000 | 0.4571 | 0.09833 | 4.65× | 0.1787 | 0.07687 | 2.33× | 5 |
| 50,000 | 2.085 | 0.8882 | 2.35× | 1.323 | 0.8487 | 1.56× | 3 |
| 100,000 | 2.169 | 2.281 | 0.95× | 1.341 | 2.271 | 0.59× | 1 |
| 200,000 | 4.506 | 6.021 | 0.75× | 2.71 | 5.981 | 0.45× | 1 |
| 500,000 | 10.48 | 20.22 | 0.52× | 6.499 | 20.14 | 0.32× | 1 |
| 1,000,000 | 20.2 | 48.84 | 0.41× | 14.35 | 48.72 | 0.29× | 1 |
| 2,000,000 | 44.41 | 117.8 | 0.38× | 31.47 | 117.5 | 0.27× | 1 |
| 5,000,000 | 161 | 378.3 | 0.43× | 120.7 | 377.5 | 0.32× | 1 |
| 10,000,000 | 360.8 | 892.5 | 0.40× | 269.7 | 892.5 | 0.30× | 1 |
| 20,000,000 | 798 | 2,099 | 0.38× | 605.3 | 2,096 | 0.29× | 1 |

### BigDecimal

BigMath uses decimal fixed-point arithmetic and HALF_EVEN rounding for division. The comparison uses GMP binary floating-point (`mpf_t`), with precision derived from the digit count and requested scale. The split runner times these operations without an equality check. A size of `1,000.100` denotes 1,000 integer digits and 100 fractional digits.

| Operation | Size | BigMath warm ms | GMP warm ms | Ratio | Warm reps |
|---|---:|---:|---:|---:|---:|
| add | 100.10 digits | 0.000084 | 0.000041 | 2.05× | 7 |
| mul | 100.10 digits | 0.000042 | 0.000041 | 1.02× | 7 |
| add | 1000.100 digits | 0.000125 | <1 tick | — | 7 |
| mul | 1000.100 digits | 0.0015 | 0.001041 | 1.44× | 7 |
| add | 5000.500 digits | 0.000333 | 0.000125 | 2.66× | 5 |
| mul | 5000.500 digits | 0.02033 | 0.01258 | 1.62× | 5 |
| add | 20000.2000 digits | 0.001166 | 0.0005 | 2.33× | 3 |
| mul | 20000.2000 digits | 0.0505 | 0.09375 | 0.54× | 3 |
| div | 100.10 / 10 dp | 0.000291 | 0.000084 | 3.46× | 7 |
| div | 1000.100 / 100 dp | 0.000916 | 0.001791 | 0.51× | 7 |
| div | 5000.500 / 500 dp | 0.009625 | 0.02475 | 0.39× | 5 |
| div | 20000.2000 / 2000 dp | 0.06571 | 0.1919 | 0.34× | 3 |
| div | 2000.200 / 0 dp | 0.0005 | 0.004875 | 0.10× | 5 |
| div | 2000.200 / 10 dp | 0.000708 | 0.004875 | 0.15× | 5 |
| div | 2000.200 / 100 dp | 0.001375 | 0.00525 | 0.26× | 5 |
| div | 2000.200 / 1000 dp | 0.007375 | 0.009041 | 0.82× | 5 |
| div | 2000.200 / 5000 dp | 0.02996 | 0.03288 | 0.91× | 5 |
| parse | 100 digits | 0.000375 | 0.000459 | 0.82× | 7 |
| parse | 1000 digits | 0.004625 | 0.004458 | 1.04× | 7 |
| parse | 10000 digits | 0.08733 | 0.1048 | 0.83× | 5 |
| parse | 50000 digits | 0.6247 | 1.001 | 0.62× | 3 |
| tostr | 100 digits | 0.000291 | 0.000458 | 0.64× | 7 |
| tostr | 1000 digits | 0.006625 | 0.004541 | 1.46× | 7 |
| tostr | 10000 digits | 0.1768 | 0.1013 | 1.74× | 5 |
| tostr | 50000 digits | 1.308 | 1.089 | 1.20× | 3 |

## Small operands

`tests/performance/small_ops_vs_gmp.cpp` times the BigInteger operators against GMP on operands of 1–16 limbs (up to about 300 decimal digits), where per-call overhead and allocation decide the cost. It ran single-threaded on the same machine against two builds: **before** is `main` at 3dafa62, where every value lives in a heap-allocated vector; **after** is the revision measured above, which stores values of up to 4 limbs inside the object and computes operands of up to 16 limbs in a stack buffer. Each figure is the fastest of seven rounds of 400,000 calls. BigMath operators return a new value each call; the GMP loop writes into an output it reuses, which is GMP's usual pattern. Division and remainder use a 2n-limb dividend and an n-limb divisor. [CSV](docs/benchmarks/2026-11-21-small-operands.csv)

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
