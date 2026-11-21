# BigMath vs GMP benchmark results

## Single-thread comparison

On November 21, 2026, BigMath and GMP 6.3.0 each used one thread on the Apple M1 Max (10 cores, 32 GiB). The full split benchmark ran the two libraries sequentially on the same 91 operand sets as the default-build run below, and both comparisons use the same GMP pass. The Release build used Apple Clang 21, native arm64 instructions and the default 64-bit limbs. The library and its benchmark executable were both compiled with `-DBIGMATH_USE_THREADS=0`. The measured revision is `main` after the Toom-3 tier (#157), which splits Karatsuba levels three ways from 168 limbs (160 for squares) and stays below the NTT entry. Of these operands, only the 5,000-digit product (about 260 limbs) and the sub-products of the smaller divisions fall in its range; the 5,000-digit product took 17.4 µs against 16.9 µs in the previous run, matching the 1–2.5% cost the Toom-3 PR measured near 256 × 256 limbs. The rest of the single-thread rows of 10 µs or longer moved by −1.2% at the median. [Measured revision](https://github.com/mmurshed/biginteger/tree/ff7df9f) · [Previous single-thread run](docs/benchmarks/2026-11-21-sbo-m1-max-single-thread.csv)

The full CSV includes the first call, fastest subsequent call, warm repetition count and arithmetic check status for all 111 rows. Most large sizes have one warm observation per library. All 24 product checks and 40 division checks matched GMP’s reference hashes. Each division mode was checked against the corresponding GMP operation. [README performance table](README.md#performance-numbers) · [Single-thread CSV](docs/benchmarks/2026-11-21-toom-m1-max-single-thread.csv) · [Single-thread figure](docs/images/bigmath_vs_gmp_single_thread.png)

To repeat the single-thread run, give the driver its own build and results directory. The figure generator uses the checked-in CSVs and requires Matplotlib.

```sh
CXXFLAGS=-DBIGMATH_USE_THREADS=0 BENCH_DIR="$PWD/.benchcache/single-thread" \
  HISTORY_MD="$PWD/.benchcache/single-thread/history.md" \
  bash tests/performance/benchsuite/run_benchmark.sh full --force-gmp
python3 docs/images/make_benchmark_plots.py --ratio-only
```

## Default build with up to eight BigMath threads

The November 21, 2026 run measured `main` after the Toom-3 tier (#157) on an Apple M1 Max, with balanced multiplication through 200 million digits per operand and skewed division through 200 million ÷ 40 million digits. Both libraries produced the requested outputs on the same inputs. All 24 product checks and 40 division checks matched the GMP reference hashes. [Full results CSV](docs/benchmarks/2026-11-21-toom-m1-max.csv) · [Measured revision](https://github.com/mmurshed/biginteger/tree/ff7df9f) · [Previous run (inline storage, #156)](docs/benchmarks/2026-11-21-sbo-m1-max.csv) · [Run before that (PR #152)](docs/benchmarks/2026-11-21-pr152-m1-max.csv)

**Run-to-run variation.** The published CSVs come from the second of two runs of both passes. The first was discarded: a malware scan and a file indexer ran during it, and GMP measured 70–85% slow on the smallest products. In the published run GMP's timings are within 5% of the previous GMP pass on all 77 rows of 10 µs or longer, and within 1.5% on 69 of them. Rows with a single warm observation still vary between runs, and the multithreaded ones most. Against the previous run, 8-thread rows that run entirely on the NTT (which the Toom-3 change does not touch) moved by up to ±18%: products of 50,000–500,000 digits were 9–12% faster, while 1,000,000 ÷ 200,000-digit division was 14–18% slower. Multithreaded parsing measured 0.57, 2.95 and 6.61 ms at 100,000, 500,000 and 1,000,000 digits, at or near the fast end of the 0.57–1.02, 2.96–4.55 and 6.4–9.9 ms spread of the three previous runs. Single-thread rows moved much less: 75 of the 77 rows of 10 µs or longer changed by under 5%, the exceptions being the 1,000,000-digit product (−6.8%) and the 5,000,000-digit equal-length quotient and remainder (−24%, a 0.07 ms row). That points at thread-pool timing rather than the arithmetic.

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
| 1,000 × 1,000 | 0.1365 | 0.1118 | 1.22× | 0.001291 | 0.000959 | 1.35× | 7 |
| 5,000 × 5,000 | 0.02354 | 0.03833 | 0.61× | 0.01692 | 0.01221 | 1.39× | 5 |
| 10,000 × 10,000 | 0.3725 | 0.03175 | 11.73× | 0.02242 | 0.02804 | 0.80× | 5 |
| 50,000 × 50,000 | 0.4671 | 0.2815 | 1.66× | 0.1147 | 0.2734 | 0.42× | 3 |
| 100,000 × 100,000 | 0.8768 | 0.6323 | 1.39× | 0.2097 | 0.588 | 0.36× | 3 |
| 500,000 × 500,000 | 2.047 | 4.277 | 0.48× | 1.146 | 4.234 | 0.27× | 1 |
| 1,000,000 × 1,000,000 | 4.741 | 8.949 | 0.53× | 2.789 | 8.847 | 0.32× | 1 |
| 2,000,000 × 2,000,000 | 10.7 | 20.31 | 0.53× | 7.213 | 20.33 | 0.35× | 1 |
| 5,000,000 × 5,000,000 | 18.64 | 66.13 | 0.28× | 11.77 | 62.79 | 0.19× | 1 |
| 10,000,000 × 10,000,000 | 41.81 | 202.8 | 0.21× | 27.5 | 199.9 | 0.14× | 1 |
| 20,000,000 × 20,000,000 | 94.47 | 266.4 | 0.35× | 68.37 | 265.3 | 0.26× | 1 |
| 50,000,000 × 50,000,000 | 555.3 | 660.8 | 0.84× | 449 | 650.7 | 0.69× | 1 |
| 100,000,000 × 100,000,000 | 1,231 | 1,393 | 0.88× | 969.8 | 1,376 | 0.70× | 1 |
| 200,000,000 × 200,000,000 | 4,209 | 2,862 | 1.47× | 2,130 | 2,828 | 0.75× | 1 |

#### Skewed lengths

| Operand digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 100,000 × 10,000 | 0.1569 | 0.3218 | 0.49× | 0.1023 | 0.3015 | 0.34× | 3 |
| 500,000 × 50,000 | 0.507 | 2.068 | 0.25× | 0.4612 | 2.039 | 0.23× | 1 |
| 1,000,000 × 100,000 | 1.206 | 4.455 | 0.27× | 1.258 | 4.429 | 0.28× | 1 |
| 2,000,000 × 200,000 | 3.23 | 9.312 | 0.35× | 3.216 | 9.244 | 0.35× | 1 |
| 5,000,000 × 500,000 | 9.547 | 30.85 | 0.31× | 8.176 | 30.84 | 0.27× | 1 |
| 10,000,000 × 1,000,000 | 22.61 | 70.04 | 0.32× | 18.5 | 69.68 | 0.27× | 1 |
| 20,000,000 × 2,000,000 | 69.28 | 158.6 | 0.44× | 65.35 | 158.2 | 0.41× | 1 |
| 50,000,000 × 5,000,000 | 195.2 | 662.8 | 0.29× | 217.1 | 660.4 | 0.33× | 1 |
| 100,000,000 × 10,000,000 | 443.4 | 996.3 | 0.45× | 455.3 | 994.4 | 0.46× | 1 |
| 200,000,000 × 20,000,000 | 957.3 | 1,893 | 0.51× | 993.7 | 1,890 | 0.53× | 1 |

### Division

Quotient-only uses BigMath `/` and GMP `mpz_tdiv_q`. Quotient and remainder uses BigMath `DivideAndRemainder` and GMP `mpz_tdiv_qr`. Both modes use the same dividend and divisor. Equal-digit operands can produce a quotient of zero or one, so those results depend on the sampled values.

#### Skewed dividends

| Dividend ÷ divisor digits | BigMath Q ms | GMP Q ms | Q ratio | BigMath QR ms | GMP QR ms | QR ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 40,000 ÷ 10,000 | 0.2859 | 0.2001 | 1.43× | 0.2856 | 0.2249 | 1.27× | 3 |
| 100,000 ÷ 10,000 | 0.8476 | 0.4268 | 1.99× | 0.8476 | 0.4409 | 1.92× | 3 |
| 200,000 ÷ 50,000 | 0.9873 | 1.507 | 0.66× | 1.01 | 1.626 | 0.62× | 1 |
| 500,000 ÷ 100,000 | 2.228 | 4.207 | 0.53× | 2.263 | 4.491 | 0.50× | 1 |
| 1,000,000 ÷ 200,000 | 5.075 | 9.476 | 0.54× | 4.995 | 10.08 | 0.50× | 1 |
| 2,000,000 ÷ 500,000 | 8.49 | 22.23 | 0.38× | 8.41 | 23.92 | 0.35× | 1 |
| 5,000,000 ÷ 1,000,000 | 24.27 | 63.21 | 0.38× | 24.61 | 67.66 | 0.36× | 1 |
| 10,000,000 ÷ 2,000,000 | 53.64 | 144.5 | 0.37× | 54.26 | 154.2 | 0.35× | 1 |
| 20,000,000 ÷ 4,000,000 | 100.9 | 319.7 | 0.32× | 103.7 | 342.2 | 0.30× | 1 |
| 50,000,000 ÷ 10,000,000 | 223.5 | 1,213 | 0.18× | 225.9 | 1,276 | 0.18× | 1 |
| 100,000,000 ÷ 20,000,000 | 516.4 | 2,288 | 0.23× | 517.4 | 2,490 | 0.21× | 1 |
| 200,000,000 ÷ 40,000,000 | 1,312 | 4,137 | 0.32× | 1,302 | 4,452 | 0.29× | 1 |

At 200,000 ÷ 50,000 digits, the in-process run measured warm ratios of 0.92× for quotient-only and 1.04× for quotient and remainder. Each suite took one warm measurement per library at that size. Both suites favored BigMath in both modes from the next sampled shape, 500,000 ÷ 100,000, through 200 million ÷ 40 million digits. [In-process output](docs/benchmarks/2026-11-21-pr152-in-process.txt)

#### Equal digit counts

| Dividend ÷ divisor digits | BigMath Q ms | GMP Q ms | Q ratio | BigMath QR ms | GMP QR ms | QR ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 ÷ 1,000 | 0.000084 | <1 tick | — | 0.000166 | 0.000042 | 3.95× | 7 |
| 5,000 ÷ 5,000 | 0.000083 | <1 tick | — | 0.000458 | 0.000167 | 2.74× | 5 |
| 10,000 ÷ 10,000 | 0.000041 | <1 tick | — | 0.000208 | 0.000042 | 4.95× | 5 |
| 50,000 ÷ 50,000 | 0.000125 | <1 tick | — | 0.003875 | 0.001666 | 2.33× | 3 |
| 100,000 ÷ 100,000 | 0.000125 | 0.000125 | 1.00× | 0.008375 | 0.003583 | 2.34× | 3 |
| 500,000 ÷ 500,000 | 0.000083 | 0.000167 | 0.50× | 0.008792 | 0.004292 | 2.05× | 1 |
| 1,000,000 ÷ 1,000,000 | 0.000083 | 0.000167 | 0.50× | 0.01512 | 0.008583 | 1.76× | 1 |
| 5,000,000 ÷ 5,000,000 | 0.000375 | 0.00025 | 1.50× | 0.09421 | 0.04288 | 2.20× | 1 |

### Decimal parse (string → BigInteger)

| Digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.005291 | 0.004542 | 1.16× | 0.002375 | 0.001541 | 1.54× | 7 |
| 10,000 | 0.07558 | 0.05375 | 1.41× | 0.06254 | 0.03771 | 1.66× | 5 |
| 50,000 | 0.569 | 0.4271 | 1.33× | 0.5257 | 0.3946 | 1.33× | 3 |
| 100,000 | 0.6699 | 1.078 | 0.62× | 0.5691 | 1.061 | 0.54× | 3 |
| 500,000 | 3.005 | 9.019 | 0.33× | 2.954 | 8.93 | 0.33× | 1 |
| 1,000,000 | 6.618 | 20.99 | 0.32× | 6.614 | 20.78 | 0.32× | 1 |
| 2,000,000 | 15.07 | 49.72 | 0.30× | 14.77 | 47.15 | 0.31× | 1 |
| 5,000,000 | 36.16 | 146.3 | 0.25× | 36.13 | 146 | 0.25× | 1 |
| 10,000,000 | 79.09 | 341.6 | 0.23× | 78.93 | 343.4 | 0.23× | 1 |
| 20,000,000 | 176 | 805.9 | 0.22× | 175.5 | 802.5 | 0.22× | 1 |
| 50,000,000 | 696.8 | 2,627 | 0.27× | 708.8 | 2,650 | 0.27× | 1 |

### ToString (BigInteger → string)

| Digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.0075 | 0.0065 | 1.15× | 0.006209 | 0.003375 | 1.84× | 7 |
| 10,000 | 0.3933 | 0.09275 | 4.24× | 0.1762 | 0.07733 | 2.28× | 5 |
| 50,000 | 1.787 | 0.888 | 2.01× | 1.285 | 0.8432 | 1.52× | 3 |
| 100,000 | 2.18 | 2.285 | 0.95× | 1.334 | 2.267 | 0.59× | 1 |
| 200,000 | 4.26 | 6.026 | 0.71× | 2.694 | 5.961 | 0.45× | 1 |
| 500,000 | 10.11 | 20.1 | 0.50× | 6.581 | 20.08 | 0.33× | 1 |
| 1,000,000 | 19.81 | 49.15 | 0.40× | 14.12 | 49.38 | 0.29× | 1 |
| 2,000,000 | 43.66 | 117.6 | 0.37× | 31.46 | 117.4 | 0.27× | 1 |
| 5,000,000 | 162.1 | 376.8 | 0.43× | 120.2 | 376.6 | 0.32× | 1 |
| 10,000,000 | 354.9 | 893.7 | 0.40× | 268.4 | 901.9 | 0.30× | 1 |
| 20,000,000 | 791.5 | 2,092 | 0.38× | 603.4 | 2,094 | 0.29× | 1 |

### BigDecimal

BigMath uses decimal fixed-point arithmetic and HALF_EVEN rounding for division. The comparison uses GMP binary floating-point (`mpf_t`), with precision derived from the digit count and requested scale. The split runner times these operations without an equality check. A size of `1,000.100` denotes 1,000 integer digits and 100 fractional digits.

| Operation | Size | BigMath warm ms | GMP warm ms | Ratio | Warm reps |
|---|---:|---:|---:|---:|---:|
| add | 100.10 digits | 0.000083 | <1 tick | — | 7 |
| mul | 100.10 digits | 0.000042 | 0.000041 | 1.02× | 7 |
| add | 1000.100 digits | 0.000125 | <1 tick | — | 7 |
| mul | 1000.100 digits | 0.001541 | 0.001042 | 1.48× | 7 |
| add | 5000.500 digits | 0.000333 | 0.000125 | 2.66× | 5 |
| mul | 5000.500 digits | 0.01917 | 0.01321 | 1.45× | 5 |
| add | 20000.2000 digits | 0.001125 | 0.0005 | 2.25× | 3 |
| mul | 20000.2000 digits | 0.04367 | 0.09292 | 0.47× | 3 |
| div | 100.10 / 10 dp | 0.000291 | 0.000083 | 3.51× | 7 |
| div | 1000.100 / 100 dp | 0.000916 | 0.001833 | 0.50× | 7 |
| div | 5000.500 / 500 dp | 0.009625 | 0.02558 | 0.38× | 5 |
| div | 20000.2000 / 2000 dp | 0.06133 | 0.1921 | 0.32× | 3 |
| div | 2000.200 / 0 dp | 0.0005 | 0.005 | 0.10× | 5 |
| div | 2000.200 / 10 dp | 0.00075 | 0.005083 | 0.15× | 5 |
| div | 2000.200 / 100 dp | 0.001375 | 0.005458 | 0.25× | 5 |
| div | 2000.200 / 1000 dp | 0.007375 | 0.009333 | 0.79× | 5 |
| div | 2000.200 / 5000 dp | 0.03 | 0.03367 | 0.89× | 5 |
| parse | 100 digits | 0.000375 | 0.000459 | 0.82× | 7 |
| parse | 1000 digits | 0.004625 | 0.004542 | 1.02× | 7 |
| parse | 10000 digits | 0.08612 | 0.1045 | 0.82× | 5 |
| parse | 50000 digits | 0.6236 | 0.9947 | 0.63× | 3 |
| tostr | 100 digits | 0.00025 | 0.0005 | 0.50× | 7 |
| tostr | 1000 digits | 0.006583 | 0.004583 | 1.44× | 7 |
| tostr | 10000 digits | 0.1764 | 0.1012 | 1.74× | 5 |
| tostr | 50000 digits | 1.289 | 1.084 | 1.19× | 3 |

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
