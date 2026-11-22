# BigMath vs GMP benchmark results

## BigInteger equality

This benchmark compares the previous implementation of `a == b` (`a.CompareTo(b) == 0`) with the specialized equality operator in the same Release binary. The Apple M1 Max ran Apple Clang 21 with the default 64-bit limbs. Each timing is the median of nine warm batches on the same two operands. The harness calibrates repetitions separately for each predicate, alternates their measurement order, and checks the result of every batch. Construction and copying are outside the timer. Nanosecond-scale timings include an indirect function call. [CSV for 1, 4, 1,000 and 1,000,000 limbs](docs/benchmarks/2026-11-21-equality-m1-max.csv)

| Comparison at 1,000,000 limbs | Previous equality | Specialized equality |
|---|---:|---:|
| Same object | 324 µs | 1.58 ns |
| Equal copy | 337 µs | 247 µs |
| Lowest limb differs | 336 µs | 2.19 ns |
| Middle limb differs | 163 µs | 99.8 µs |
| Second-highest limb differs | 3.14 ns | 4.72 ns |
| Highest limb differs | 2.73 ns | 1.96 ns |
| Opposite signs | 1.24 ns | 1.55 ns |

Run the benchmark from a configured Release build.

```sh
cmake --build build --target equality_bench
./build/equality_bench
```

## Single-thread comparison

The November 21, 2026 full timing run measured revision `c945bdc` on the Apple M1 Max. At 200M digits per operand the single-thread build took 2.92 s and GMP took 2.87 s. At 50M and 100M digits BigMath used 0.98× and 1.00× GMP’s time. [Measured revision](https://github.com/mmurshed/biginteger/tree/c945bdc)

The full CSV includes the first call, fastest subsequent call, warm repetition count and arithmetic check status for all 111 rows. Most large sizes have one warm observation per library. All 24 product checks and 40 division checks matched GMP’s reference hashes. Each division mode was checked against the corresponding GMP operation. [README performance table](README.md#performance-numbers) · [Single-thread CSV](docs/benchmarks/2026-11-21-crt-ceiling-full-m1-max-single-thread.csv) · [Single-thread figure](docs/images/bigmath_vs_gmp_single_thread.png)

The timing figures also show balanced products from 200M to 1B digits measured in fresh processes on the same revision. The two 200M samples use the same operands and show how the long-lived and fresh-process runs compare. [Isolated timings, peaks and hashes](docs/benchmarks/2026-11-21-crt-ceiling-full-m1-max-memory.csv)

To repeat the single-thread run, give the driver its own build and results directory. The figure generator uses the checked-in CSVs and requires Matplotlib.

```sh
CXXFLAGS=-DBIGMATH_USE_THREADS=0 BENCH_DIR="$PWD/.benchcache/single-thread" \
  HISTORY_MD="$PWD/.benchcache/single-thread/history.md" \
  bash tests/performance/benchsuite/run_benchmark.sh full --force-gmp
python3 docs/images/make_benchmark_plots.py --ratio-only
```

## Default build with up to eight BigMath threads

With up to eight BigMath threads, equal-length products at 50M, 100M and 200M digits took 182.6, 346.1 and 720.8 ms. The 200M-digit product used 0.25× GMP’s time. [Default-build CSV](docs/benchmarks/2026-11-21-crt-ceiling-full-m1-max.csv) · [GMP reference and hashes](docs/benchmarks/2026-11-21-crt-ceiling-full-m1-max-gmp.csv) · [Previous full run](docs/benchmarks/2026-11-21-blocked-crt-m1-max.csv)

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

To measure time and peak RSS in fresh processes, including the 500M, 750M and 1B products, run the extended profile.

```sh
bash tests/performance/benchsuite/run_benchmark.sh full --extended
```

To regenerate all figures from their published inputs, run these commands from the repository root. The historical before-and-after figures use the measurements from their respective experiments.

```sh
python3 docs/images/make_benchmark_plots.py
python3 docs/images/make_memory_plots.py
python3 docs/images/make_memory_fix_plot.py
python3 docs/images/make_division_balanced_plot.py
python3 docs/images/make_mfa_fusion_plot.py
python3 docs/algorithms-guide/figures/make_figures.py
(cd docs/algorithms-guide && latexmk -pdf algorithms-guide.tex)
```

A separate run of `./build/bench_vs_gmp` on PR #152 measured both libraries in one process on different seeded operands. It completed 151 measurements. Its single 200-million-digit balanced product took 3,721 ms in BigMath and 3,057 ms in GMP (1.22×). [In-process output](docs/benchmarks/2026-11-21-pr152-in-process.txt)

### Multiplication

The equal-length series measures two independently generated operands. The skewed series holds their digit-length ratio at 10:1.

#### Equal lengths

| Operand digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 × 1,000 | 0.07617 | 0.009125 | 8.35× | 0.001834 | 0.000916 | 2.00× | 7 |
| 5,000 × 5,000 | 0.03692 | 0.01983 | 1.86× | 0.024 | 0.01179 | 2.04× | 5 |
| 10,000 × 10,000 | 0.6702 | 0.03208 | 20.89× | 0.05008 | 0.029 | 1.73× | 5 |
| 50,000 × 50,000 | 0.6565 | 0.2757 | 2.38× | 0.2172 | 0.2712 | 0.80× | 3 |
| 100,000 × 100,000 | 1.155 | 0.6298 | 1.83× | 0.3307 | 0.6104 | 0.54× | 3 |
| 500,000 × 500,000 | 2.48 | 4.388 | 0.57× | 1.399 | 4.298 | 0.33× | 1 |
| 1,000,000 × 1,000,000 | 5.016 | 9.21 | 0.54× | 2.897 | 9.11 | 0.32× | 1 |
| 2,000,000 × 2,000,000 | 11.28 | 22.8 | 0.49× | 7.559 | 22.36 | 0.34× | 1 |
| 5,000,000 × 5,000,000 | 11.19 | 63.99 | 0.17× | 7.24 | 62.12 | 0.12× | 1 |
| 10,000,000 × 10,000,000 | 27.13 | 209.3 | 0.13× | 18.83 | 209.5 | 0.09× | 1 |
| 20,000,000 × 20,000,000 | 51.94 | 274.7 | 0.19× | 40.17 | 267 | 0.15× | 1 |
| 50,000,000 × 50,000,000 | 222.7 | 675.4 | 0.33× | 182.6 | 658.5 | 0.28× | 1 |
| 100,000,000 × 100,000,000 | 425.6 | 1,402 | 0.30× | 346.1 | 1,394 | 0.25× | 1 |
| 200,000,000 × 200,000,000 | 944.4 | 2,909 | 0.32× | 720.8 | 2,867 | 0.25× | 1 |

#### Skewed lengths

| Operand digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 100,000 × 10,000 | 0.1501 | 0.3344 | 0.45× | 0.1098 | 0.328 | 0.33× | 3 |
| 500,000 × 50,000 | 0.4857 | 2.232 | 0.22× | 0.4454 | 2.204 | 0.20× | 1 |
| 1,000,000 × 100,000 | 1.307 | 4.627 | 0.28× | 1.086 | 4.558 | 0.24× | 1 |
| 2,000,000 × 200,000 | 3.028 | 9.505 | 0.32× | 2.958 | 9.392 | 0.31× | 1 |
| 5,000,000 × 500,000 | 8.241 | 31.65 | 0.26× | 7.439 | 31.56 | 0.24× | 1 |
| 10,000,000 × 1,000,000 | 21.24 | 75.58 | 0.28× | 14.99 | 75.79 | 0.20× | 1 |
| 20,000,000 × 2,000,000 | 40.43 | 159.1 | 0.25× | 35.76 | 158.5 | 0.23× | 1 |
| 50,000,000 × 5,000,000 | 81.72 | 687.2 | 0.12× | 79.81 | 667.6 | 0.12× | 1 |
| 100,000,000 × 10,000,000 | 178.1 | 1,007 | 0.18× | 171.2 | 994.3 | 0.17× | 1 |
| 200,000,000 × 20,000,000 | 370 | 1,911 | 0.19× | 367.3 | 1,906 | 0.19× | 1 |

### Division

Quotient-only uses BigMath `/` and GMP `mpz_tdiv_q`. Quotient and remainder uses BigMath `DivideAndRemainder` and GMP `mpz_tdiv_qr`. Both modes use the same dividend and divisor. Equal-digit operands can produce a quotient of zero or one, so those results depend on the sampled values.

#### Skewed dividends

| Dividend ÷ divisor digits | BigMath Q ms | GMP Q ms | Q ratio | BigMath QR ms | GMP QR ms | QR ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 40,000 ÷ 10,000 | 0.2889 | 0.1976 | 1.46× | 0.2887 | 0.2191 | 1.32× | 3 |
| 100,000 ÷ 10,000 | 0.8487 | 0.4295 | 1.98× | 0.8503 | 0.4389 | 1.94× | 3 |
| 200,000 ÷ 50,000 | 0.9872 | 1.525 | 0.65× | 0.9931 | 1.665 | 0.60× | 1 |
| 500,000 ÷ 100,000 | 2.614 | 4.234 | 0.62× | 2.383 | 4.553 | 0.52× | 1 |
| 1,000,000 ÷ 200,000 | 4.338 | 9.746 | 0.45× | 4.689 | 10.27 | 0.46× | 1 |
| 2,000,000 ÷ 500,000 | 7.988 | 22.24 | 0.36× | 8.272 | 25.85 | 0.32× | 1 |
| 5,000,000 ÷ 1,000,000 | 24.66 | 62.95 | 0.39× | 24.88 | 67.22 | 0.37× | 1 |
| 10,000,000 ÷ 2,000,000 | 53.57 | 143.3 | 0.37× | 54.9 | 151.4 | 0.36× | 1 |
| 20,000,000 ÷ 4,000,000 | 86.75 | 319.9 | 0.27× | 87.64 | 341 | 0.26× | 1 |
| 50,000,000 ÷ 10,000,000 | 153.7 | 1,225 | 0.13× | 163.1 | 1,414 | 0.12× | 1 |
| 100,000,000 ÷ 20,000,000 | 332.7 | 2,319 | 0.14× | 328.8 | 2,505 | 0.13× | 1 |
| 200,000,000 ÷ 40,000,000 | 692.1 | 4,195 | 0.16× | 716.2 | 4,440 | 0.16× | 1 |

At 200,000 ÷ 50,000 digits the current full run measured 0.65× GMP’s time for the quotient and 0.60× for quotient and remainder. The earlier in-process run on PR #152 measured 0.92× and 1.04× on different operands. Both suites took one warm measurement per library at that size. [In-process output](docs/benchmarks/2026-11-21-pr152-in-process.txt)

#### Equal digit counts

| Dividend ÷ divisor digits | BigMath Q ms | GMP Q ms | Q ratio | BigMath QR ms | GMP QR ms | QR ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 ÷ 1,000 | 0.000125 | <1 tick | — | 0.000166 | 0.000041 | 4.05× | 7 |
| 5,000 ÷ 5,000 | 0.000083 | <1 tick | — | 0.000417 | 0.000167 | 2.50× | 5 |
| 10,000 ÷ 10,000 | 0.000041 | <1 tick | — | 0.000209 | 0.000083 | 2.52× | 5 |
| 50,000 ÷ 50,000 | 0.000125 | 0.000042 | 2.98× | 0.003875 | 0.001666 | 2.33× | 3 |
| 100,000 ÷ 100,000 | 0.000125 | 0.000166 | 0.75× | 0.008083 | 0.003542 | 2.28× | 3 |
| 500,000 ÷ 500,000 | 0.000125 | 0.000166 | 0.75× | 0.009583 | 0.004333 | 2.21× | 1 |
| 1,000,000 ÷ 1,000,000 | 0.000083 | 0.000167 | 0.50× | 0.01588 | 0.008542 | 1.86× | 1 |
| 5,000,000 ÷ 5,000,000 | 0.000042 | 0.000209 | 0.20× | 0.1084 | 0.04358 | 2.49× | 1 |

### Decimal parse (string → BigInteger)

| Digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.0055 | 0.004375 | 1.26× | 0.002666 | 0.001625 | 1.64× | 7 |
| 10,000 | 0.07892 | 0.05237 | 1.51× | 0.06375 | 0.039 | 1.63× | 5 |
| 50,000 | 0.8451 | 0.4853 | 1.74× | 0.5142 | 0.4011 | 1.28× | 3 |
| 100,000 | 1.113 | 1.116 | 1.00× | 0.548 | 1.062 | 0.52× | 3 |
| 500,000 | 3.018 | 8.816 | 0.34× | 2.96 | 8.714 | 0.34× | 1 |
| 1,000,000 | 6.644 | 21 | 0.32× | 6.869 | 20.39 | 0.34× | 1 |
| 2,000,000 | 17.39 | 47.68 | 0.36× | 14.82 | 46.77 | 0.32× | 1 |
| 5,000,000 | 38.03 | 151.6 | 0.25× | 35.84 | 148 | 0.24× | 1 |
| 10,000,000 | 78.87 | 345.9 | 0.23× | 76.48 | 344.5 | 0.22× | 1 |
| 20,000,000 | 162.1 | 808.4 | 0.20× | 166.2 | 805.6 | 0.21× | 1 |
| 50,000,000 | 528.3 | 2,644 | 0.20× | 519.1 | 2,663 | 0.19× | 1 |

### ToString (BigInteger → string)

| Digits | BigMath first ms | GMP first ms | First ratio | BigMath warm ms | GMP warm ms | Warm ratio | Warm reps |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | 0.008334 | 0.006541 | 1.27× | 0.00625 | 0.003583 | 1.74× | 7 |
| 10,000 | 0.1562 | 0.1015 | 1.54× | 0.1118 | 0.07767 | 1.44× | 5 |
| 50,000 | 1.351 | 0.9076 | 1.49× | 0.99 | 0.8577 | 1.15× | 3 |
| 100,000 | 1.826 | 2.314 | 0.79× | 1.122 | 2.296 | 0.49× | 1 |
| 200,000 | 3.816 | 6.042 | 0.63× | 2.322 | 6.027 | 0.39× | 1 |
| 500,000 | 9.13 | 20.52 | 0.44× | 5.933 | 20.5 | 0.29× | 1 |
| 1,000,000 | 20.52 | 49.72 | 0.41× | 13.21 | 49.66 | 0.27× | 1 |
| 2,000,000 | 39.67 | 119.7 | 0.33× | 28.12 | 119.5 | 0.24× | 1 |
| 5,000,000 | 147.2 | 381.9 | 0.39× | 107.2 | 381.1 | 0.28× | 1 |
| 10,000,000 | 314.9 | 901.3 | 0.35× | 233.9 | 900.8 | 0.26× | 1 |
| 20,000,000 | 668 | 2,127 | 0.31× | 498.7 | 2,120 | 0.24× | 1 |

### BigDecimal

BigMath uses decimal fixed-point arithmetic and HALF_EVEN rounding for division. The comparison uses GMP binary floating-point (`mpf_t`), with precision derived from the digit count and requested scale. The split runner times these operations without an equality check. A size of `1,000.100` denotes 1,000 integer digits and 100 fractional digits.

| Operation | Size | BigMath warm ms | GMP warm ms | Ratio | Warm reps |
|---|---:|---:|---:|---:|---:|
| add | 100.10 digits | 0.000083 | <1 tick | — | 7 |
| mul | 100.10 digits | 0.000042 | 0.000041 | 1.02× | 7 |
| add | 1000.100 digits | 0.000125 | <1 tick | — | 7 |
| mul | 1000.100 digits | 0.0015 | 0.001083 | 1.39× | 7 |
| add | 5000.500 digits | 0.000333 | 0.000125 | 2.66× | 5 |
| mul | 5000.500 digits | 0.01892 | 0.01279 | 1.48× | 5 |
| add | 20000.2000 digits | 0.001125 | 0.000458 | 2.46× | 3 |
| mul | 20000.2000 digits | 0.04562 | 0.09671 | 0.47× | 3 |
| div | 100.10 / 10 dp | 0.000291 | 0.000125 | 2.33× | 7 |
| div | 1000.100 / 100 dp | 0.000875 | 0.001916 | 0.46× | 7 |
| div | 5000.500 / 500 dp | 0.009625 | 0.02621 | 0.37× | 5 |
| div | 20000.2000 / 2000 dp | 0.05987 | 0.1993 | 0.30× | 3 |
| div | 2000.200 / 0 dp | 0.0005 | 0.005041 | 0.10× | 5 |
| div | 2000.200 / 10 dp | 0.000708 | 0.005042 | 0.14× | 5 |
| div | 2000.200 / 100 dp | 0.001375 | 0.005625 | 0.24× | 5 |
| div | 2000.200 / 1000 dp | 0.007375 | 0.009625 | 0.77× | 5 |
| div | 2000.200 / 5000 dp | 0.03008 | 0.03483 | 0.86× | 5 |
| parse | 100 digits | 0.000375 | 0.0005 | 0.75× | 7 |
| parse | 1000 digits | 0.004583 | 0.004583 | 1.00× | 7 |
| parse | 10000 digits | 0.08496 | 0.1046 | 0.81× | 5 |
| parse | 50000 digits | 0.61 | 1.024 | 0.60× | 3 |
| tostr | 100 digits | 0.000292 | 0.000458 | 0.64× | 7 |
| tostr | 1000 digits | 0.006625 | 0.00475 | 1.39× | 7 |
| tostr | 10000 digits | 0.1124 | 0.1057 | 1.06× | 5 |
| tostr | 50000 digits | 0.9761 | 1.095 | 0.89× | 3 |

## Memory optimization benchmarks

Each comparison uses Release builds on the Apple M1 Max, identical operands and a fresh process for each measurement. Tables give peak RSS through the warm call and warm time. BigMath cells are medians of three processes per build, reversing the run order for the second pair. GMP 6.3.0 ran once per shape. Each process made one first call and one warm call. Products were hashed against GMP.

### Threaded skew dispatch

Compared with `9a3f89a`, skewed products with an available fused matrix transform now bypass truncation at transform lengths of at least 2²¹ and a limb ratio of at least eight. The 5M case exercises the original path. Standard shapes use the full benchmark manifest. Intermediate shapes use the seeded binary generator from `run_memory.py`. [All observations](docs/benchmarks/memory-fix-skew-default.csv)

| Shape | Before | After | GMP |
|---|---:|---:|---:|
| mul 5M × 500k | 87.9 MiB · 8.4 ms | 87.9 MiB · 8.7 ms | 12.2 MiB · 31.8 ms |
| mul 10M × 1.24M | 173.6 MiB · 19.0 ms | 110.8 MiB · 15.8 ms | 23.6 MiB · 69.8 ms |
| mul 10M × 1M | 173.3 MiB · 18.0 ms | 109.9 MiB · 15.8 ms | 18.5 MiB · 69.3 ms |
| mul 10M × 625k | 169.1 MiB · 18.0 ms | 110.5 MiB · 16.7 ms | 17.3 MiB · 78.3 ms |
| mul 12M × 1.2M | 180.4 MiB · 24.9 ms | 112.2 MiB · 17.7 ms | 26.3 MiB · 65.0 ms |

### Compact single-thread weights

The large serial truncated transform now factors each weight into a block and offset. The baseline is the single-thread path from `9a3f89a`, also used by `dbe3b06`. The 10M case exercises the smaller transform. The 162M inputs use the seeded generator. [All observations](docs/benchmarks/memory-fix-weights-single.csv)

| Shape | Before | After | GMP |
|---|---:|---:|---:|
| mul 10M × 10M | 107.4 MiB · 110.9 ms | 107.4 MiB · 111.1 ms | 75.4 MiB · 207.4 ms |
| mul 50M × 50M | 919.9 MiB · 632.7 ms | 345.2 MiB · 640.7 ms | 280.6 MiB · 653.8 ms |
| mul 162M × 162M | 3,410.6 MiB · 2,400.2 ms | 1,112.5 MiB · 2,419.3 ms | 932.3 MiB · 2,369.8 ms |
| mul 200M × 200M | 3,640.6 MiB · 2,860.5 ms | 1,342.5 MiB · 2,917.0 ms | 1,116.8 MiB · 2,872.3 ms |

### Reuse threaded transform buffers

Compared with `4e13510`, fused products of length 2²⁴ and above process each prime in parallel over rows, then reuse both scratch planes for the next prime. This cuts four full-size planes. The 160M and 162M cases straddle a transform doubling. [All observations](docs/benchmarks/memory-fix-planes-default.csv)

| Shape | Before | After | GMP |
|---|---:|---:|---:|
| mul 50M × 50M | 766.4 MiB · 168.4 ms | 533.7 MiB · 175.2 ms | 281.8 MiB · 675.0 ms |
| mul 160M × 160M | 1,658.7 MiB · 365.1 ms | 1,181.5 MiB · 377.7 ms | 895.5 MiB · 2,498.3 ms |
| mul 162M × 162M | 2,878.4 MiB · 693.5 ms | 1,898.3 MiB · 715.8 ms | 932.3 MiB · 2,409.0 ms |
| mul 200M × 200M | 2,966.6 MiB · 713.9 ms | 1,988.5 MiB · 720.4 ms | 1,117.4 MiB · 2,912.1 ms |

### Share conversion transform storage

The baseline is `c49f04c`. A separate instrumented run counted table capacities at construction and writable vector capacities at thread exit. Parsing 50M digits built 648 MiB of leaf plans across workers for 144 MiB of distinct tables. Formatting 10M digits built 198 MiB for 72 MiB of distinct tables and retained 33 MiB of cyclic-product buffers in addition to linear-product buffers. Capacity counts exclude allocator overhead. The timing runs use uninstrumented binaries. [Capacity measurements](docs/benchmarks/memory-conversion-retained-capacity.csv) [All observations](docs/benchmarks/memory-fix-conversion-default.csv)

| Shape | Before | After | GMP |
|---|---:|---:|---:|
| parse 1M | 34.2 MiB · 6.9 ms | 28.5 MiB · 7.9 ms | 8.1 MiB · 21.1 ms |
| parse 50M | 1,588.4 MiB · 566.5 ms | 1,084.4 MiB · 543.6 ms | 307.0 MiB · 2,664.6 ms |
| tostr 10M | 525.5 MiB · 242.0 ms | 372.2 MiB · 236.7 ms | 49.1 MiB · 914.5 ms |
| tostr 20M | 1,052.1 MiB · 525.9 ms | 743.8 MiB · 509.4 ms | 102.3 MiB · 2,142.8 ms |

## Peak memory and large-product time

The November 21, 2026 run measured revision `c945bdc` on an Apple M1 Max with 32 GiB of RAM. Each case ran in its own process for each library. The 111 standard cases share operands with the timing tables above. Seventeen additional balanced products use seeded inputs up to one billion decimal digits per operand. All 81 integer arithmetic results per BigMath build matched GMP. [Times, RSS and hashes](docs/benchmarks/2026-11-21-crt-ceiling-full-m1-max-memory.csv) · [Previous full sweep](docs/benchmarks/2026-11-21-blocked-crt-memory-m1-max.csv) · [Benchmark instructions](tests/performance/benchsuite/README.md#isolated-peak-memory-and-time)

**How to read the figures.** RSS is the operating system’s high-water mark for resident process memory. We record it with the inputs ready, after the first operation, and after the last warm call, before hashing. Arithmetic panels show the last high-water mark. The warm BigMath product assigns a new result to the previous output variable. GMP passes the same destination back to `mpz_mul`. The measurements include the storage each public operation uses. Conversion panels show both readings. A ratio of 4 means BigMath reached four times GMP’s peak. Small operands include roughly 2 MiB of process startup memory. A fresh process at every size excludes caches built by preceding cases. The ordinary timing tables use a long-lived benchmark process.

![Default BigMath versus GMP peak RSS across four operation families](docs/images/bigmath_vs_gmp_memory.png)

![Single-threaded BigMath versus GMP peak RSS across four operation families](docs/images/bigmath_vs_gmp_memory_single_thread.png)

### Memory spikes to investigate

- Multiplication at 5M × 500k digits had the largest threaded peak RSS relative to GMP in the integer sweep, 8.38× (88.2 versus 10.5 MiB). The single-thread maximum was 8.35× for multiplication at 10M × 1M digits (155.1 versus 18.6 MiB).
- Formatting 10M digits on the threaded build used 373.3 MiB versus 47.4 MiB for GMP. The 10M × 1M product used 111.8 versus 18.6 MiB. [Conversion controls](#share-conversion-transform-storage) · [Skewed multiplication controls](#threaded-skew-dispatch)
- From 160M to 162M digits per operand, the threaded product rose from 1.15 to 1.85 GiB as its transform length doubled. At 200M digits the threaded and single-thread builds used 1.94 and 1.31 GiB. GMP used 1.09 GiB.

### Controlled comparisons

These historical controls measured the implementation before the four memory optimizations. Each cell gives the median peak RSS through warm calls and the median warm time across three independent processes, using the same operands and reversing the build order for the second repeat. The transform control sets `BIGMATH_NTT_TRUNC=0`. The conversion control raises the separate parsing and formatting parallel thresholds above the tested sizes. Other build flags match the corresponding stock build. To reproduce, select entries 11, 19, 54 and 64 of the full manifest and run the isolated benchmark on the stock and control builds. [Control measurements](docs/benchmarks/2026-11-21-memory-controls-m1-max.csv) · [CRT plans and workspaces](include/biginteger/algorithms/multiplication/NTTMultiplicationCrt.h) · [Conversion implementation](src/common/Parser.cpp)

| Operation | Build | Stock | Control |
|---|---|---:|---:|
| 10M × 1M multiply | Up to 8 threads | 173 MiB · 20.0 ms | 110 MiB · 15.6 ms |
| 50M × 50M multiply | 1 thread | 920 MiB · 624.4 ms | 735 MiB · 1,055.0 ms |
| Format 10M digits | Up to 8 threads | 527 MiB · 234.6 ms | 324 MiB · 439.0 ms |
| Parse 50M digits | Up to 8 threads | 1,595 MiB · 554.8 ms | 874 MiB · 1,190.0 ms |

The controls motivated the threaded dispatch change, compact serial weights and conversion cache sharing. [Before-and-after measurements](#memory-optimization-benchmarks)

### Balanced products through one billion digits

Each cell reports one warm product time and the process peak RSS through that call. The 17 additional sizes use seeded binary operands, and the standard sizes use operands from the full manifest. All three libraries ran on matching inputs. [Full isolated measurements](docs/benchmarks/2026-11-21-crt-ceiling-full-m1-max-memory.csv)

| Digits per operand | BigMath, up to 8 threads | BigMath, 1 thread | GMP, 1 thread |
|---:|---:|---:|---:|
| 10M | 0.017 s · 0.12 GiB | 0.109 s · 0.10 GiB | 0.205 s · 0.07 GiB |
| 10.2M | 0.035 s · 0.19 GiB | 0.110 s · 0.08 GiB | 0.209 s · 0.07 GiB |
| 20M | 0.042 s · 0.22 GiB | 0.222 s · 0.13 GiB | 0.278 s · 0.11 GiB |
| 20.3M | 0.080 s · 0.38 GiB | 0.250 s · 0.14 GiB | 0.342 s · 0.13 GiB |
| 40M | 0.083 s · 0.43 GiB | 0.487 s · 0.25 GiB | 0.554 s · 0.22 GiB |
| 41M | 0.168 s · 0.50 GiB | 0.544 s · 0.29 GiB | 0.746 s · 0.27 GiB |
| 50M | 0.171 s · 0.52 GiB | 0.658 s · 0.34 GiB | 0.672 s · 0.28 GiB |
| 80M | 0.193 s · 0.59 GiB | 1.074 s · 0.49 GiB | 1.221 s · 0.44 GiB |
| 82M | 0.352 s · 0.97 GiB | 1.192 s · 0.56 GiB | 1.450 s · 0.51 GiB |
| 100M | 0.399 s · 1.01 GiB | 1.412 s · 0.66 GiB | 1.404 s · 0.55 GiB |
| 160M | 0.379 s · 1.15 GiB | 2.328 s · 0.98 GiB | 2.557 s · 0.87 GiB |
| 162M | 0.720 s · 1.85 GiB | 2.475 s · 1.09 GiB | 2.421 s · 0.91 GiB |
| 200M | 0.706 s · 1.94 GiB | 2.906 s · 1.31 GiB | 2.906 s · 1.09 GiB |
| 320M | 0.778 s · 2.22 GiB | 4.856 s · 1.95 GiB | 5.076 s · 1.75 GiB |
| 324M | 1.119 s · 2.24 GiB | 6.527 s · 2.13 GiB | 5.078 s · 1.82 GiB |
| 480M | 1.392 s · 3.03 GiB | 8.232 s · 2.91 GiB | 8.202 s · 2.62 GiB |
| 484M | 1.383 s · 3.05 GiB | 8.368 s · 2.93 GiB | 8.282 s · 2.63 GiB |
| 486M | 1.378 s · 3.06 GiB | 8.287 s · 2.94 GiB | 8.278 s · 2.73 GiB |
| 500M | 1.502 s · 3.14 GiB | 8.991 s · 3.02 GiB | 8.640 s · 2.73 GiB |
| 520M | 1.683 s · 3.25 GiB | 9.419 s · 3.13 GiB | 8.865 s · 2.84 GiB |
| 750M | 6.990 s · 3.61 GiB | 43.454 s · 3.33 GiB | 12.695 s · 4.12 GiB |
| 1,000M | 12.290 s · 4.19 GiB | 75.265 s · 3.91 GiB | 18.578 s · 5.45 GiB |

The rerun crossed the old 485M timing boundary smoothly. From 484M to 486M digits, warm time changed from 1.383 to 1.378 s threaded and from 8.368 to 8.287 s single-threaded. At 750M and 1B, both BigMath builds used less peak memory than GMP.

**Why the old curve jumps.** A single CRT transform holds 2²⁶ coefficients, enough for 323.23M decimal digits per balanced operand. Just above that limit, the old dispatcher split each operand in half and ran four full block products. Around 484.84M, the single-thread build additionally doubled the shorter convolution inside each of those four products. From 484M to 486M the previous single-thread time rose from 14.32 to 18.93 s. The threaded build used a full matrix transform for each block on both sides and took 2.95 and 3.08 s.

**What changed.** In the default configurations, balanced inputs in the next octave use two convolution residues. A short product of the low input coefficients resolves the 485M boundary. At 500M the three lengths are 2²⁶, 2²⁵ and 2²³ per prime. In the earlier three-process comparison, median warm time fell from 3.13 to 1.50 s threaded and from 18.84 to 8.85 s single-threaded. Threaded peak RSS went from 3.03 to 3.14 GiB and single-threaded peak RSS from 2.75 to 3.02 GiB. Above about 646M digits the dispatcher returns to blocked multiplication. It computes nine block products at 750M and sixteen at 1B. [Paired measurements before and after PR #175](docs/benchmarks/2026-11-21-crt-ceiling-m1-max.csv) · [Convolution and reconstruction](docs/MULTIPLICATION.md#balanced-products-above-the-transform-limit)

Reproduce the focused run with the isolated benchmark. Run each build pair three times at 500M for the median.

```sh
python3 tests/performance/benchsuite/run_memory.py DATASET DEFAULT_BUILD SINGLE_BUILD boundary.csv \
  --balanced-digits 320000000 324000000 480000000 484000000 486000000 500000000 520000000
```

![Previous and revised warm multiplication times for both thread modes across the 323M and 485M boundaries](docs/images/bigmath_crt_ceiling.png)

![Measured warm product time and peak memory from ten million to one billion digits per operand](docs/images/bigmath_vs_gmp_large_products.png)

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
bit-exact.* [Original measurements for PR #78](https://github.com/mmurshed/biginteger/blob/d6752bc15f68a498eb7c951f3734eeceaf3ba75e/BENCHMARK.md#mfa-transpose-fusion-pr-78-2026-05-31)

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
