# Why decimal conversion uses more memory than GMP

This investigation measured BigMath at `c861c9f` and GMP 6.3.0 on an Apple M1 Max with 32 GiB of RAM. Both libraries began with essentially the same input footprint. BigMath’s three-prime number-theoretic transforms (NTTs) accounted for most of the large allocations still live at the peak. Formatting also keeps powers of ten and Newton divisors for recursive quotient-and-remainder splits. [Measured revision](https://github.com/mmurshed/biginteger/tree/c861c9f) · [Conversion implementation](../../src/common/Parser.cpp)

## Peak RSS and warm time

Each cell reports the median process peak through the warm call and the median warm time across three fresh processes. Each process makes one first and one warm conversion on the same operand. Both BigMath configurations and GMP use the full benchmark’s input files. The serial-subtree control keeps threaded arithmetic available and turns off parallel conversion subtrees. [Individual measurements](2026-11-21-conversion-memory-controls.csv)

| Case | GMP, 1 thread | BigMath, 1 thread | BigMath, up to 8 threads | BigMath, serial subtrees |
| --- | ---: | ---: | ---: | ---: |
| Parse 20M | 111 MiB · 808 ms | 242 MiB · 912 ms | 324 MiB · 165 ms | 257 MiB · 418 ms |
| Parse 50M | 305 MiB · 2,670 ms | 553 MiB · 2,717 ms | 1,080 MiB · 542 ms | 873 MiB · 1,206 ms |
| Format 10M | 48 MiB · 911 ms | 288 MiB · 1,123 ms | 375 MiB · 240 ms | 309 MiB · 467 ms |
| Format 20M | 102 MiB · 2,192 ms | 520 MiB · 2,547 ms | 756 MiB · 515 ms | 634 MiB · 967 ms |

## Where the allocations go

The following figures count simultaneously live allocation requests of at least 32 KiB in one additional process per case. The input and output column includes the input, the result string and worker substring buffers when live. The workspace column includes writable coefficient planes and transpose buffers. Immutable tables include roots and precomputed multipliers. The remaining column includes powers of ten, divisors, reciprocals, limb results and other arithmetic temporaries. The last column is the live total immediately after the warm result was destroyed, with the input still present. Units are MiB. [Allocation totals for all three builds](2026-11-21-conversion-memory-allocations.csv)

| Threaded case | Workspace | Tables | Input and output | Other | Total at peak | After return |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Parse 20M | 140.7 | 44.6 | 19.1 | 42.1 | 246.5 | 218.7 |
| Parse 50M | 447.0 | 170.6 | 47.7 | 107.0 | 772.4 | 703.1 |
| Format 10M | 150.4 | 80.6 | 22.2 | 49.6 | 302.8 | 259.3 |
| Format 20M | 288.0 | 161.6 | 44.5 | 98.7 | 592.8 | 507.2 |

- **Parsing.** Combining the two halves multiplies the upper half by a power of ten. Between 20M and 50M digits, the top product’s padded transform length grows from 2²¹ to 2²³. At 50M, the caller allocated nine 32 MiB planes for this product alone. At the warm peak, the caller and workers accounted for 305 and 142 MiB of writable transform storage respectively. The input string was 47.7 MiB. [Parser joins](../../src/common/Parser.cpp) · [Fused transform buffers](../../include/biginteger/algorithms/multiplication/NTTMultiplicationCrt.h)
- **Formatting.** Splitting by powers of ten requires divisions. Newton division multiplies by cached reciprocals, computes remainders, and uses both linear and cyclic products. At the same 20M digit count, formatting used 288 MiB of transform workspace versus 141 MiB for parsing, and 162 MiB of transform tables versus 45 MiB. Formatting’s entire input and output buffer category was 44.5 MiB. [Power-of-ten divider chain](../../src/common/Parser.cpp) · [Newton division](../../include/biginteger/algorithms/division/NewtonDivision.h)
- **Between calls.** The thread pool retains writable buffers for reuse. Workers share immutable plans by size and prime, and the parser also caches powers of ten by exponent. In the 50M parse and 20M format cases, GMP’s tracked live allocations fell back to the input alone, 47.7 and 7.9 MiB respectively. BigMath retained 703 and 507 MiB including input. With one BigMath thread the corresponding totals were 339 and 347 MiB. These measurements cover repeated calls at one size.

Physical peak RSS also includes smaller allocations, allocator-managed pages, stacks and runtime memory. The allocation and RSS peaks can occur at different times. The allocation breakdown therefore explains ownership without assigning every RSS byte to a source line.

## Two independent build experiments

The first experiment moved the fused multiplication scratch-reuse cutoff from 2²⁴ to 2²⁰ transform coefficients. Processing the three primes in sequence allows three output planes to share two scratch planes. The second kept the stock scratch cutoff and lowered `BIGMATH_NTT_TRUNC_MAX` from 2²⁰ to 2¹⁷, limiting the range of truncated products that use dense root tables. Each experiment has its own interleaved stock reference. [Raw first and warm observations](2026-11-21-conversion-memory-controls.csv)

| Case | Stock | Earlier five-plane reuse | Stock, second run | Lower dense truncation limit |
| --- | ---: | ---: | ---: | ---: |
| Parse 20M | 324 MiB · 165 ms | 286 MiB · 167 ms | 325 MiB · 172 ms | 323 MiB · 170 ms |
| Parse 50M | 1,080 MiB · 542 ms | 891 MiB · 556 ms | 1,078 MiB · 577 ms | 1,210 MiB · 592 ms |
| Format 10M | 375 MiB · 240 ms | 329 MiB · 247 ms | 372 MiB · 236 ms | 331 MiB · 235 ms |
| Format 20M | 756 MiB · 515 ms | 649 MiB · 515 ms | 746 MiB · 512 ms | 610 MiB · 509 ms |

Earlier scratch reuse lowered peak RSS in all four conversions. Parsing 50M digits used 17.5% less memory with a 2.5% increase in median warm time. Formatting 20M digits used 14.2% less memory at essentially the same warm time. Limiting dense truncation improved formatting 20M by 18.2% in peak RSS, but increased the 50M parse peak by 12.2%.

The five-plane experiment also affected ordinary multiplication. A 5M-by-5M product rose from 7.66 to 9.22 ms (20.4%), and a 40M-by-40M product rose from 82.19 to 87.13 ms (6.0%). These results motivated the conversion-specific five-plane path. Both the caller and workers apply it while converting decimal values. [Implemented change and paired results](../../BENCHMARK.md#reuse-two-scratch-buffers-during-decimal-conversion)

| Digits per operand | Stock | Earlier five-plane reuse |
| --- | ---: | ---: |
| 5M | 64 MiB · 8 ms | 54 MiB · 9 ms |
| 10M | 120 MiB · 18 ms | 98 MiB · 17 ms |
| 20M | 224 MiB · 38 ms | 171 MiB · 38 ms |
| 40M | 435 MiB · 82 ms | 326 MiB · 87 ms |

These standalone controls used matched operand pairs, one warm product per process and three processes per build. Cases 8–10 came from the full manifest, and the 40M case used the extended sweep’s seeded binary operands. All products matched GMP’s limb hash.

## Reproduction

Use Release builds with `-O3 -march=native` and the repository’s `bench_run_bigmath` and `bench_run_gmp` executables. Select cases 53, 54, 64 and 65 of the full manifest, in that order. Keep the original parse text files. For formatting, use `bench_run_gmp --export-integer FULL_DATASET CASE_ID 0 RAW_PATH` to prepare each binary fixture before launching the measured process. Invoke `RUNNER --memory-case FIXTURE_DIR INDEX` for indices 0 through 3, three times per build with varied build order. The serial-subtree control sets `BIGMATH_PARSE_PARALLEL_THRESHOLD` and `BIGMATH_TOSTR_PARALLEL_THRESHOLD` to 2147483647. The one-thread build uses `BIGMATH_USE_THREADS=0`.

The separate allocation probes recorded allocation and deallocation sizes via C++ new/delete and GMP’s memory hooks. The BigMath probe symbolized allocation stacks and took snapshots of live bytes by allocation site during the warm call. Reported wall times and RSS come from uninstrumented builds. The two experimental builds passed multiplication correctness. Across the four conversion fixtures, parsed limb values matched GMP exactly, and formatted strings matched the original decimal operands byte-for-byte.
