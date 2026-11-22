# Split BigMath-vs-GMP benchmark suite

`bench_vs_gmp.cpp` regenerates random operands and reruns GMP on every
invocation, which is slow. This suite splits that into three cacheable stages so
repeated BigMath runs only pay for BigMath.

## Usage

```
./run_benchmark.sh [profile] [--force-dataset] [--force-gmp] [--memory | --extended]
```

`profile` = `quick` | `default` | `full` (default `default`). The driver runs
CMake on every invocation so library and workload dependencies are tracked.
Set `HISTORY_MD` to direct local test runs to a separate history file.

## Stages

1. **Dataset generation** (`gen_dataset`) — materializes random operands once
   from deterministic per-entry seeds, writes `manifest.csv` + one decimal file
   per operand, and the driver zips it to `dataset_<profile>.zip`. Skipped if the
   zip or an extracted copy already exists.
2. **GMP run** (`run_gmp`) — times GMP on the dataset and records a
   base-independent value hash for each checked output. The driver reuses
   `gmp_v2_<machine>_<profile>.csv` until its runner or manifest changes.
3. **BigMath run** (`run_bigmath`) — times BigMath on the same dataset, verifies
   mul/div against the cached GMP hashes, and prints BigMath/GMP ratios. Saved
   under an incrementing run id (`bigmath_<machine>_<profile>_<id>.csv`); the
   `historical`-flagged rows are appended to `HISTORY.md`.

## Layout

- Cache (gitignored): `$BENCH_DIR` (default `<repo>/.benchcache`)
  - `dataset/<profile>/` + `dataset_<profile>.zip` — operands
  - `results/` — per-machine GMP + per-run BigMath CSVs
  - `cmake/` — CMake build tree and stage binaries
- Committed: `HISTORY.md` — long-lived comparison table

## Notes

- Operands are reproducible on any machine (fixed seeds), so the dataset is
  machine-independent; only timings/results are machine-keyed.
- Each division entry supplies identical operands for `div_q` and `div_qr`.
  Both sides compute the requested outputs. QR records both value hashes.
  CSVs record one first call and the fastest of the following `warm_reps`
  calls. A first call is per operation and case within the running process.
- The value hash is FNV-1a over the integer's little-endian 32-bit word stream
  (one word per 32-bit limb or two per 64-bit limb, trimmed to match GMP export), so
  correctness is checked without the slow decimal `ToString` path.
- Workload tables live in `workload.h`; the on-disk `manifest.csv` is the source
  of truth at run time.
- CMake also exposes `bench_gen_dataset`, `bench_run_gmp`, `bench_run_bigmath`
  for IDE/CI parity; the shell driver is the primary interface.

## Isolated peak memory and time

On macOS or Linux, the memory mode runs GMP, the default BigMath build and a
second BigMath build compiled with `-DBIGMATH_USE_THREADS=0`. Use a threaded
default build for the first comparison. The coordinator requires Python 3.9
or later. Run either command from the repository root.

```sh
bash tests/performance/benchsuite/run_benchmark.sh quick --memory
bash tests/performance/benchsuite/run_benchmark.sh full --extended
```

Each library and operation gets a fresh executable process. The CSV records
the operating system's peak resident set size (RSS) in bytes
after input setup, after the first call, and after all warm calls. Each peak
includes earlier phases, the resident inputs, output, allocator and thread
workspaces. Result hashing follows the last reading. Subtracting the setup peak gives an increase in the high-water mark. For
small inputs the runtime itself can account for most of the total. The
`warm_ms` column is the fastest warm call, with repetitions from the manifest.

Integer arithmetic and formatting use the original decimal fixtures converted
to GMP raw integers in a separate preparation process. Both libraries read the
same bytes directly into integer storage. Parsing still starts from decimal
text. The full sweep includes the existing BigDecimal cases against GMP's
`mpf`; those operations have different decimal/rounding semantics and are
timed without an equality check. Every completed integer multiplication and
division is checked against the GMP result hash.

`--extended` adds balanced products around transform boundaries and at 500M,
750M and 1B decimal digits **per operand**, including 484M and 486M on opposite sides of the short-transform boundary. These extra fixtures use seeded
binary magnitudes with exactly the specified decimal length, generated outside
the measurement process. Balanced inputs between about 323M and 646M digits use two convolution residues in the default configuration. Larger inputs use block products.
Historical executables that rejected those sizes reported `crt_length_limit` with empty timing and memory fields.
Other exceptions, hash mismatches and failed processes fail the run. Raw inputs
are discarded after each case. The driver prints the unique CSV path under
`$BENCH_DIR/results/`.

Use the Python driver directly with existing builds and an extracted dataset.

```sh
python3 tests/performance/benchsuite/run_memory.py DATASET DEFAULT_BUILD SINGLE_BUILD results.csv --extended
python3 tests/performance/benchsuite/run_memory.py DATASET DEFAULT_BUILD SINGLE_BUILD boundary.csv --balanced-digits 484000000 486000000 500000000
python3 docs/images/make_memory_plots.py --csv results.csv
```
