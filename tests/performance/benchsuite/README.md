# Split BigMath-vs-GMP benchmark suite

`bench_vs_gmp.cpp` regenerates random operands and reruns GMP on every
invocation, which is slow. This suite splits that into three cacheable stages so
repeated BigMath runs only pay for BigMath.

## Usage

```
./run_benchmark.sh [profile] [--force-dataset] [--force-gmp]
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
