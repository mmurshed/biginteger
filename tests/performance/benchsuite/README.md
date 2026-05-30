# Split BigMath-vs-GMP benchmark suite

`bench_vs_gmp.cpp` regenerates random operands and reruns GMP on every
invocation, which is slow. This suite splits that into three cacheable stages so
repeated BigMath runs only pay for BigMath.

## Usage

```
./run_benchmark.sh [profile] [--force-dataset] [--force-gmp]
```

`profile` = `quick` | `default` | `full` (default `default`). The driver runs
only the stages that are stale.

## Stages

1. **Dataset generation** (`gen_dataset`) — materializes random operands once
   from deterministic per-entry seeds, writes `manifest.csv` + one decimal file
   per operand, and the driver zips it to `dataset_<profile>.zip`. Skipped if the
   zip or an extracted copy already exists.
2. **GMP run** (`run_gmp`) — times GMP on the dataset and records a
   base-independent value hash for the exact-integer ops. Saved to
   `gmp_<machine>_<profile>.csv` and reused across BigMath edits, since GMP
   timings don't change. Runs once per machine per profile.
3. **BigMath run** (`run_bigmath`) — times BigMath on the same dataset, verifies
   mul/div against the cached GMP hashes, and prints BigMath/GMP ratios. Saved
   under an incrementing run id (`bigmath_<machine>_<profile>_<id>.csv`); the
   `historical`-flagged rows are appended to `HISTORY.md`.

## Layout

- Cache (gitignored): `$BENCH_DIR` (default `<repo>/.benchcache`)
  - `dataset/<profile>/` + `dataset_<profile>.zip` — operands
  - `results/` — per-machine GMP + per-run BigMath CSVs
  - `bin/` — compiled stage binaries
- Committed: `HISTORY.md` — long-lived comparison table

## Notes

- Operands are reproducible on any machine (fixed seeds), so the dataset is
  machine-independent; only timings/results are machine-keyed.
- The value hash is FNV-1a over the integer's little-endian 32-bit word stream
  (low 32 bits of each base-2^32 BigMath limb; `mpz_export` for GMP), so
  correctness is checked without the slow decimal `ToString` path.
- Workload tables live in `workload.h`; the on-disk `manifest.csv` is the source
  of truth at run time.
- CMake also exposes `bench_gen_dataset`, `bench_run_gmp`, `bench_run_bigmath`
  for IDE/CI parity; the shell driver is the primary interface.
