#!/usr/bin/env bash
# Driver for the split BigMath-vs-GMP benchmark.
#
#   ./run_benchmark.sh [profile] [--force-dataset] [--force-gmp] [--memory | --extended]
#
#   profile = quick | default | full   (default: "default")
#
# It runs only the stages that are stale:
#   1. dataset   — generated + zipped once per profile (skipped if the zip or an
#                  extracted copy already exists)
#   2. GMP run   — run once per machine per profile (skipped if its result CSV
#                  exists), since GMP timings don't change between BigMath edits
#   3. BigMath   — always run; saved under an incrementing run id and the
#                  flagged rows appended to HISTORY.md
#
# Cache lives in $BENCH_DIR (default <repo>/.benchcache); HISTORY.md is committed.
set -euo pipefail

# ── locations ────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/../../.." && pwd)"
PROFILE="default"
FORCE_DATASET=0
FORCE_GMP=0
MEMORY=0
EXTENDED=0
for arg in "$@"; do
    case "$arg" in
        quick|default|full) PROFILE="$arg" ;;
        --memory)           MEMORY=1 ;;
        --extended)         MEMORY=1; EXTENDED=1 ;;
        --force-dataset)    FORCE_DATASET=1 ;;
        --force-gmp)        FORCE_GMP=1 ;;
        *) echo "unknown arg: $arg" >&2; exit 2 ;;
    esac
done

BENCH_DIR="${BENCH_DIR:-$REPO/.benchcache}"
BIN_DIR="$BENCH_DIR/cmake"
DATASET_ROOT="$BENCH_DIR/dataset"
DATA_DIR="$DATASET_ROOT/$PROFILE"
ZIP="$BENCH_DIR/dataset_$PROFILE.zip"
RESULTS_DIR="$BENCH_DIR/results"
HISTORY_MD="${HISTORY_MD:-$SCRIPT_DIR/HISTORY.md}"
mkdir -p "$BIN_DIR" "$DATASET_ROOT" "$RESULTS_DIR"

# ── machine identity ─────────────────────────────────────────────────────────
if [[ "$(uname)" == "Darwin" ]]; then
    CPU="$(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo unknown)"
else
    CPU="$(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2 || echo unknown)"
fi
MACHINE="$(printf '%s_%s' "$CPU" "$(uname -m)" | tr -cs 'A-Za-z0-9' '_' | sed 's/_*$//')"
echo "[driver] machine=$MACHINE profile=$PROFILE bench_dir=$BENCH_DIR"

# CMake tracks all library and workload dependencies, including generated
# build flags. Reconfigure on every invocation to pick up caller overrides.
cmake -S "$REPO" -B "$BIN_DIR" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER="${CXX:-c++}" -DCMAKE_CXX_FLAGS="${CXXFLAGS:-}" \
    -DBIGMATH_BUILD_TESTS=ON
cmake --build "$BIN_DIR" --target bench_gen_dataset bench_run_gmp bench_run_bigmath --parallel

# ── stage 1: dataset ─────────────────────────────────────────────────────────
if [[ $FORCE_DATASET -eq 1 ]]; then rm -rf "$DATA_DIR" "$ZIP"; fi
if [[ ! -f "$DATA_DIR/manifest.csv" ]]; then
    if [[ -f "$ZIP" ]]; then
        echo "[dataset] extracting cached $ZIP"
        unzip -q -o "$ZIP" -d "$DATASET_ROOT"
    else
        echo "[dataset] generating ($PROFILE)"
        mkdir -p "$DATA_DIR"
        "$BIN_DIR/bench_gen_dataset" "$DATA_DIR" "$PROFILE"
        echo "[dataset] zipping -> $ZIP"
        ( cd "$DATASET_ROOT" && zip -q -r -1 "$ZIP" "$PROFILE" )
    fi
else
    echo "[dataset] using existing $DATA_DIR"
fi

# RSS is a process high-water mark. Run each library/case in its own process
# and use the same isolated GMP run for both BigMath configurations.
if [[ $MEMORY -eq 1 ]]; then
    SINGLE_BIN="$BENCH_DIR/cmake-memory-single"
    cmake -S "$REPO" -B "$SINGLE_BIN" -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_COMPILER="${CXX:-c++}" \
        -DCMAKE_CXX_FLAGS="${CXXFLAGS:-} -DBIGMATH_USE_THREADS=0" -DBIGMATH_BUILD_TESTS=ON
    cmake --build "$SINGLE_BIN" --target bench_run_bigmath --parallel
    memory_args=()
    memory_profile="$PROFILE"
    if [[ $EXTENDED -eq 1 ]]; then
        memory_args+=(--extended)
        memory_profile="${PROFILE}_extended"
    fi
    out=$(mktemp "$RESULTS_DIR/memory_${MACHINE}_${memory_profile}_XXXXXX")
    trap 'rm -f "$out"' EXIT
    python3 "$SCRIPT_DIR/run_memory.py" "$DATA_DIR" "$BIN_DIR" "$SINGLE_BIN" "$out" ${memory_args[@]+"${memory_args[@]}"}
    mv "$out" "$out.csv"
    echo "[driver] memory and timing results: $out.csv"
    exit 0
fi

# ── stage 2: GMP reference (per machine) ─────────────────────────────────────
GMP_CSV="$RESULTS_DIR/gmp_v2_${MACHINE}_${PROFILE}.csv"
if [[ $FORCE_GMP -eq 1 ]]; then rm -f "$GMP_CSV"; fi
if [[ ! -f "$GMP_CSV" || "$BIN_DIR/bench_run_gmp" -nt "$GMP_CSV" || "$DATA_DIR/manifest.csv" -nt "$GMP_CSV" || $FORCE_DATASET -eq 1 ]]; then
    echo "[gmp] running (no cached result for this machine/profile)"
    reference=$(mktemp "$RESULTS_DIR/gmp.XXXXXX")
    trap 'rm -f "$reference"' EXIT
    "$BIN_DIR/bench_run_gmp" "$DATA_DIR" "$reference" "$MACHINE"
    mv "$reference" "$GMP_CSV"
else
    echo "[gmp] using cached $GMP_CSV"
fi

# ── stage 3: BigMath (always; incrementing run id) ───────────────────────────
RUN_ID=1
shopt -s nullglob
for f in "$RESULTS_DIR"/bigmath_${MACHINE}_${PROFILE}_*.csv; do
    n="${f##*_}"; n="${n%.csv}"
    [[ "$n" =~ ^[0-9]+$ ]] && (( n >= RUN_ID )) && RUN_ID=$((n + 1))
done
shopt -u nullglob
OUT_CSV="$RESULTS_DIR/bigmath_${MACHINE}_${PROFILE}_${RUN_ID}.csv"
echo "[bigmath] run_id=$RUN_ID -> $OUT_CSV"
"$BIN_DIR/bench_run_bigmath" "$DATA_DIR" "$GMP_CSV" "$OUT_CSV" "$RUN_ID" "$MACHINE" "$PROFILE" "$HISTORY_MD"

echo "[driver] done. results: $OUT_CSV   history: $HISTORY_MD"
