#!/usr/bin/env bash
# Driver for the split BigMath-vs-GMP benchmark.
#
#   ./run_benchmark.sh [profile] [--force-dataset] [--force-gmp]
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
for arg in "$@"; do
    case "$arg" in
        quick|default|full) PROFILE="$arg" ;;
        --force-dataset)    FORCE_DATASET=1 ;;
        --force-gmp)        FORCE_GMP=1 ;;
        *) echo "unknown arg: $arg" >&2; exit 2 ;;
    esac
done

BENCH_DIR="${BENCH_DIR:-$REPO/.benchcache}"
BIN_DIR="$BENCH_DIR/bin"
DATASET_ROOT="$BENCH_DIR/dataset"
DATA_DIR="$DATASET_ROOT/$PROFILE"
ZIP="$BENCH_DIR/dataset_$PROFILE.zip"
RESULTS_DIR="$BENCH_DIR/results"
HISTORY_MD="$SCRIPT_DIR/HISTORY.md"
mkdir -p "$BIN_DIR" "$DATASET_ROOT" "$RESULTS_DIR"

# ── machine identity ─────────────────────────────────────────────────────────
if [[ "$(uname)" == "Darwin" ]]; then
    CPU="$(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo unknown)"
else
    CPU="$(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2 || echo unknown)"
fi
MACHINE="$(printf '%s_%s' "$CPU" "$(uname -m)" | tr -cs 'A-Za-z0-9' '_' | sed 's/_*$//')"
echo "[driver] machine=$MACHINE profile=$PROFILE bench_dir=$BENCH_DIR"

# ── compiler setup ───────────────────────────────────────────────────────────
CXX="${CXX:-c++}"
CXXFLAGS="-std=c++20 -O3 -march=native -DNDEBUG -I$REPO/include -I$SCRIPT_DIR"
GMP_PREFIX=""
for p in /opt/homebrew /usr/local; do
    [[ -f "$p/include/gmp.h" ]] && GMP_PREFIX="$p" && break
done
GMP_FLAGS="-lgmp"
[[ -n "$GMP_PREFIX" ]] && GMP_FLAGS="-I$GMP_PREFIX/include -L$GMP_PREFIX/lib -lgmp"

# The library is not header-only: Parse/ToString/operator*//, BigDecimal, etc.
# are defined under src/ and bigdecimal/. run_bigmath must compile + link them
# (gen_dataset and run_gmp reference no BigMath symbols, so they don't need this).
LIB_SRCS="$(find "$REPO/src" -name '*.cpp'; echo "$REPO/bigdecimal/BigDecimal.cpp")"
LIB_SRCS="$(echo $LIB_SRCS)"  # collapse newlines to spaces

build() { # <out> <src> [extra-flags-string]
    local out="$BIN_DIR/$1"; local src="$SCRIPT_DIR/$2"; local extra="${3:-}"
    if [[ ! -x "$out" || "$src" -nt "$out" || "$SCRIPT_DIR/bench_common.h" -nt "$out" ]]; then
        echo "[build] $1"
        # $extra intentionally unquoted so multiple flags word-split.
        $CXX $CXXFLAGS "$src" -o "$out" $extra
    fi
}

# ── stage 1: dataset ─────────────────────────────────────────────────────────
if [[ $FORCE_DATASET -eq 1 ]]; then rm -rf "$DATA_DIR" "$ZIP"; fi
if [[ ! -f "$DATA_DIR/manifest.csv" ]]; then
    if [[ -f "$ZIP" ]]; then
        echo "[dataset] extracting cached $ZIP"
        unzip -q -o "$ZIP" -d "$DATASET_ROOT"
    else
        echo "[dataset] generating ($PROFILE)"
        build gen_dataset gen_dataset.cpp
        mkdir -p "$DATA_DIR"
        "$BIN_DIR/gen_dataset" "$DATA_DIR" "$PROFILE"
        echo "[dataset] zipping -> $ZIP"
        ( cd "$DATASET_ROOT" && zip -q -r -1 "$ZIP" "$PROFILE" )
    fi
else
    echo "[dataset] using existing $DATA_DIR"
fi

# ── stage 2: GMP reference (per machine) ─────────────────────────────────────
GMP_CSV="$RESULTS_DIR/gmp_${MACHINE}_${PROFILE}.csv"
if [[ $FORCE_GMP -eq 1 ]]; then rm -f "$GMP_CSV"; fi
if [[ ! -f "$GMP_CSV" ]]; then
    echo "[gmp] running (no cached result for this machine/profile)"
    build run_gmp run_gmp.cpp "$GMP_FLAGS"
    "$BIN_DIR/run_gmp" "$DATA_DIR" "$GMP_CSV" "$MACHINE"
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
build run_bigmath run_bigmath.cpp "$LIB_SRCS"
"$BIN_DIR/run_bigmath" "$DATA_DIR" "$GMP_CSV" "$OUT_CSV" "$RUN_ID" "$MACHINE" "$PROFILE" "$HISTORY_MD"

echo "[driver] done. results: $OUT_CSV   history: $HISTORY_MD"
