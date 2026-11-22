#!/usr/bin/env bash
# Run from the source root, with the platform compiler already on PATH.
set -euo pipefail
source .github/scripts/configure-platform.sh

out="$PWD/build/tuning-results"
candidate="$PWD/build/tune-candidate"
baseline="$PWD/build/tune-baseline"
baseline_source="$PWD/build/tune-baseline-source"
mkdir -p "$out"
exec > >(tee "$out/run.log") 2>&1

tuner=dispatch_tuner
[[ "$PLATFORM_KEY" != windows-* ]] || tuner+=.exe
profile=$(platform_profile "$PWD")
candidate_commit="${GITHUB_SHA:-$(git rev-parse HEAD)}"
{
    echo "platform=$PLATFORM_KEY"
    echo "candidate_commit=$candidate_commit"
    echo "stage=$TUNING_STAGE"
    echo "mode=$TUNING_MODE"
    echo "release_cpu_flags=$BIGMATH_CPU"
    uname -a
    "$CXX" --version
    case "$PLATFORM_KEY" in
        macos-*) sysctl -n machdep.cpu.brand_string ;;
        linux-*) lscpu ;;
        windows-*) echo "${PROCESSOR_IDENTIFIER:-unknown CPU}" ;;
    esac
} | tee "$out/environment.txt"

record_profile() {
    local header="$1" name="$2"
    if [[ "$header" == none ]]; then
        printf 'Generic defaults\n' > "$out/$name-profile.txt"
    else
        cp "$header" "$out/$name-profile.txt"
    fi
}

record_profile "$profile" candidate
configure_platform "$candidate" "$profile"
if [[ "$TUNING_STAGE" == survey ]]; then
    full_flag=
    if [[ "$TUNING_MODE" == full ]]; then
        full_flag=--full
    fi
    cmake --build "$candidate" --parallel 3 --target dispatch_tuner
    for sample in 1 2 3; do
        "$candidate/$tuner" ${full_flag:+"$full_flag"} \
            --emit-header "$out/survey-$sample.h" | tee "$out/survey-$sample.log"
    done
    exit
fi

base_profile=$(platform_profile "$baseline_source")
baseline_commit="${TUNING_BASELINE_SHA:-$(git -C "$baseline_source" rev-parse HEAD)}"
printf '%s\n' "$baseline_commit" > "$out/baseline-commit.txt"
record_profile "$base_profile" baseline
configure_platform "$baseline" "$base_profile"
cmake --build "$baseline" --parallel 3 --target dispatch_tuner
cmake --build "$candidate" --parallel 3 \
    --target dispatch_tuner mult_correctness div_correctness unit_tests pool_temperature mfa_roundtrip
ctest --test-dir "$candidate" --output-on-failure --timeout 600 | tee "$out/ctest.log"

for sample in 1 2 3; do
    if (( sample % 2 )); then
        variants=(baseline candidate)
    else
        variants=(candidate baseline)
    fi
    for variant in "${variants[@]}"; do
        "$PWD/build/tune-$variant/$tuner" --holdout | tee "$out/$variant-$sample.csv"
    done
done
