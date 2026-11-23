# Platform tuning

BigMath chooses multiplication, squaring, and division algorithms using
size and shape thresholds. The right switch points depend on the CPU,
compiler, and whether worker threads are active. A platform profile overrides
individual thresholds while using the library defaults for the rest.

## Survey the crossovers

Open [Actions → Tune platform dispatch](https://github.com/mmurshed/biginteger/actions/workflows/tune.yml).
Run the workflow with `stage=survey` and `mode=full`. A quick survey uses fewer
sizes. The workflow builds for the same five OS, architecture, and compiler
combinations as the [release workflow](RELEASING.md), with the same CPU instruction
baseline.

```sh
gh workflow run tune.yml --ref main -f stage=survey -f mode=full
gh run list --workflow tune.yml --limit 5
gh run download RUN_ID --pattern 'survey-*' --dir survey
```

Each `survey-<platform>` artifact contains three sweeps (`survey-*.log`), the
suggested headers from each sweep (`survey-*.h`), the profile used for the build,
and compiler, runner, and commit metadata. Times in the survey are in
milliseconds. The logs show competing kernels around multiplication,
squaring, and division crossovers. The tuner compares burst and parked
worker states for multiplication and squaring. Across the current 64-bit
division sweep it checks both quotient and remainder and requires a 5%
challenger lead at two adjacent sizes before proposing a new floor.

Addition and subtraction use separate warm and parked cutoffs for
equal-length inputs. The survey alternates serial and parallel trials for
both random operands and full carry or borrow chains. Each survey
requires a lead of at least 5% for both shapes at two adjacent sizes
and every larger sampled size. The workflow selects the largest cutoff
supported by all three surveys. Zero keeps that operation and worker
state serial, including when the runner cannot reliably establish that
state for timing. The `selected-add-sub-profiles` artifact contains the
selected macro definitions for each platform. Copy those definitions
into the matching profile header while preserving the other overrides.

Use `mode=add-sub` with either stage to run the focused survey or holdout.
Holdout timing uses the public integer operators and different lengths
and operands from the survey. Parked timings include waking workers in
an already initialized pool. The first operation to start the pool also
pays its one-time startup cost.

## Choose and validate a profile

1. Compare all three sweeps. Keep changes with a repeatable crossover, and
   preserve the current values when the measurements disagree. Avoid moving
   unrelated thresholds together so you can tell which change caused a result.
2. On a branch, copy the relevant proposal to
   `include/biginteger/build/platform/<platform>.h`. Keep the limb-width and
   `#ifndef` guards, and remove definitions for thresholds you want to leave
   at their defaults. Commit the selected files and push the branch.
3. Run `stage=validate` on the branch. This builds both the proposed profile
   and the profile from the default branch at dispatch time, using the proposed
   branch's source code and the same runner for both. The candidate must pass
   all five correctness suites. The benchmark also checks results against
   schoolbook multiplication and FastDivision on the held-out inputs.

   ```sh
   gh workflow run tune.yml --ref your-profile-branch -f stage=validate
   gh run list --workflow tune.yml --limit 5
   gh run download RUN_ID --pattern 'validate-*' --dir validation
   gh run download RUN_ID --name validation-comparison --dir validation-summary
   ```

4. Read the workflow summary and the per-shape CSVs. Speedup is baseline time
   divided by candidate time, so values above 1 favor the candidate. The
   comparison uses the median of three paired runs and includes the slowest
   case in each operation and worker state. Review any regressions, including
   when the average improves. If the evidence is mixed, adjust the branch and
   repeat validation. Link the runs and summarize the changes in the PR before
   merging.

Hosted runners can vary by CPU and load even with the same runner label. Treat
the results as a baseline for the published packages. For deployment-specific
thresholds, run the tuner and validate on the deployment hardware.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target dispatch_tuner --parallel
./build/dispatch_tuner --full --emit-header /tmp/bigmath-candidate.h
./build/dispatch_tuner --holdout
```

Use `dispatch_tuner.exe` in the Windows MSYS2 UCRT64 shell. To build locally
with a specific profile, configure a clean build with
`-DBIGMATH_PLATFORM_PROFILE=/absolute/path/to/profile.h`, then rebuild the
tuner, library and tests. Use `-DBIGMATH_PLATFORM_PROFILE=none` for the generic
baseline. The native optimization default in local Release builds can be
disabled with `-DBIGMATH_NATIVE_OPTIMIZATION=OFF`.

## How builds pick up the profile

A release package embeds the profile matching its OS, architecture and
compiler at the tagged commit. Its exported `bigmath::bigmath` target uses
the same profile for library headers seen by applications. If no matching
header exists, the package embeds the generic defaults. Use the exported
CMake target when consuming release binaries.

A source build without an explicit CMake profile auto-selects a checked-in
header for the compiler's target. Command-line threshold macros take
precedence, and absent settings fall through to
[`DispatchThresholds.h`](../include/biginteger/build/DispatchThresholds.h).
For source builds, `BIGMATH_PLATFORM_OVERRIDE` can force a profile header and
`BIGMATH_PLATFORM_NONE` disables profile selection.
