# Platform tuning

Dispatch thresholds (when multiplication switches Classic → Karatsuba → Toom →
NTT, when division switches Fast → Burnikel-Ziegler → Newton, etc.) have
crossover points that depend on the CPU and compiler. BigMath ships portable
defaults and lets each platform override them with a tuned profile.

## How resolution works

`common/Constants.h` includes, in order (first definition wins; every threshold
macro is `#ifndef`-guarded):

1. command-line `-D` overrides
2. `build/PlatformConfig.h` → the per-platform tuned profile, if one exists
3. `build/DispatchThresholds.h` → generic portable defaults

`PlatformConfig.h` auto-selects a profile with `__has_include`, keyed on the host
OS / arch / compiler:

```
include/biginteger/build/platform/<os>-<arch>-<compiler>.h
```

e.g. `linux-x86_64-gcc.h`, `macos-arm64-clang.h`, `windows-x86_64-clang.h`. If
the matching file isn't committed, the defaults apply unchanged — no profile is
required to build. A profile only `#define`s the macros it tuned, so it overrides
individual defaults without disabling the fallback for the rest.

Escape hatches: `-DBIGMATH_PLATFORM_OVERRIDE='"path/to/profile.h"'` forces a
specific profile; `-DBIGMATH_PLATFORM_NONE` disables auto-selection.

## Generating a profile

### On real hardware (recommended)

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release        # adds -march=native
cmake --build build --target dispatch_tuner -j
./build/dispatch_tuner --full --emit-header \
    include/biginteger/build/platform/<os>-<arch>-<compiler>.h
```

`--full` runs a wider size sweep (slower, more accurate). Drop it for a quick
pass. Commit the emitted file; `PlatformConfig.h` finds it automatically on the
next build.

### Via CI

The **tune** workflow (`.github/workflows/tune.yml`, manual `workflow_dispatch`)
runs the tuner across the CI matrix (Linux gcc/clang, macOS clang, Windows
clang) and opens a single PR with the regenerated profiles.

> ⚠️ **CI values are not canonical.** GitHub-hosted runners are shared VMs whose
> host CPU varies between runs, and the build uses `-march=native`. A profile
> generated on CI reflects *that runner's* CPU. Use the workflow for convenience
> and as a baseline; for values you rely on, tune on the deployment hardware.

## Adding a new platform

1. Add a matrix entry (and, if a new OS, a job) to `tune.yml` with a `key`
   matching the `<os>-<arch>-<compiler>.h` naming.
2. Add the corresponding auto-selection branch to `PlatformConfig.h`.

The Windows path uses MSYS2 + Clang because the library requires a GCC/Clang
toolchain (`__int128`, `__builtin_*_overflow`); MSVC cannot compile it.
