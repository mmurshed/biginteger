# Platform tuning profiles

Auto-generated dispatch-threshold profiles, one per `<os>-<arch>-<compiler>`:

| File | Host |
|------|------|
| `linux-x86_64-gcc.h`     | Linux / x86-64 / GCC |
| `linux-x86_64-clang.h`   | Linux / x86-64 / Clang |
| `linux-arm64-gcc.h`      | Linux / ARM64 / GCC |
| `linux-arm64-clang.h`    | Linux / ARM64 / Clang |
| `macos-arm64-clang.h`    | macOS / Apple Silicon / Clang |
| `macos-x86_64-clang.h`   | macOS / Intel / Clang |
| `windows-x86_64-clang.h` | Windows / x86-64 / Clang (MSYS2) |

`../PlatformConfig.h` includes the matching file automatically via
`__has_include` — if a profile for the host isn't present, the generic defaults
in `../DispatchThresholds.h` apply. Each profile only `#define`s the macros it
tuned (every one `#ifndef`-guarded), so it overrides the defaults without
disabling the fallback.

## Generating / refreshing a profile

These files are produced by `tests/performance/dispatch_tuner`. Don't hand-edit.

- On real target hardware:
  ```
  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
  cmake --build build --target dispatch_tuner -j
  ./build/dispatch_tuner --full --emit-header \
      include/biginteger/build/platform/<os>-<arch>-<compiler>.h
  ```
- Or trigger the **Tune** workflow (`.github/workflows/tune.yml`,
  `workflow_dispatch`) to run the tuner on the CI matrix and open a PR with the
  regenerated profiles.

See [docs/PLATFORM_TUNING.md](../../../../docs/PLATFORM_TUNING.md).

> **Note on CI-tuned values.** GitHub-hosted runners are shared VMs whose host
> CPU varies between runs, and `-march=native` builds tune for that CPU.
> Profiles generated on CI are a baseline, not canonical — run the tuner on the
> actual deployment hardware for values you depend on.
