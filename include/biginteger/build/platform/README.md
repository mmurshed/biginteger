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

## Updating a profile

Run the dispatch tuner or the **Tune platform dispatch** workflow to produce
candidate headers. Review the benchmark logs, retain the supported overrides,
and validate the proposed profile before merging it.

See [Platform tuning](../../../../docs/PLATFORM_TUNING.md) for the complete workflow.
