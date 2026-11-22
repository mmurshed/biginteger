# BigMath releases

## Use a prebuilt library

Download the archive for your platform from [GitHub Releases](https://github.com/mmurshed/biginteger/releases). Each archive contains the static `libbigmath.a`, C++ headers, the CMake package, and the license. BigMath requires C++20 and the compiler extensions described in [BUILDING.md](../BUILDING.md).

| Archive suffix | Requirements |
| --- | --- |
| `linux-x86_64-gcc.tar.gz` | Linux x86-64, GCC and GNU libstdc++ compatible with Ubuntu 22.04 |
| `linux-x86_64-clang.tar.gz` | Linux x86-64, Clang and GNU libstdc++ compatible with Ubuntu 24.04 |
| `macos-arm64-clang.tar.gz` | macOS 12 or newer, Apple silicon and Apple Clang |
| `macos-x86_64-clang.tar.gz` | macOS 12 or newer, Intel and Apple Clang |
| `windows-x86_64-clang.zip` | Windows x86-64 and MSYS2 UCRT64 Clang |

Use the same C++ standard library and a compatible compiler and runtime when
building your application. On Windows, run the build from the MSYS2 UCRT64
shell with Clang. For MSVC, a different standard library, or a different
CPU architecture, check the source-build requirements in [BUILDING.md](../BUILDING.md).

Compare the archive's SHA-256 with its entry in `SHA256SUMS`. Use
`shasum -a 256 <archive>` on macOS, `sha256sum <archive>` on Linux, or
`Get-FileHash <archive> -Algorithm SHA256` in PowerShell. Extract the archive
and use the extracted directory containing `include/` and `lib/` as the
CMake prefix.

```cmake
find_package(bigmath CONFIG REQUIRED)
target_link_libraries(my_app PRIVATE bigmath::bigmath)
```

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH="/absolute/path/to/extracted/bigmath"
cmake --build build
```

## Build a snapshot

Open [Actions → Release libraries](https://github.com/mmurshed/biginteger/actions/workflows/release.yml),
choose **Run workflow**, and select a branch. Download the `release-*` artifacts
from the completed run. The same workflow verifies packaging changes in pull
requests. Both paths build, run the four correctness suites, and compile and run
a separate CMake consumer against the extracted archive.

```sh
gh workflow run release.yml --ref main
gh run list --workflow release.yml --limit 5
gh run download RUN_ID --pattern 'release-*'
```

## Publish a version

1. Update the project `VERSION` in [CMakeLists.txt](../CMakeLists.txt) and the
   changelog, then merge those changes into `main`.
2. Pull `main` and push a new annotated tag matching that version. For example,
   after setting the project version to `13.0.1`:

   ```sh
   git pull --ff-only origin main
   git tag -a v13.0.1 -m "BigMath 13.0.1"
   git push origin refs/tags/v13.0.1
   ```

3. Check the tag run in Actions. Once every package passes, GitHub publishes
   the platform archives and `SHA256SUMS` together. For a transient build or
   upload failure, use **Re-run failed jobs** on the tag run. If a failed upload
   left a draft, inspect it in Releases before retrying publication.

The release workflow uses baseline CPU instructions. Source builds continue to
default to native optimization. Dispatch thresholds follow the checked-in
platform profiles at the tagged commit. See [Platform tuning](PLATFORM_TUNING.md).
