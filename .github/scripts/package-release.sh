#!/usr/bin/env bash
# Run from the source root, with the platform compiler already on PATH.
set -euo pipefail

package="bigmath-${PACKAGE_VERSION:?}-${PLATFORM_KEY:?}"
build="$PWD/build/release"
stage="$build/stage"
dist="$PWD/dist"
unpacked="$build/unpacked"

source .github/scripts/configure-platform.sh
profile=$(platform_profile "$PWD")
configure_platform "$build" "$profile"
cmake --build "$build" --parallel 3 --target mult_correctness div_correctness unit_tests pool_temperature mfa_roundtrip
ctest --test-dir "$build" --output-on-failure --timeout 600
cmake --install "$build" --prefix "$stage/$package"

mkdir -p "$dist" "$unpacked"
if [[ "$PLATFORM_KEY" == windows-* ]]; then
    archive="$dist/$package.zip"
    (cd "$stage" && cmake -E tar cf "$archive" --format=zip "$package")
else
    archive="$dist/$package.tar.gz"
    (cd "$stage" && cmake -E tar czf "$archive" "$package")
fi
(cd "$unpacked" && cmake -E tar xf "$archive")
test -f "$unpacked/$package/share/bigmath/LICENSE"
if [[ "$profile" == none ]]; then
    profile=cmake/DefaultPlatform.h
fi
cmake -E compare_files "$profile" "$unpacked/$package/include/biginteger/build/PackagedPlatform.h"

smoke="$build/package-smoke"
cmake -S .github/package-smoke -B "$smoke" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release "-DCMAKE_PREFIX_PATH=$unpacked/$package"
cmake --build "$smoke"
ctest --test-dir "$smoke" --output-on-failure
