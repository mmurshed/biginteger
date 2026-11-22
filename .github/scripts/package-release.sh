#!/usr/bin/env bash
# Run from the source root, with the platform compiler already on PATH.
set -euo pipefail

package="bigmath-${PACKAGE_VERSION:?}-${PACKAGE_KEY:?}"
build="$PWD/build/release"
stage="$build/stage"
dist="$PWD/dist"
unpacked="$build/unpacked"

platform_flags=()
if [[ "$PACKAGE_KEY" == macos-* ]]; then
    platform_flags+=(-DCMAKE_OSX_DEPLOYMENT_TARGET=12.0)
fi

cmake -S . -B "$build" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release -DBIGMATH_NATIVE_OPTIMIZATION=OFF \
    "-DCMAKE_CXX_FLAGS_RELEASE=-O3 ${BIGMATH_CPU:?} -DNDEBUG" "${platform_flags[@]}"
cmake --build "$build" --parallel 3 --target mult_correctness div_correctness unit_tests mfa_roundtrip
ctest --test-dir "$build" --output-on-failure --timeout 600
cmake --install "$build" --prefix "$stage/$package"

mkdir -p "$dist" "$unpacked"
if [[ "$PACKAGE_KEY" == windows-* ]]; then
    archive="$dist/$package.zip"
    (cd "$stage" && cmake -E tar cf "$archive" --format=zip "$package")
else
    archive="$dist/$package.tar.gz"
    (cd "$stage" && cmake -E tar czf "$archive" "$package")
fi
(cd "$unpacked" && cmake -E tar xf "$archive")
test -f "$unpacked/$package/share/bigmath/LICENSE"

smoke="$build/package-smoke"
cmake -S .github/package-smoke -B "$smoke" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release "-DCMAKE_PREFIX_PATH=$unpacked/$package"
cmake --build "$smoke"
ctest --test-dir "$smoke" --output-on-failure
