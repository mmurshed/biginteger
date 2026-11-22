#!/usr/bin/env bash
# Common CMake settings for release builds and the dispatch tuner.
configure_platform() {
    local build="$1" profile="$2"
    local deployment_target=
    if [[ "$PLATFORM_KEY" == macos-* ]]; then
        deployment_target=-DCMAKE_OSX_DEPLOYMENT_TARGET=12.0
    fi
    cmake -S . -B "$build" -G Ninja \
        -DCMAKE_BUILD_TYPE=Release -DBIGMATH_NATIVE_OPTIMIZATION=OFF \
        "-DBIGMATH_PLATFORM_PROFILE=$profile" \
        "-DCMAKE_CXX_FLAGS_RELEASE=-O3 ${BIGMATH_CPU:?} -DNDEBUG" ${deployment_target:+"$deployment_target"}
}

platform_profile() {
    local header="$1/include/biginteger/build/platform/$PLATFORM_KEY.h"
    if [[ -f "$header" ]]; then
        printf '%s\n' "$header"
    else
        printf '%s\n' none
    fi
}
