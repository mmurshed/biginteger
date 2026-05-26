cmake_minimum_required(VERSION 3.16)
project(bigmath
    VERSION 11.0.0
    DESCRIPTION "Arbitrary-precision integer library with NTT multiplication and Newton division"
    LANGUAGES CXX)

# ─── language requirements ────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)   # rely on standard C++20; library uses __uint128_t
                                # and __builtin_*_overflow which Clang/GCC provide
                                # as extensions but compile under -std=c++20

# ─── default build type ───────────────────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# ─── compiler flags ───────────────────────────────────────────────────────────
# -march=native is critical: enables UMULH/CSEL optimal codegen on ARM64 and
# AVX2/BMI2 on x86. The library specifically targets it. Without it, the
# Karatsuba leaf and NTT butterfly regress noticeably.
if(NOT CMAKE_CXX_FLAGS_RELEASE MATCHES "march=")
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG")
    endif()
endif()

# ─── options ──────────────────────────────────────────────────────────────────
option(BIGMATH_BUILD_TESTS "Build correctness and performance tests" ON)
option(BIGMATH_BUILD_SHARED "Also build a shared-library variant" OFF)

# ─── library target ───────────────────────────────────────────────────────────
file(GLOB_RECURSE BIGMATH_SOURCES CONFIGURE_DEPENDS src/*.cpp)
file(GLOB_RECURSE BIGMATH_HEADERS CONFIGURE_DEPENDS include/biginteger/*.h)

add_library(bigmath STATIC ${BIGMATH_SOURCES} ${BIGMATH_HEADERS})
add_library(bigmath::bigmath ALIAS bigmath)

target_include_directories(bigmath
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

target_compile_features(bigmath PUBLIC cxx_std_20)

# Shared library variant (optional; mirrors object code).
if(BIGMATH_BUILD_SHARED)
    add_library(bigmath_shared SHARED ${BIGMATH_SOURCES} ${BIGMATH_HEADERS})
    add_library(bigmath::bigmath_shared ALIAS bigmath_shared)
    target_include_directories(bigmath_shared
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>)
    target_compile_features(bigmath_shared PUBLIC cxx_std_20)
    set_target_properties(bigmath_shared PROPERTIES
        OUTPUT_NAME bigmath
        SOVERSION ${PROJECT_VERSION_MAJOR}
        VERSION ${PROJECT_VERSION})
endif()

# ─── tests ────────────────────────────────────────────────────────────────────
if(BIGMATH_BUILD_TESTS)
    enable_testing()

    function(bigmath_add_test name source)
        add_executable(${name} ${source})
        target_link_libraries(${name} PRIVATE bigmath::bigmath)
        add_test(NAME ${name} COMMAND ${name})
    endfunction()

    # Correctness harnesses
    bigmath_add_test(mult_correctness tests/mult_correctness.cpp)
    bigmath_add_test(div_correctness  tests/div_correctness.cpp)

    # Aggregated unit-test binary (one TU per area, single REGISTER_TEST registry)
    file(GLOB BIGMATH_UNIT_TEST_SOURCES CONFIGURE_DEPENDS tests/unit/*.cpp)
    add_executable(unit_tests ${BIGMATH_UNIT_TEST_SOURCES})
    target_link_libraries(unit_tests PRIVATE bigmath::bigmath)
    target_include_directories(unit_tests PRIVATE tests/unit)
    add_test(NAME unit_tests COMMAND unit_tests)

    # Performance harnesses (no `add_test` because they print, not assert)
    add_executable(multperf_simple tests/multperf_simple.cpp)
    target_link_libraries(multperf_simple PRIVATE bigmath::bigmath)

    # Interactive calculator REPL (not registered as a ctest target).
    add_executable(calculator calculator/calculator.cpp)
    target_link_libraries(calculator PRIVATE bigmath::bigmath)

    add_executable(divperf_simple tests/divperf_simple.cpp)
    target_link_libraries(divperf_simple PRIVATE bigmath::bigmath)

    # GMP-comparison bench (requires libgmp)
    find_library(GMP_LIBRARY gmp)
    find_path(GMP_INCLUDE_DIR gmp.h)
    if(GMP_LIBRARY AND GMP_INCLUDE_DIR)
        add_executable(bench_vs_gmp tests/performance/bench_vs_gmp.cpp)
        target_link_libraries(bench_vs_gmp PRIVATE bigmath::bigmath ${GMP_LIBRARY})
        target_include_directories(bench_vs_gmp PRIVATE ${GMP_INCLUDE_DIR})
    else()
        message(STATUS "GMP not found; bench_vs_gmp target skipped")
    endif()
endif()

# ─── install ──────────────────────────────────────────────────────────────────
include(GNUInstallDirs)
install(TARGETS bigmath
    EXPORT bigmath-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/biginteger DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
