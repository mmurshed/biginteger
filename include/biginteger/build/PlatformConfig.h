#ifndef BIGMATH_PLATFORM_CONFIG
#define BIGMATH_PLATFORM_CONFIG

// Platform-specific tuned dispatch thresholds.
//
// Pulled in by common/Constants.h *before* build/DispatchThresholds.h. A tuned
// profile only #defines the threshold macros it measured (each #ifndef-guarded),
// so it wins over the generic defaults without disabling the fallback for any
// macro it did not set.
//
// Profiles live under include/biginteger/build/platform/ and are named
//
//   <os>-<arch>-<compiler>.h      e.g. linux-x86_64-gcc.h, macos-arm64-clang.h
//
// Selection is automatic: the branches below match the host OS/arch/compiler and
// include the matching profile *if it exists* (via __has_include). With no
// profile committed for the host, the defaults in DispatchThresholds.h apply
// unchanged — no profile is required to build.
//
// Generate a profile with the tuner (see docs/PLATFORM_TUNING.md and
// .github/workflows/tune.yml):
//
//   dispatch_tuner --full --emit-header \
//     include/biginteger/build/platform/<os>-<arch>-<compiler>.h
//
// Escape hatches:
//   -DBIGMATH_PLATFORM_OVERRIDE='"path/to/profile.h"'  force a specific profile
//   -DBIGMATH_PLATFORM_NONE                            disable auto-selection

#if defined(BIGMATH_PLATFORM_OVERRIDE)

#include BIGMATH_PLATFORM_OVERRIDE

#elif !defined(BIGMATH_PLATFORM_NONE)

// ─── architecture detection ───────────────────────────────────────────────────
#if defined(__x86_64__) || defined(_M_X64)
#define BIGMATH_PLATFORM_X86_64 1
#elif defined(__aarch64__) || defined(_M_ARM64)
#define BIGMATH_PLATFORM_ARM64 1
#endif

// ─── auto-selection (clang checked before gcc: clang also defines __GNUC__) ────
#if defined(__APPLE__) && defined(BIGMATH_PLATFORM_ARM64) && defined(__clang__)
#if __has_include("biginteger/build/platform/macos-arm64-clang.h")
#include "biginteger/build/platform/macos-arm64-clang.h"
#endif

#elif defined(__APPLE__) && defined(BIGMATH_PLATFORM_X86_64) && defined(__clang__)
#if __has_include("biginteger/build/platform/macos-x86_64-clang.h")
#include "biginteger/build/platform/macos-x86_64-clang.h"
#endif

#elif defined(_WIN32) && defined(BIGMATH_PLATFORM_X86_64) && defined(__clang__)
#if __has_include("biginteger/build/platform/windows-x86_64-clang.h")
#include "biginteger/build/platform/windows-x86_64-clang.h"
#endif

#elif defined(__linux__) && defined(BIGMATH_PLATFORM_X86_64) && defined(__clang__)
#if __has_include("biginteger/build/platform/linux-x86_64-clang.h")
#include "biginteger/build/platform/linux-x86_64-clang.h"
#endif

#elif defined(__linux__) && defined(BIGMATH_PLATFORM_X86_64) && defined(__GNUC__) && !defined(__clang__)
#if __has_include("biginteger/build/platform/linux-x86_64-gcc.h")
#include "biginteger/build/platform/linux-x86_64-gcc.h"
#endif

#elif defined(__linux__) && defined(BIGMATH_PLATFORM_ARM64) && defined(__clang__)
#if __has_include("biginteger/build/platform/linux-arm64-clang.h")
#include "biginteger/build/platform/linux-arm64-clang.h"
#endif

#elif defined(__linux__) && defined(BIGMATH_PLATFORM_ARM64) && defined(__GNUC__) && !defined(__clang__)
#if __has_include("biginteger/build/platform/linux-arm64-gcc.h")
#include "biginteger/build/platform/linux-arm64-gcc.h"
#endif

#endif

#endif // selection

#endif
