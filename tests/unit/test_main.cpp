#include "unit_test_framework.h"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <string>

using namespace bigmath_ut;

namespace
{
  bool MatchFilter(const TestCase &tc, const char *filter)
  {
    if (!filter || !*filter)
      return true;
    std::string full = std::string(tc.suite) + "." + tc.name;
    return std::strstr(full.c_str(), filter) != nullptr;
  }
}

int main(int argc, char **argv)
{
  const char *filter = nullptr;
  for (int i = 1; i < argc; ++i)
  {
    if (std::strncmp(argv[i], "--filter=", 9) == 0)
      filter = argv[i] + 9;
    else if (std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0)
    {
      std::printf("usage: %s [--filter=substr]\n", argv[0]);
      return 0;
    }
  }

  auto &reg = Registry();
  int passed = 0, failed = 0, skipped = 0;
  auto t0 = std::chrono::steady_clock::now();

  for (const auto &tc : reg)
  {
    if (!MatchFilter(tc, filter))
    {
      ++skipped;
      continue;
    }
    auto ts = std::chrono::steady_clock::now();
    try
    {
      tc.fn();
      auto te = std::chrono::steady_clock::now();
      double ms = std::chrono::duration<double, std::milli>(te - ts).count();
      std::printf("  [PASS] %s.%s  (%.2f ms)\n", tc.suite, tc.name, ms);
      ++passed;
    }
    catch (const AssertionFailure &af)
    {
      std::printf("  [FAIL] %s.%s\n    %s\n", tc.suite, tc.name, af.message.c_str());
      ++failed;
    }
    catch (const std::exception &e)
    {
      std::printf("  [FAIL] %s.%s\n    uncaught std::exception: %s\n", tc.suite, tc.name, e.what());
      ++failed;
    }
    catch (...)
    {
      std::printf("  [FAIL] %s.%s\n    uncaught non-std exception\n", tc.suite, tc.name);
      ++failed;
    }
  }

  auto t1 = std::chrono::steady_clock::now();
  double total_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

  std::printf("\n");
  std::printf("===========================================================\n");
  std::printf(" Tests: %d passed, %d failed", passed, failed);
  if (skipped) std::printf(", %d skipped", skipped);
  std::printf("   (%.1f ms total)\n", total_ms);
  std::printf("===========================================================\n");

  return failed == 0 ? 0 : 1;
}
