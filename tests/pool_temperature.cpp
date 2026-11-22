#include "unit_test_framework.h"
#include "biginteger/common/Parallel.h"

#include <atomic>
#include <chrono>
#include <cstdint>
#include <latch>
#include <thread>

namespace
{
  using Clock = std::chrono::steady_clock;
  std::atomic<std::int64_t> elapsedMs{0};
#ifdef BIGMATH_POOL_SPIN_US
  constexpr bool spinEnabled = BIGMATH_POOL_SPIN_US > 0;
#else
  constexpr bool spinEnabled = true;
#endif

  bool WaitForTemperature(bool warm)
  {
    const auto deadline = Clock::now() + std::chrono::seconds(5);
    while (BigMath::ParallelWarm() != warm)
    {
      if (Clock::now() >= deadline)
        return false;
      // Only advance time when asking workers to park. After a prewarm,
      // their spin window stays open until the test has observed it.
      if (!warm)
        elapsedMs.fetch_add(1000);
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    return true;
  }
}

namespace BigMath
{
  std::chrono::steady_clock::time_point PoolSpinNow() noexcept
  {
    return Clock::time_point{} + std::chrono::milliseconds(elapsedMs.load());
  }
}

REGISTER_TEST(ParallelPool, TemperatureTracksParkAndPrewarm)
{
  using namespace BigMath;
  const auto threads = static_cast<Int>(ParallelNumThreads());
  if (threads == 1 || !spinEnabled)
  {
    ParallelPrewarm();
    ASSERT_TRUE(!ParallelWarm());
    return;
  }

  // Every worker must have started before initial parking. A worker that
  // starts late could otherwise be mistaken for one woken by the prewarm.
  std::latch started(threads);
  ParallelDo(threads, [&](Int, Int) { started.arrive_and_wait(); });
  ASSERT_TRUE(ParallelWarm());
  ASSERT_TRUE(WaitForTemperature(false));
  ParallelPrewarm();
  ASSERT_TRUE(WaitForTemperature(true));
  ASSERT_TRUE(WaitForTemperature(false));
}
