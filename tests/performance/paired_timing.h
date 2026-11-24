#ifndef BIGMATH_PAIRED_TIMING
#define BIGMATH_PAIRED_TIMING

#include <algorithm>
#include <array>
#include <chrono>
#include <utility>

namespace bench
{
  template <class FirstCall, class SecondCall>
  std::pair<double, double> MeasureTimes(unsigned calls, FirstCall first, SecondCall second)
  {
    using Clock = std::chrono::steady_clock;
    constexpr std::size_t rounds = 7;
    std::array<double, rounds> firstTimes{}, secondTimes{};
    const auto time = [calls](auto &run) {
      const auto start = Clock::now();
      for (unsigned i = 0; i < calls; ++i)
        run();
      return std::chrono::duration<double, std::micro>(Clock::now() - start).count() / calls;
    };
    for (std::size_t i = 0; i < rounds; ++i)
    {
      if (i % 2 == 0)
      {
        firstTimes[i] = time(first);
        secondTimes[i] = time(second);
      }
      else
      {
        secondTimes[i] = time(second);
        firstTimes[i] = time(first);
      }
    }
    std::sort(firstTimes.begin(), firstTimes.end());
    std::sort(secondTimes.begin(), secondTimes.end());
    return {firstTimes[rounds / 2], secondTimes[rounds / 2]};
  }
}

#endif
