/**
 * Experimental 64-bit limb addition and subtraction with the BigMath pool.
 * Forces the serial or shared parallel magnitude kernel for comparison.
 */

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdio>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <string_view>
#include <thread>
#include <vector>

#include "biginteger/algorithms/Addition.h"
#include "biginteger/algorithms/Subtraction.h"
#include "biginteger/algorithms/ParallelAddSub.h"
#include "biginteger/common/Parallel.h"

using namespace BigMath;

namespace
{
  using Clock = std::chrono::steady_clock;
  using namespace BigMath::AddSubDetail;

  template <bool Subtract>
  void Check(Limbs a, Limbs b, SizeT workers)
  {
    const auto size = static_cast<SizeT>(a.size());
    std::vector<DataT> expected(size), actual(size);
    const DataT flag = SerialBlock<Subtract>(expected.data(), a, b, 0, size);
    if (ParallelMagnitude<Subtract>(actual.data(), a, b, workers) != flag || actual != expected)
      throw std::runtime_error("parallel magnitude differs from serial");
    std::copy(a.begin(), a.end(), actual.begin());
    if (ParallelMagnitude<Subtract>(actual.data(), actual, b, workers) != flag || actual != expected)
      throw std::runtime_error("parallel magnitude differs with left alias");
    std::copy(b.begin(), b.end(), actual.begin());
    if (ParallelMagnitude<Subtract>(actual.data(), a, {actual.data(), b.size()}, workers) != flag ||
        actual != expected)
      throw std::runtime_error("parallel magnitude differs with right alias");
  }

  void Verify(SizeT workers)
  {
    std::mt19937_64 random(7);
    for (unsigned trial = 0; trial < 240; ++trial)
    {
      std::vector<DataT> a(random() % 4097), b(random() % (a.size() + 1));
      for (auto &limb : a)
        limb = random();
      for (auto &limb : b)
        limb = random();
      for (SizeT tasks = 1; tasks <= workers; ++tasks)
      {
        Check<false>(a, b, tasks);
        Check<true>(a, b, tasks);
      }
    }

    // Cross every block boundary and force the final flag. Then terminate
    // the chain at different positions inside and at the edge of a block.
    for (SizeT size : {43u, 65'537u, 524'288u})
    {
      std::vector<DataT> ones(size, std::numeric_limits<DataT>::max()), zeros(size), one{1};
      for (SizeT tasks = 1; tasks <= workers; ++tasks)
      {
        Check<false>(ones, one, tasks);
        Check<true>(zeros, one, tasks);
        for (SizeT pos : {size / tasks, size / 2, size - 1})
        {
          if (pos == size)
            continue;
          --ones[pos];
          ++zeros[pos];
          Check<false>(ones, one, tasks);
          Check<true>(zeros, one, tasks);
          ++ones[pos];
          --zeros[pos];
        }
      }
    }
  }

  void PrepareWorkers(bool warm)
  {
    const auto deadline = Clock::now() + std::chrono::milliseconds(100);
    if (warm)
    {
      do
      {
        ParallelPrewarm();
        if (ParallelWarm())
          return;
        std::this_thread::yield();
      } while (Clock::now() < deadline);
    }
    else
    {
      // Leave time for workers to park, then verify the pool reports cold.
      do
      {
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        if (!ParallelWarm())
          return;
      } while (Clock::now() < deadline);
    }
    throw std::runtime_error("could not establish benchmark worker state");
  }

  template <class Serial, class Parallel>
  void Measure(const char *op, SizeT limbs, SizeT workers, std::string_view shape,
               const char *layer, bool warm, Serial serial, Parallel parallel)
  {
    constexpr std::size_t Repetitions = 31;
    std::array<double, Repetitions> serialTimes, parallelTimes, ratios;
    const auto time = [warm](auto &call) {
      PrepareWorkers(warm);
      const auto start = Clock::now();
      call();
      return std::chrono::duration<double, std::micro>(Clock::now() - start).count();
    };
    for (std::size_t i = 0; i < Repetitions; ++i)
    {
      if (i % 2 == 0)
      {
        serialTimes[i] = time(serial);
        parallelTimes[i] = time(parallel);
      }
      else
      {
        parallelTimes[i] = time(parallel);
        serialTimes[i] = time(serial);
      }
      ratios[i] = serialTimes[i] / parallelTimes[i];
    }
    const auto median = [](auto &values) {
      std::sort(values.begin(), values.end());
      return values[values.size() / 2];
    };
    std::printf("%s,%.*s,%s,%s,%u,%u,%.3f,%.3f,%.4f\n", op, static_cast<int>(shape.size()),
                shape.data(), layer, warm ? "warm" : "parked", limbs, workers,
                median(serialTimes), median(parallelTimes), median(ratios));
  }

  template <bool Subtract>
  void Benchmark(SizeT size, SizeT workers, std::string_view shape)
  {
    const char *op = Subtract ? "sub" : "add";
    std::mt19937_64 random(size);
    std::vector<DataT> a(size), b(shape == "skewed" ? size / 4 : size);
    if (shape == "chain")
    {
      if constexpr (Subtract)
        a.back() = 1;
      else
        std::fill(a.begin(), a.end(), std::numeric_limits<DataT>::max());
      b.front() = 1;
    }
    else
    {
      for (auto &limb : a)
        limb = random();
      for (auto &limb : b)
        limb = random();
      // Force nonzero top limbs and a >= b for magnitude subtraction.
      a.back() |= DataT{1} << 63;
      b.back() = (b.back() >> 1) | 1;
    }

    const auto serialResult = [&] {
      if constexpr (Subtract)
        return BigMath::Subtract(a, b, Base2_64);
      else
        return BigMath::Add(a, b, Base2_64);
    };
    const auto expected = serialResult();
    if (ParallelResult<Subtract>(a, b, workers) != expected)
      throw std::runtime_error("allocating prototype differs from the library");

    std::vector<DataT> serialOut(size), parallelOut(size), serialVector, parallelVector;
    DataT serialFlag = 0, parallelFlag = 0;
    for (bool warm : {true, false})
    {
      Measure(op, size, workers, shape, "kernel", warm,
              [&] { serialFlag = SerialBlock<Subtract>(serialOut.data(), a, b, 0, size); },
              [&] { parallelFlag = ParallelMagnitude<Subtract>(parallelOut.data(), a, b, workers); });
      if (serialOut != parallelOut || serialFlag != parallelFlag)
        throw std::runtime_error("timed raw results differ");
      Measure(op, size, workers, shape, "magnitude", warm,
              [&] { serialVector = serialResult(); },
              [&] { parallelVector = ParallelResult<Subtract>(a, b, workers); });
      if (serialVector != expected || parallelVector != expected)
        throw std::runtime_error("timed allocating results differ");
    }
  }
}

int main(int argc, char **argv)
{
  if (argc > 2 || (argc == 2 && std::string_view(argv[1]) != "--check"))
  {
    std::fputs("Usage: parallel_add_sub_bench [--check]\n", stderr);
    return 2;
  }
  const SizeT workers = ParallelNumThreads();
  Verify(workers);
  if (argc == 2)
  {
    std::printf("Matched the serial kernels with 1 through %u blocks, including exact aliases and full carry/borrow chains.\n",
                workers);
    return 0;
  }
  if (workers < 2 || LimbBits != 64)
  {
    std::fputs("The timing sweep needs 64-bit limbs and a pool with at least two execution threads.\n",
               stderr);
    return 2;
  }

  std::puts("operation,shape,layer,pool,limbs,blocks,serial_us,parallel_us,speedup");
  SizeT last = 0;
  for (SizeT blocks : {2u, std::min(4u, workers), workers})
  {
    if (blocks == last)
      continue;
    last = blocks;
    for (SizeT size : {2'048u, 4'096u, 6'144u, 8'192u, 12'288u, 16'384u,
                       24'576u, 32'768u, 49'152u, 65'536u, 98'304u, 131'072u,
                       262'144u, 524'288u, 1'048'576u, 2'097'152u, 4'194'304u})
    {
      Benchmark<false>(size, blocks, "balanced");
      Benchmark<true>(size, blocks, "balanced");
    }
    for (SizeT size : {65'536u, 524'288u, 2'097'152u})
      for (std::string_view shape : {"skewed", "chain"})
      {
        Benchmark<false>(size, blocks, shape);
        Benchmark<true>(size, blocks, shape);
      }
  }
}
