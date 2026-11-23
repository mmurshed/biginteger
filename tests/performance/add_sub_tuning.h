#ifndef BIGMATH_ADD_SUB_TUNING
#define BIGMATH_ADD_SUB_TUNING

#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <thread>

#include "biginteger/algorithms/ParallelAddSub.h"
#include "biginteger/ops/Operations.h"

namespace BigMath::AddSubTuning
{
  struct Cutoffs
  {
    SizeT addWarm = 0, addParked = 0, subWarm = 0, subParked = 0;
  };

  using Clock = std::chrono::steady_clock;

  inline bool WarmPool()
  {
    if (ParallelWarm())
      return true;
    const auto workers = static_cast<Int>(ParallelNumThreads());
    const auto deadline = Clock::now() + std::chrono::milliseconds(250);
    do
    {
      // On oversubscribed runners a wake-up alone may leave too little time
      // for enough workers to spin simultaneously. Give them an untimed job.
      ParallelDo(workers, [](Int, Int) {
        const auto until = Clock::now() + std::chrono::microseconds(150);
        while (Clock::now() < until) {}
      });
      if (ParallelWarm())
        return true;
    } while (Clock::now() < deadline);
    return false;
  }

  inline void PreparePool(bool warm)
  {
    if (!warm)
    {
      do { std::this_thread::sleep_for(std::chrono::milliseconds(1)); } while (ParallelWarm());
      return;
    }
    if (!WarmPool())
      throw std::runtime_error("could not warm pool for add/sub tuning");
  }

  // Only recommend a cutoff if both shapes continue to win at every larger
  // sampled size, with at least two adjacent wins. Zero leaves the path off.
  struct Crossover
  {
    SizeT first = 0, wins = 0;
    void Observe(SizeT limbs, double speedup)
    {
      if (speedup < 1.05)
      {
        first = 0;
        wins = 0;
        return;
      }
      if (wins++ == 0)
        first = limbs;
    }
    SizeT Cutoff() const { return wins >= 2 ? first : 0; }
  };

  template <bool Subtract>
  auto Operands(SizeT size, bool chain, std::mt19937_64 &random)
  {
    std::vector<DataT> a(size), b(size);
    if (chain)
    {
      if constexpr (!Subtract)
        std::fill(a.begin(), a.end(), LimbMask);
      a.back() = 2;
      b.back() = b.front() = 1;
    }
    else
    {
      for (auto &word : a) word = random() & LimbMask;
      for (auto &word : b) word = random() & LimbMask;
      a.back() |= DataT{1} << (LimbBits - 1);
      b.back() = (b.back() >> 1) | 1;
    }
    return std::pair(std::move(a), std::move(b));
  }

  template <bool Subtract>
  auto Reference(std::vector<DataT> const &a, std::vector<DataT> const &b)
  {
    if constexpr (Subtract)
      return BigMath::Subtract(a, b, CurrentBase);
    else
      return BigMath::Add(a, b, CurrentBase);
  }

  template <class Call>
  double TimeUs(Call call)
  {
    const auto start = Clock::now();
    call();
    return std::chrono::duration<double, std::micro>(Clock::now() - start).count();
  }

  template <std::size_t Count>
  double Median(std::array<double, Count> samples)
  {
    std::sort(samples.begin(), samples.end());
    return samples[Count / 2];
  }

  template <bool Subtract>
  double Measure(SizeT size, SizeT blocks, bool warm, bool chain, std::mt19937_64 &random)
  {
    const auto [a, b] = Operands<Subtract>(size, chain, random);
    const auto expected = Reference<Subtract>(a, b);
    std::vector<DataT> serial, parallel;
    const auto timeSerial = [&] {
      PreparePool(warm);
      return TimeUs([&] { serial = Reference<Subtract>(a, b); });
    };
    const auto timeParallel = [&] {
      PreparePool(warm);
      return TimeUs([&] { parallel = AddSubDetail::ParallelResult<Subtract>(a, b, blocks); });
    };
    std::array<double, 11> serialTimes, parallelTimes, ratios;
    for (std::size_t i = 0; i < ratios.size(); ++i)
    {
      if (i % 2 == 0)
      {
        serialTimes[i] = timeSerial();
        parallelTimes[i] = timeParallel();
      }
      else
      {
        parallelTimes[i] = timeParallel();
        serialTimes[i] = timeSerial();
      }
      ratios[i] = serialTimes[i] / parallelTimes[i];
      if (serial != expected || parallel != expected)
        throw std::runtime_error("add/sub survey result mismatch");
    }

    const double speedup = Median(ratios);
    std::cout << (Subtract ? "sub" : "add") << ',' << (warm ? "warm" : "parked") << ','
              << (chain ? "chain" : "balanced") << ',' << size << ',' << blocks << ','
              << std::fixed << std::setprecision(6) << Median(serialTimes) << ','
              << Median(parallelTimes) << ',' << speedup << '\n';
    return speedup;
  }

  template <bool Subtract>
  SizeT Survey(std::vector<SizeT> const &sizes, SizeT workers, bool warm)
  {
    Crossover crossover;
    std::mt19937_64 random(0xADD5AB);
    for (SizeT size : sizes)
    {
      const auto blocks = AddSubDetail::BlockCount(size, workers);
      if (blocks <= 1)
        continue;
      const double balanced = Measure<Subtract>(size, blocks, warm, false, random);
      const double chain = Measure<Subtract>(size, blocks, warm, true, random);
      crossover.Observe(size, std::min(balanced, chain));
    }
    return crossover.Cutoff();
  }

  inline Cutoffs Tune(bool full)
  {
    const SizeT workers = ParallelNumThreads();
    if (LimbBits != 64 || workers <= 1)
      return {};

    std::vector<SizeT> sizes{4096, 6144, 8192, 12288, 16384, 24576, 32768,
                             49152, 65536, 98304, 131072, 262144, 524288, 1048576};
    if (full)
      sizes.insert(sizes.end(), {2097152, 4194304});

    std::cout << "\n== Parallel Addition and Subtraction ==\n"
              << "operation,state,shape,limbs,blocks,serial_us,parallel_us,paired_speedup\n";
    Cutoffs cutoffs;
    if (WarmPool())
    {
      cutoffs.addWarm = Survey<false>(sizes, workers, true);
      cutoffs.subWarm = Survey<true>(sizes, workers, true);
    }
    cutoffs.addParked = Survey<false>(sizes, workers, false);
    cutoffs.subParked = Survey<true>(sizes, workers, false);
    return cutoffs;
  }

  template <bool Subtract>
  void HoldoutCase(SizeT size, bool chain, bool warm, std::mt19937_64 &random)
  {
    const auto [a, b] = Operands<Subtract>(size, chain, random);
    const BigInteger left(a, false), right(b, false);
    const BigInteger expected(Reference<Subtract>(a, b), false);
    BigInteger result;
    std::array<double, 11> samples;
    for (double &sample : samples)
    {
      PreparePool(warm);
      sample = TimeUs([&] {
        if constexpr (Subtract)
          result = left - right;
        else
          result = left + right;
      });
      if (result != expected)
        throw std::runtime_error("add/sub holdout result mismatch");
    }
    std::cout << (Subtract ? "sub" : "add") << (chain ? "-chain" : "-balanced") << ','
              << (warm ? "burst" : "parked") << ',' << size << ',' << size << ','
              << std::fixed << std::setprecision(6) << Median(samples) << '\n';
  }

  inline void Holdout()
  {
    if (LimbBits != 64 || ParallelNumThreads() <= 1)
      return;
    const bool haveWarm = WarmPool();
    std::mt19937_64 random(0xD15CA7C4);
    for (SizeT size : {3073u, 5279u, 9461u, 20483u, 36857u, 57347u, 81937u, 155581u,
                       349519u, 786433u, 1572869u, 3145739u})
      for (bool chain : {false, true})
      {
        if (haveWarm)
        {
          HoldoutCase<false>(size, chain, true, random);
          HoldoutCase<true>(size, chain, true, random);
        }
        HoldoutCase<false>(size, chain, false, random);
        HoldoutCase<true>(size, chain, false, random);
      }
  }
}

#endif
