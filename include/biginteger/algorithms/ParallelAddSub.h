#ifndef BIGMATH_PARALLEL_ADD_SUB
#define BIGMATH_PARALLEL_ADD_SUB

#include <algorithm>
#include <array>
#include <limits>
#include <span>
#include <vector>

#include "../common/LimbKernels.h"
#include "../common/Parallel.h"

namespace BigMath::AddSubDetail
{
  using Limbs = std::span<const DataT>;

  inline SizeT BlockCount(SizeT limbs, SizeT workers)
  {
    constexpr SizeT MinBlockLimbs = 2048;
    return std::max<SizeT>(1, std::min(workers, limbs / MinBlockLimbs));
  }

  template <bool Subtract>
  inline SizeT DispatchBlocks(SizeT a, SizeT b)
  {
    constexpr SizeT warm = Subtract ? BIGMATH_SUB_PARALLEL_WARM_LIMBS : BIGMATH_ADD_PARALLEL_WARM_LIMBS;
    constexpr SizeT parked = Subtract ? BIGMATH_SUB_PARALLEL_PARKED_LIMBS : BIGMATH_ADD_PARALLEL_PARKED_LIMBS;
    if constexpr (BIGMATH_LIMB_64 && BIGMATH_USE_THREADS && (warm != 0 || parked != 0))
    {
      const bool allowParked = parked != 0 && a >= parked;
      if (a == b && (allowParked || (warm != 0 && a >= warm)))
        return BlockCount(a, ParallelAvailableThreads(allowParked));
    }
    return 1;
  }

  template <bool Subtract>
  inline DataT SerialBlock(DataT *out, Limbs a, Limbs b, SizeT first, SizeT end)
  {
    const SizeT shared = std::max(first, std::min(end, static_cast<SizeT>(b.size())));
    DataT flag = 0;
    if (shared != first)
    {
      if constexpr (Subtract)
        flag = LimbKernel::SubN(out + first, a.data() + first, b.data() + first, shared - first);
      else
        flag = LimbKernel::AddN(out + first, a.data() + first, b.data() + first, shared - first);
    }
    if (shared != end)
    {
      if constexpr (Subtract)
        flag = LimbKernel::SubBorrow(out + shared, a.data() + shared, end - shared, flag);
      else
        flag = LimbKernel::AddCarry(out + shared, a.data() + shared, end - shared, flag);
    }
    return flag;
  }

  // out spans a.size() limbs and may alias either input exactly. a is at
  // least as long as b, and the requested block count fits in the pool.
  // Returns the carry or borrow beyond the last limb.
  template <bool Subtract>
  inline DataT ParallelMagnitude(DataT *out, Limbs a, Limbs b, SizeT workers)
  {
    const auto count = static_cast<SizeT>(a.size());
    const SizeT tasks = std::min(count, workers);
    if (tasks <= 1)
      return SerialBlock<Subtract>(out, a, b, 0, count);

    constexpr DataT propagatingLimb = Subtract ? 0 : std::numeric_limits<DataT>::max();
    struct Block
    {
      SizeT first;
      SizeT end;
      SizeT propagateEnd;
      DataT outgoing;
      DataT incoming;
    };
    std::array<Block, BIGMATH_MAX_THREADS> blocks;
    ParallelDo(static_cast<Int>(tasks), [&](Int begin, Int end) {
      for (Int i = begin; i < end; ++i)
      {
        auto &block = blocks[i];
        block.first = static_cast<SizeT>(ULong{count} * static_cast<UInt>(i) / tasks);
        block.end = static_cast<SizeT>(ULong{count} * static_cast<UInt>(i + 1) / tasks);
        block.outgoing = SerialBlock<Subtract>(out, a, b, block.first, block.end);
        block.propagateEnd = block.first;
        // A carry crosses an all-ones result. A borrow crosses an all-zero
        // result. In either case only the prefix can be affected by a fixup.
        if (i != 0)
          while (block.propagateEnd < block.end && out[block.propagateEnd] == propagatingLimb)
            ++block.propagateEnd;
      }
    });

    DataT flag = 0;
    SizeT longFixups = 0;
    for (SizeT i = 0; i < tasks; ++i)
    {
      auto &block = blocks[i];
      block.incoming = flag;
      flag = block.outgoing | (flag & (block.propagateEnd == block.end));
      if (block.incoming && block.propagateEnd - block.first >= 16'384)
        ++longFixups;
    }

    const auto fix = [&](Int begin, Int end) {
      for (Int i = begin; i < end; ++i)
      {
        auto const &block = blocks[i];
        if (!block.incoming)
          continue;
        std::fill(out + block.first, out + block.propagateEnd, ~propagatingLimb);
        if (block.propagateEnd != block.end)
        {
          if constexpr (Subtract)
            --out[block.propagateEnd];
          else
            ++out[block.propagateEnd];
        }
      }
    };
    // Random inputs need at most a few limb writes here. Only dispatch a
    // second time if multiple blocks have enough fixup work to share.
    if (longFixups > 1)
      ParallelDo(static_cast<Int>(tasks), fix);
    else
      fix(0, static_cast<Int>(tasks));
    return flag;
  }

  template <bool Subtract>
  inline std::vector<DataT> ParallelResult(Limbs a, Limbs b, SizeT workers)
  {
    std::vector<DataT> result(a.size() + 1);
    const DataT flag = ParallelMagnitude<Subtract>(result.data(), a, b, workers);
    if constexpr (!Subtract)
      result.back() = flag;
    TrimZeros(result);
    return result;
  }

}

#endif
