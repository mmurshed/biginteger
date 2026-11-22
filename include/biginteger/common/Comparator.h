/**
 * BigMath: Limb-vector comparison primitives.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_COMPARATOR
#define BIGINTEGER_COMPARATOR

#include <span>
#include <vector>

#include "Constants.h"
#include "Util.h"

namespace BigMath
{
  namespace CompareDetail
  {
    inline constexpr std::size_t BlockLimbs = 16;

    inline Int CompareScalarRange(const DataT *a, const DataT *b, std::size_t size)
    {
      while (size > 0)
      {
        --size;
        if (a[size] != b[size])
          return a[size] > b[size] ? 1 : -1;
      }
      return 0;
    }

    inline Int CompareBlocks(const DataT *a, const DataT *b, std::size_t size)
    {
      if (a == b)
        return 0;
      // Reduce the differences in a fixed-size block so the compiler can
      // vectorize the scan. A differing block is then ordered by whole limbs.
      while (size >= BlockLimbs)
      {
        size -= BlockLimbs;
        DataT differences = 0;
        for (std::size_t i = 0; i < BlockLimbs; ++i)
          differences |= a[size + i] ^ b[size + i];
        if (differences != 0)
          return CompareScalarRange(a + size, b + size, BlockLimbs);
      }
      return CompareScalarRange(a, b, size);
    }

    inline Int CompareEqualLength(const DataT *a, const DataT *b, std::size_t size)
    {
      if (size == 0)
        return 0;
      --size;
      if (a[size] != b[size])
        return a[size] > b[size] ? 1 : -1;
      if (size == 0)
        return 0;
      --size;
      if (a[size] != b[size])
        return a[size] > b[size] ? 1 : -1;
      if (size >= BlockLimbs)
        return CompareBlocks(a, b, size);
      return CompareScalarRange(a, b, size);
    }

    inline bool HasNonzeroHigh(std::span<const DataT> high)
    {
      auto size = high.size();
      while (size >= BlockLimbs)
      {
        size -= BlockLimbs;
        DataT nonzero = 0;
        for (std::size_t i = 0; i < BlockLimbs; ++i)
          nonzero |= high[size + i];
        if (nonzero != 0)
          return true;
      }
      return !IsZero(high.first(size));
    }
  }

  // Spans may include high zero padding, even for otherwise equal values.
  inline Int Compare(std::span<const DataT> a, std::span<const DataT> b)
  {
    const auto aSize = a.size();
    const auto bSize = b.size();
    if (aSize == bSize)
      return CompareDetail::CompareEqualLength(a.data(), b.data(), aSize);
    if (aSize > bSize)
    {
      if (a.back() != 0 || CompareDetail::HasNonzeroHigh(a.subspan(bSize, aSize - bSize - 1)))
        return 1;
      return CompareDetail::CompareEqualLength(a.data(), b.data(), bSize);
    }
    if (b.back() != 0 || CompareDetail::HasNonzeroHigh(b.subspan(aSize, bSize - aSize - 1)))
      return -1;
    return CompareDetail::CompareEqualLength(a.data(), b.data(), aSize);
  }

  inline Int Compare(std::vector<DataT> const &a, DataT b)
  {
    if (a.size() > 1 &&
        (a.back() != 0 || (a.size() > 2 && CompareDetail::HasNonzeroHigh(
            std::span<const DataT>(a).subspan(1, a.size() - 2)))))
      return 1;
    const DataT low = a.empty() ? 0 : a.front();
    if (low == b)
      return 0;
    return low > b ? 1 : -1;
  }

  // Unsigned magnitude compare. Returns 0/+/− for ==, >, <.
  // O(n) time, O(1) space. High zero padding does not affect the result.
  inline Int Compare(std::vector<DataT> const &a, std::vector<DataT> const &b)
  {
    return Compare(std::span<const DataT>(a.data(), a.size()),
                   std::span<const DataT>(b.data(), b.size()));
  }
}

#endif
