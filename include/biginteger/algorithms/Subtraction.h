/**
 * BigMath: Multi-precision unsigned subtraction primitives.
 *
 * All operations are on unsigned magnitudes; sign handling lives in the
 * BigInteger operators on top. Subtract assumes a ≥ b.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef SUBTRACTION
#define SUBTRACTION

#include <vector>

#include "../common/Util.h"

namespace BigMath
{
  void SubtractFrom(std::vector<DataT> &a, SizeT aStart, SizeT aEnd,
                    ULong b, BaseT base);

  void SubtractFrom(std::vector<DataT> &a, ULong b, BaseT base);

  // result[rStart..] = a[aStart..aEnd] - b[bStart..bEnd]. Assumes a >= b.
  void Subtract(std::vector<DataT> const &a, SizeT aStart, SizeT aEnd,
                std::vector<DataT> const &b, SizeT bStart, SizeT bEnd,
                std::vector<DataT> &result, SizeT rStart,
                BaseT base);

  void SubtractFrom(std::vector<DataT> &a, std::vector<DataT> const &b, BaseT base);

  void SubtractFrom(std::vector<DataT> &a, SizeT aStart, SizeT aEnd,
                    std::vector<DataT> const &b, SizeT bStart, SizeT bEnd,
                    BaseT base);

  // result = a - b (allocates).
  std::vector<DataT> Subtract(std::vector<DataT> const &a,
                              std::vector<DataT> const &b,
                              BaseT base);

  // Target-buffer overload: caller owns `result`. Reuses storage when large enough.
  void Subtract(std::vector<DataT> &result,
                std::vector<DataT> const &a,
                std::vector<DataT> const &b,
                BaseT base);
}

#endif
