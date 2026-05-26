/**
 * BigMath: Multi-precision unsigned addition primitives.
 *
 * All operations are on unsigned magnitudes; sign handling lives in the
 * BigInteger operators on top.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef ADDITION
#define ADDITION

#include <vector>

#include "../common/Util.h"

namespace BigMath
{
  // Add scalar b into a[aStart..aEnd], extending a if the carry exceeds the range.
  void AddTo(std::vector<DataT> &a, SizeT aStart, SizeT aEnd,
             ULong b, BaseT base);

  // result[rStart..] = a[aStart..aEnd] + b[bStart..bEnd]. Caller pre-sizes result.
  void Add(std::vector<DataT> const &a, SizeT aStart, SizeT aEnd,
           std::vector<DataT> const &b, SizeT bStart, SizeT bEnd,
           std::vector<DataT> &result, SizeT rStart,
           BaseT base);

  // a += b (whole-vector in-place).
  void AddTo(std::vector<DataT> &a, std::vector<DataT> const &b, BaseT base);

  // a[aStart..aEnd] += b[bStart..bEnd] (range-in-place).
  void AddTo(std::vector<DataT> &a, SizeT aStart, SizeT aEnd,
             std::vector<DataT> const &b, SizeT bStart, SizeT bEnd,
             BaseT base);

  // result = a + b (allocates the result vector).
  std::vector<DataT> Add(std::vector<DataT> const &a,
                         std::vector<DataT> const &b,
                         BaseT base);

  // a += b (scalar in-place); a may grow if b's contribution carries past a's top.
  void AddTo(std::vector<DataT> &a, ULong b, BaseT base);
}

#endif
