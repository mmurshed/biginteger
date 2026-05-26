/**
 * BigMath: Generic base conversion implementation.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/algorithms/BaseConversion.h"

namespace BigMath
{
  std::vector<DataT> ConvertBase(std::vector<DataT> const &bigIntB1,
                                 SizeT start, SizeT end,
                                 BaseT base1,
                                 BaseT base2)
  {
    if (base1 == base2)
      return std::vector<DataT>(bigIntB1);

    std::vector<DataT> bigIntB2(1, 0);

    for (Int i = (Int)end; i >= (Int)start; i--)
    {
      ClassicMultiplication::MultiplyTo(bigIntB2, base1, base2);
      AddTo(bigIntB2, bigIntB1[i], base2);
    }

    return bigIntB2;
  }

  std::vector<DataT> ConvertBase(std::vector<DataT> const &bigIntB1,
                                 BaseT base1,
                                 BaseT base2)
  {
    return ConvertBase(bigIntB1, 0, bigIntB1.size() - 1, base1, base2);
  }
}
