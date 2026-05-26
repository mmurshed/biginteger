/**
 * BigMath: Generic base conversion (Knuth Vol. 2 §4.4).
 *
 * Used by tests and legacy I/O paths. Production decimal I/O routes through
 * common/Parser.h instead, which uses the optimized D&C path with cached
 * powers of 10.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BASE_CONVERSION
#define BASE_CONVERSION

#include <vector>

#include "../common/Util.h"
#include "Addition.h"
#include "multiplication/ClassicMultiplication.h"

namespace BigMath
{
  std::vector<DataT> ConvertBase(std::vector<DataT> const &bigIntB1,
                                 SizeT start, SizeT end,
                                 BaseT base1,
                                 BaseT base2);

  std::vector<DataT> ConvertBase(std::vector<DataT> const &bigIntB1,
                                 BaseT base1,
                                 BaseT base2);
}

#endif
