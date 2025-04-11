/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_SHIFT
#define BIGINTEGER_SHIFT

#include <utility>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "../algorithms/classic/CommonAlgorithms.h"

namespace BigMath
{
  BigInteger operator<<(BigInteger const &a, SizeT b)
  {
    vector<DataT> result = CommonAlgorithms::ShiftLeft(a.GetInteger(), b);
    return BigInteger(result, a.IsNegative());
  }

  BigInteger operator>>(BigInteger const &a, SizeT b)
  {
    vector<DataT> result = CommonAlgorithms::ShiftRight(a.GetInteger(), b);
    return BigInteger(result, a.IsNegative());
  }
}

#endif