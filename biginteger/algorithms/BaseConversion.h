/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BASE_CONVERSION
#define BASE_CONVERSION

#include <vector>
#include <string>
using namespace std;

#include "../common/Util.h"
#include "Addition.h"
#include "multiplication/ClassicMultiplication.h"
#include "division/KnuthDivision.h"

namespace BigMath
{
  // Convert an integer from base1 to base2
  // Donald E. Knuth 4.4
  vector<DataT> ConvertBase(
      vector<DataT> const &bigIntB1, SizeT start, SizeT end,
      BaseT base1,
      BaseT base2)
  {
    if (base1 == base2)
    {
      vector<DataT> bigIntB2(bigIntB1);
      return bigIntB2;
    }

    vector<DataT> bigIntB2(1, 0);

    for (Int i = (Int)end; i >= (Int)start; i--)
    {
      ClassicMultiplication::MultiplyTo(
          bigIntB2,
          base1,
          base2);
      AddTo(
          bigIntB2,
          bigIntB1[i],
          base2);
    }

    return bigIntB2;
  }

  vector<DataT> ConvertBase(
      vector<DataT> const &bigIntB1,
      BaseT base1,
      BaseT base2)
  {
    return ConvertBase(
        bigIntB1, 0, bigIntB1.size() - 1,
        base1,
        base2);
  }

}

#endif
