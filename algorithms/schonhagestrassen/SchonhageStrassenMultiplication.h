
/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef TOOMCOOK_MULTIPLICATION
#define TOOMCOOK_MULTIPLICATION

#include <vector>
#include <stack>
using namespace std;

#include "../../BigInteger.h"
#include "../BigIntegerUtil.h"
#include "../classic/ClassicAddition.h"
#include "../classic/ClassicSubtraction.h"
#include "../classic/ClassicMultiplication.h"
#include "../classic/ClassicDivision.h"

namespace BigMath
{
  class SchonhageStrassenMultiplication
  {
    vector<DataT> Multiply(
      vector<DataT> const& a,
      vector<DataT> const& b,
      ULong base,
      bool trim=true)
    {
      vector<DataT> w = ClassicMultiplication::Multiply(a, b, base);

      if(trim)
      {
        BigIntegerUtil::TrimZeros(w);
      }
      
      return w;

    }
   };
}

#endif
