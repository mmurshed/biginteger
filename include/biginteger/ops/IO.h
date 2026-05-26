/**
 * BigMath: BigInteger stream I/O.
 *
 *   operator<<(ostream, vector<DataT>)  — diagnostic dump of internal limbs
 *   operator<<(ostream, BigInteger)     — decimal text
 *   operator>>(istream, BigInteger)     — parse decimal from streambuf
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_IO
#define BIGINTEGER_IO

#include <iostream>
#include <vector>

#include "../BigInteger.h"
#include "../common/Util.h"
#include "../common/Constants.h"
#include "../common/Parser.h"

namespace BigMath
{
  std::ostream &operator<<(std::ostream &stream, std::vector<DataT> const &out);
  std::ostream &operator<<(std::ostream &stream, BigInteger const &out);
  std::istream &operator>>(std::istream &stream, BigInteger &in);
}

#endif
