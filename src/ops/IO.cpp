/**
 * BigMath: BigInteger stream I/O implementation.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/ops/IO.h"
#include "biginteger/algorithms/Addition.h"
#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"

#include <iomanip>

namespace BigMath
{
  // Limb-level diagnostic dump: pretty-prints the raw limb vector in a fixed-width
  // format. Useful for debugging; not the canonical text form (see ToString).
  std::ostream &operator<<(std::ostream &stream, std::vector<DataT> const &out)
  {
    if (out.empty())
    {
      stream << "0";
      return stream;
    }

    int pad = Base100_Zeroes;

    auto it = out.rbegin();
    stream << *it;
    ++it;

    for (; it != out.rend(); ++it)
      stream << std::setw(pad) << std::setfill('0') << *it;

    return stream;
  }

  std::ostream &operator<<(std::ostream &stream, BigInteger const &out)
  {
    stream << ToString(out);
    return stream;
  }

  // Fused extract + parse. Reads digits directly from the streambuf and folds 18 at
  // a time into a base-2^32 limb vector. Avoids materializing the full digit string
  // and the second scan that Parse would do on a captured std::string.
  std::istream &operator>>(std::istream &stream, BigInteger &in)
  {
    bool isNegative = false;

    stream >> std::ws;

    auto *buf = stream.rdbuf();
    int c = buf->sgetc();
    if (c == '-')
    {
      isNegative = true;
      buf->sbumpc();
      c = buf->sgetc();
    }
    else if (c == '+')
    {
      buf->sbumpc();
      c = buf->sgetc();
    }

    while (c == '0')
    {
      buf->sbumpc();
      c = buf->sgetc();
    }

    std::vector<DataT> r;
    r.reserve(64);
    r.push_back(0);

    ULong chunk = 0;
    ULong chunkMul = 1;
    bool anyDigit = false;

    while (c >= '0' && c <= '9')
    {
      chunk = chunk * 10 + (ULong)(c - '0');
      chunkMul *= 10;
      anyDigit = true;
      buf->sbumpc();
      c = buf->sgetc();

      if (chunkMul == Base10_18)
      {
        ClassicMultiplication::MultiplyTo(r, Base10_18, Base2_32);
        AddTo(r, chunk, Base2_32);
        chunk = 0;
        chunkMul = 1;
      }
    }

    if (chunkMul > 1)
    {
      ClassicMultiplication::MultiplyTo(r, chunkMul, Base2_32);
      AddTo(r, chunk, Base2_32);
    }

    if (c == EOF)
      stream.setstate(std::ios::eofbit);

    in = anyDigit ? BigInteger(r, isNegative) : BigInteger();
    return stream;
  }
}
