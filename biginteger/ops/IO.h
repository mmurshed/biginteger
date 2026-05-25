/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_IO
#define BIGINTEGER_IO

#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

#include "../BigInteger.h"
#include "../common/Util.h"
#include "../common/Constants.h"
#include "../common/Parser.h"

namespace BigMath
{
  ostream &operator<<(ostream &stream, const vector<DataT> &out)
  {
    if (out.empty())
    {
      stream << "0";
      return stream;
    }

    int pad = Base100_Zeroes;

    // Print in reverse order (most-significant digit first)
    auto it = out.rbegin();

    // First digit: print without leading zeros (optional – if you want to avoid 0005 printing as 0005 at the front)
    stream << *it;
    ++it;

    // Remaining digits: pad to 4 digits with leading zeros
    for (; it != out.rend(); ++it)
    {
      stream << setw(pad) << setfill('0') << *it;
    }

    return stream;
  }

  // Inserter
  ostream &operator<<(ostream &stream, BigInteger const &out)
  {
    stream << ToString(out);
    return stream;
  }

  // Extracter — fused extract + parse. Reads digits directly from the streambuf
  // and folds 18 at a time into a base-2^32 limb vector. Avoids materializing
  // the full digit string and the second scan that Parse would do.
  istream &operator>>(istream &stream, BigInteger &in)
  {
    bool isNegative = false;

    stream >> ws;

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

    // Skip leading zeros.
    while (c == '0')
    {
      buf->sbumpc();
      c = buf->sgetc();
    }

    vector<DataT> r;
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
      stream.setstate(ios::eofbit);

    if (!anyDigit)
    {
      in = BigInteger();
    }
    else
    {
      in = BigInteger(r, isNegative);
    }

    return stream;
  }

}

#endif