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
#include "../BigIntegerUtil.h"
#include "../BigIntegerParser.h"

namespace BigMath
{
// Inserter for vector
#include <iomanip> // for setw and setfill

  ostream &operator<<(ostream &stream, const vector<DataT> &out)
  {
    if (out.empty())
    {
      stream << "0";
      return stream;
    }

    int pad = BigIntegerUtil::Base100_Zeroes; // Number of digits to pad with leading zeros

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
    stream << BigIntegerParser::ToString(out);
    return stream;
  }

  // Extracter
  istream &operator>>(istream &stream, BigInteger &in)
  {
    SizeT SIZE = 10000;
    char *data = new char[SIZE];

    bool isNegative = false;

    stream >> ws;

    int input = stream.get();

    if (input == '-')
      isNegative = true;
    else if (input == '+')
      isNegative = false;
    else
      stream.putback(input);

    SizeT len = 0;

    while (!stream.eof())
    {
      input = stream.get();
      if (stream.eof())
        break;

      if (isdigit(input))
        data[len++] = input;
      else
      {
        stream.putback(input);
        break;
      }

      // Resize the input buffer
      if (len >= SIZE)
      {
        SIZE += SIZE;
        char *p = new char[SIZE];
        strncpy(p, data, SIZE);
        delete[] data;
        data = p;
      }
    }
    data[len] = 0;

    in = BigIntegerParser::Parse(data);

    if (!in.IsZero() && isNegative)
      in.SetSign(isNegative);

    delete[] data;

    return stream;
  }

}

#endif