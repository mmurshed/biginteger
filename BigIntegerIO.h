/**
 * BigInteger Class
 * Version 8.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGERIO_H
#define BIGINTEGERIO_H

#include <iostream>
using namespace std;

#include "BigIntegerUtil.h"
#include "BigInteger.h"
#include "BigIntegerParser.h"

namespace BigMath
{
// Inserter
  ostream& operator<<(ostream& stream,  BigInteger const& out)
  {
    stream << BigIntegerParser::ToString(out);
    return stream;
  }

  // Extracter
  istream& operator>>(istream& stream, BigInteger& in)
  {
    SizeT SIZE = 10000;
    char *data = new char[SIZE];

    bool isNegative = false;

    stream >> ws;

    int input = stream.get();

    if(input == '-')
      isNegative = true;
    else if(input == '+')
      isNegative = false;
    else
      stream.putback(input);
    
    SizeT len = 0;
    
    while(!stream.eof())
    {
      input = stream.get();
      if(stream.eof())
        break;

      if(isdigit(input))
        data[len++] = input;
      else
      {
        stream.putback(input);
        break;
      }

      // Resize the input buffer
      if(len >= SIZE)
      {
        SIZE += SIZE;
        char *p = new char [SIZE];
        strcpy(p, data);
        delete [] data;
        data = p;
      }
    }
    data[len] = 0;
    
    in = BigIntegerParser::Parse(data);
    
    if(!in.IsZero() && isNegative)
        in.SetSign(isNegative);

    delete [] data;

    return stream;
  }

}

#endif