/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Thu Aug 29 02:15:55 BDT 2002
    copyright            : (C) 2002 by S. M. Mahbub Murshed
    email                : murshed@gmail.com
 ***************************************************************************/
 
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <cmath>
#include <cstring>
#include <ctime>
#include <strstream>
#include <string>
#include <stdexcept>

using namespace std;

#include "../BigInteger.h"
#include "../BigIntegerIO.h"
#include "../BigIntegerOperations.h"
#include "../BigIntegerParser.h"
#include "../algorithms/toomcookalt/ToomCookMultiplicationAlt.h"

using namespace BigMath;

int main(int argc, char *argv[])
{
  BigInteger a = BigIntegerParser::Parse("10");
  BigInteger b = BigIntegerParser::Parse("728");

  ToomCookMultiplicationAlt tcm;

  vector<DataT> mult = tcm.Multiply(a.GetInteger(), b.GetInteger(), BigInteger::Base());
  string str = BigIntegerParser::ToString(mult);
  bool cmp = (str == "7280");

  return 0;
}
