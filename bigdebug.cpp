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

#include "BigInteger.h"
#include "BigIntegerIO.h"
#include "BigIntegerOperations.h"
#include "BigIntegerParser.h"
#include "ToomCookAlgorithm.h"

using namespace BigMath;

int main(int argc, char *argv[])
{
  // BigInteger zero1;
  // bool z1 = zero1.IsZero(); // should be true

  // BigInteger zero2(5, false, 0);
  // bool z2 = zero2.IsZero(); // should be true

  // vector<DataT> a(zero2.GetInteger());
  // a[0] = 200;
  // BigInteger zero3(a, false);
  // bool z3 = zero3.IsZero(); // should be false

  // BigInteger zero1;
  // SizeT z1 = zero1.TrimZeros(); // should be 0

  // BigInteger zero2(5, false, 0);
  // vector<DataT> a(zero2.GetInteger());

  // SizeT z2 = zero2.TrimZeros(); // should be 5

  // a[0] = 200;
  // BigInteger zero3(a, false);
  // SizeT z3 = zero3.TrimZeros(); // should be 4

  // BigInteger zero1;
  // SizeT z1 = BigIntegerUtil::FindNonZeroByte(zero1.GetInteger()); // should be 0

  // BigInteger zero2(5, false, 0);
  // vector<DataT> a(zero2.GetInteger());

  // SizeT z2 = BigIntegerUtil::FindNonZeroByte(zero1.GetInteger()); // should be 0

  // a[0] = 200;
  // SizeT z3 = BigIntegerUtil::FindNonZeroByte(a); // should be 1

  BigInteger a = BigIntegerParser::Parse("12032759530727187026955295177759279936988003947824720956292");
  BigInteger b = BigIntegerParser::Parse("552521883657916592106");

  // BigInteger& a = BigIntegerParser::Parse("14795936");
  // BigInteger& b = BigIntegerParser::Parse("33202659");

  // BigInteger& add = a + b;
  // string& str = BigIntegerParser::ToString(add);
  // bool cmp = (str == "4898444632449261293846336479321434");

  // BigInteger& sub = a - b;
  // str = sub.ToString();
  // cmp = (str == "123456789012345678901230000000");

  ToomCookAlgorithm tca;
  vector<DataT> mult = tca.MultiplyUnsigned(a.GetInteger(), b.GetInteger(), BigInteger::Base());
  string str = BigIntegerParser::ToString(mult);
  bool cmp = (str == "6648362961520133879513534129214498836071414231023154552276131181220034018230952");

  // BigInteger& div = a / b;
  // string& str = div.ToString();
  // bool cmp = (str == "44557413417469525");

  // BigInteger& rem = a % b;
  // string& str = rem.ToString();
  // bool cmp = (str == "1494561306735");

  return 0;
}