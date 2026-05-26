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

#include "biginteger/BigInteger.h"
#include "biginteger/ops/IO.h"
#include "biginteger/ops/Operations.h"
#include "biginteger/common/Parser.h"
#include "biginteger/algorithms/multiplication/ToomCookMultiplication.h"

#include "biginteger/algorithms/multiplication/KaratsubaMultiplication.h"

#include "biginteger/algorithms/multiplication/NTTMultiplication.h"

using namespace BigMath;

int main(int argc, char *argv[])
{
  BigInteger a = Parse("123458");
  BigInteger b = Parse("234112");

  vector<DataT> mult = ToomCookMultiplication::Multiply(a.GetInteger(), b.GetInteger(), BigInteger::Base());
  string str = ToString(mult);
  bool cmp = (str == "28902999296");
  cout << "Toom-Cook result: " << str << " Cmp: " << cmp << endl;

  vector<DataT> kara = KaratsubaMultiplication::Multiply(a.GetInteger(), b.GetInteger(), BigInteger::Base());
  string str_kara = ToString(kara);
  bool cmp_kara = (str_kara == "28902999296");
  cout << "Karatsuba result: " << str_kara << " Cmp: " << cmp_kara << endl;

  vector<DataT> ntt_res = NTTMultiplication::Multiply(a.GetInteger(), b.GetInteger(), BigInteger::Base());
  string str_ntt = ToString(ntt_res);
  bool cmp_ntt = (str_ntt == "28902999296");
  cout << "NTT result: " << str_ntt << " Cmp: " << cmp_ntt << endl;

  if (str != "28902999296" || str_kara != "28902999296" || str_ntt != "28902999296") {
    cerr << "Verification FAILED!" << endl;
    return 1;
  }

  cout << "Verification succeeded!" << endl;
  return 0;
}
