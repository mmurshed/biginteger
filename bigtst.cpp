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

using namespace BigMath;

#define DEBUG 0

// BigInteger& Factorial(unsigned int n)
// {
//   BigInteger& r = *new BigInteger(1l);
//   for(unsigned int i = 2; i <= n; i++)
//     r = r * i;

//   return r;
// }

int main(int argc, char *argv[])
{
  if(argc < 3)
    return 1;

  if(!freopen(argv[1],"rt",stdin))
    return 1;
  if(!freopen(argv[2],"wt",stdout))
    return 1;

  FILE *timeFile = fopen("time.csv", "wt");
  if(!timeFile)
    return 1;

  long DATA = 1;

  if(DEBUG)
  {
    freopen("test.txt","wt",stdout);
    clock_t start = clock();
    //cout << Factorial(5000) << endl;
    clock_t end = clock();
    cerr.setf(ios::showpoint);
    cerr << (double)(end-start)/CLOCKS_PER_SEC << endl;
  }
  else
  {
    char op;
    bool rp;
    bool un;

    bool first = true;

    clock_t startTotal = clock();

    while(true)
    {
      BigInteger m, n;
      BigInteger r;
      if(!first)
        cout << endl;
      
      first = false;
      cin >> m >> op >> n;
      un = false;

      clock_t start = clock();

      try {
        switch (op)
        {
          case '+': r = m + n; break;
          case '-': r = m - n; break;
          case '*': r = m * n; break;
          case '/': r = m / n; break;
          case '%': r = m % n; break;
          // case '^': r = m.Power(n.toLong());
          case '=': rp = m.CompareTo(n) == 0; un = true; break;
          case '>': rp = m.CompareTo(n) > 0; un = true; break;
          case '<': rp = m.CompareTo(n) < 0; un = true; break;
          default: break;
        }
      } catch(...) {}

      clock_t end = clock();

      cout << "Data : " << DATA << " ";

      if(un) {
        if(rp) cout << "true";
        else cout << "false";
      }
      else cout << r;
      
      double timeTaken = (double)(end-start)/CLOCKS_PER_SEC;
      cerr << "Data : " << DATA << "  Time : " ;      
      
      cerr.setf(ios::showpoint);
      cerr << timeTaken << endl;
      fprintf(timeFile, "%d,%f\n", r.size(), timeTaken);

      DATA++;
    
      if(cin.eof())
        break;
    }

    clock_t endTotal = clock();
    cerr << "Total  Time : " ;
    cerr.setf(ios::showpoint);
    cerr << (double)(endTotal-startTotal)/CLOCKS_PER_SEC << endl;

    fclose(timeFile);

  }

  return EXIT_SUCCESS;
}

// int main(int argc, char *argv[])
// {
//   // BigInteger zero1;
//   // bool z1 = zero1.IsZero(); // should be true

//   // BigInteger zero2(5, false, 0);
//   // bool z2 = zero2.IsZero(); // should be true

//   // vector<DataT> a(zero2.GetInteger());
//   // a[0] = 200;
//   // BigInteger zero3(a, false);
//   // bool z3 = zero3.IsZero(); // should be false

//   // BigInteger zero1;
//   // SizeT z1 = zero1.TrimZeros(); // should be 0

//   // BigInteger zero2(5, false, 0);
//   // vector<DataT> a(zero2.GetInteger());

//   // SizeT z2 = zero2.TrimZeros(); // should be 5

//   // a[0] = 200;
//   // BigInteger zero3(a, false);
//   // SizeT z3 = zero3.TrimZeros(); // should be 4

//   // BigInteger zero1;
//   // SizeT z1 = BigIntegerUtil::FindNonZeroByte(zero1.GetInteger()); // should be 0

//   // BigInteger zero2(5, false, 0);
//   // vector<DataT> a(zero2.GetInteger());

//   // SizeT z2 = BigIntegerUtil::FindNonZeroByte(zero1.GetInteger()); // should be 0

//   // a[0] = 200;
//   // SizeT z3 = BigIntegerUtil::FindNonZeroByte(a); // should be 1

//   BigInteger& a = BigIntegerParser::Parse("12032759530727187026955295177759279936988003947824720956292");
//   BigInteger& b = BigIntegerParser::Parse("552521883657916592106");

//   // BigInteger& a = BigIntegerParser::Parse("14795936");
//   // BigInteger& b = BigIntegerParser::Parse("33202659");

//   // BigInteger& add = a + b;
//   // string& str = BigIntegerParser::ToString(add);
//   // bool cmp = (str == "4898444632449261293846336479321434");

//   // BigInteger& sub = a - b;
//   // str = sub.ToString();
//   // cmp = (str == "123456789012345678901230000000");

//   BigInteger& mult = a * b;
//   string& str = BigIntegerParser::ToString(mult);
//   bool cmp = (str == "6648362961520133879513534129214498836071414231023154552276131181220034018230952");

//   // BigInteger& div = a / b;
//   // string& str = div.ToString();
//   // bool cmp = (str == "44557413417469525");

//   // BigInteger& rem = a % b;
//   // string& str = rem.ToString();
//   // bool cmp = (str == "1494561306735");

//   return 0;
// }
