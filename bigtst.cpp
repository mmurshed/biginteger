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
#include "KaratsubaAlgorithm.h"
#include "BigIntegerOperations.h"
#include "BigIntegerParser.h"

using namespace BigMath;

#define DEBUG 0

BigInteger& Factorial(unsigned int n)
{
  BigInteger& r = *new BigInteger(1l);
  for(unsigned int i = 2; i <= n; i++)
    r = r * i;

  return r;
}

int main(int argc, char *argv[])
{
  if(argc < 3)
    return 1;

  if(!freopen(argv[1],"rt",stdin))
    return 1;
  if(!freopen(argv[2],"wt",stdout))
    return 1;

  long DATA = 1;

  if(DEBUG)
  {
    freopen("test.txt","wt",stdout);
    clock_t start = clock();
    cout << Factorial(5000) << endl;
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

    while(true)
    {
      BigInteger m, n, r;
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
          // case '*': r = m * n; break;
          case '*': r = KaratsubaAlgorithm::Multiply(m, n); break;
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
      
      cerr << "Data : " << DATA << "  Time : " ;
      
      cerr.setf(ios::showpoint);
      cerr << (double)(end-start)/CLOCKS_PER_SEC << endl;

      DATA++;
    
      if(cin.eof())
        break;
    }
  }

  return EXIT_SUCCESS;
}


// int main(int argc, char *argv[])
// {
//   BigInteger a = BigIntegerParser::Parse("147959365715964420810928396835147959365715964420810928396835");
//   BigInteger b = BigIntegerParser::Parse("3320645306084332064530608433206453060843320645306084");

//   // BigInteger a = BigIntegerParser::Parse("14795936");
//   // BigInteger b = BigIntegerParser::Parse("33202");


//   // BigInteger& add = a + b;
//   // string& str = BigIntegerParser::ToString(add);
//   // bool cmp = (str == "147959365715964424131573702919");

//   // BigInteger& sub = a - b;
//   // str = sub.ToString();
//   // cmp = (str == "123456789012345678901230000000");

//   BigInteger& mult = KaratsubaAlgorithm2::Multiply(a, b);
//   string& str = BigIntegerParser::ToString(mult);
//   bool cmp = (str == "491320573255932302006284699031057940887298988964037639872639444167514409048506171355173608877547200365991844140");

//   // BigInteger& div = a / b;
//   // string& str = div.ToString();
//   // bool cmp = (str == "44557413417469525");

//   // BigInteger& rem = a % b;
//   // string& str = rem.ToString();
//   // bool cmp = (str == "1494561306735");
// }
