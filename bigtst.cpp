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

      // cout << "Data : " << DATA << " ";

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

