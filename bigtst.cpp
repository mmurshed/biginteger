/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Thu Aug 29 02:15:55 BDT 2002
    copyright            : (C) 2002 by Mahbub Murshed Suman
    email                : suman@bttb.net.bd
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <malloc.h>
#include <cmath>
#include <cstring>
#include <ctime>
#include <strstream>
#include <string>
#include <stdexcept>

using namespace std;

#include "BigInteger.h"

using namespace BigMath;

#define DEBUG 0

BigInteger& Factorial(unsigned int n)
{
  BigInteger& r = *new BigInteger(1l);
  for(unsigned int i=2;i<=n;i++)
    r = r.Multiply(i);

  return r;
}

int main(int argc, char *argv[])
{
  if(freopen("test2.in","rt",stdin)==0) exit(1);
  freopen("test2.out","wt",stdout);
  // freopen("test2.err","wt",stderr);
  // ofstream ferr("test2.err");
  // ofstream foul("test2.in2");

  long DATA = 1;
  cerr << BigIntegerAbout << endl;

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
    BigInteger m, n, r;
    char op;
  bool rp;
  bool un;
  string const& sep = *new string("|");
  while(true)
  {
    cin >> m >> op >> n;
	if(cin.eof())
		break;
	/*
	m.DigitSeperatedPrint (fout,sep);
	fout << endl;
	n.DigitSeperatedPrint (fout,sep);
	fout << endl;
	*/
    
	if(op=='=' || op =='<' || op =='>') un = true;
    else un = false;
    // if(m.isZero() && n.isZero()) break;
    clock_t start = clock();
	// cerr = ferr;
    try {
      switch (op)
      {
      case '+': r = m + n; break;
      case '-': r = m - n; break;
      case '*': r = m*n; break;
      case '/': r = m / n; break;
      case '%': r = m%n; break;
      case '^': r = m.Power(n.toLong());
      case '=': rp = m.CompareTo(n) == 0;break;
      case '>': rp = m.CompareTo(n) > 0;break;
      case '<': rp = m.CompareTo(n) < 0;break;
      }
    } catch(...) {}
	// ferr.close ();
	// fclose(stderr);

    cout << "Data : " << DATA << " ";

    // cout << m << ' ' << op << ' ' << n << " = ";

	if(un) {
      if(rp) cout << "true";
      else cout << "false";
    }
    else cout << r;

	cout << endl;
    clock_t end = clock();
    cerr << "Data : " << DATA << "  Time : " ;
    cerr.setf(ios::showpoint);
    cerr << (double)(end-start)/CLOCKS_PER_SEC << endl;
    DATA ++;
  }
  }

  return EXIT_SUCCESS;
}
