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


int main(int argc, char *argv[])
{
  if(argc < 3)
    return 1;

  if(!freopen(argv[1],"rt",stdin))
    return 1;

  FILE *timeFile = fopen(argv[2], "wt");
  if(!timeFile)
    return 1;

  long DATA = 1;

  char op;

  cerr << "Data,Karatsuba,Classical,Results Digit" << endl;

  while(true)
  {
    BigInteger a, b;
   
    cin >> a >> op >> b;
    
    clock_t start = clock();

    vector<DataT> rKarat = KaratsubaAlgorithm::MultiplyUnsigned(
          a.GetInteger(),
          b.GetInteger(),
          BigInteger::Base());

    clock_t end = clock();
   
    double timeTakenKarat = (double)(end - start) / CLOCKS_PER_SEC;

    start = clock();

    vector<DataT> rClassical = ClassicalAlgorithms::MultiplyUnsigned(
          a.GetInteger(),
          b.GetInteger(),
          BigInteger::Base());

    end = clock();
   
    double timeTakenClassical = (double)(end - start) / CLOCKS_PER_SEC;

    Int cmp = ClassicalAlgorithms::UnsignedCompareTo(rKarat, rClassical);
    if(cmp != 0)
    {
      cerr << "Result mismatch" << endl;
      return cmp;
    }
   
    cerr.setf(ios::showpoint);
    cerr << DATA << "," << timeTakenKarat << "," << timeTakenClassical << "," << rKarat.size() << endl;
    fprintf(timeFile, "%f,%f,%lu\n", timeTakenKarat, timeTakenClassical, rKarat.size());

    DATA++;
  
    if(cin.eof())
      break;
  }

  fclose(timeFile);

  return EXIT_SUCCESS;
}

