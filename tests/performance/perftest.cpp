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

#include "../../BigInteger.h"
#include "../../BigIntegerIO.h"
#include "../../BigIntegerOperations.h"
#include "../../BigIntegerParser.h"
#include "../../BigIntegerComparator.h"

using namespace BigMath;


int main(int argc, char *argv[])
{
  if(argc < 4)
  {
    cerr << "Usage: perftst [INPUT] [ANSWER] [CSVFILE]" << endl;
    return 1;
  }

  if(!freopen(argv[1],"rt",stdin))
    return 1;

  ifstream ansFile(argv[2]);

  FILE *timeFile = fopen(argv[3], "wt");
  if(!timeFile)
    return 1;

  long DATA = 1;

  char op;

  fprintf(timeFile, "Karatsuba,Classical,Toom-Cook,Results Digit");
  fflush(timeFile);

  cerr << "Data,Karatsuba,Classical,Toom-Cook,Results Digit" << endl;

  while(true)
  {
    BigInteger a, b;
   
    cin >> a >> op >> b;

    string line;
    std::getline(ansFile, line);
    BigInteger ans = BigIntegerParser::Parse(line.c_str());

    clock_t start = clock();
    vector<DataT> rKarat = KaratsubaMultiplication::Multiply(
          a.GetInteger(),
          b.GetInteger(),
          BigInteger::Base());
    clock_t end = clock();
    double timeTakenKarat = (double)(end - start) / CLOCKS_PER_SEC;

    Int cmp = BigIntegerComparator::CompareTo(rKarat, ans.GetInteger());
    if(cmp != 0)
    {
      cerr << "Karatsuba algorithm failed." << endl;
      return cmp;
    }
   
    start = clock();
    vector<DataT> rClassical = ClassicMultiplication::Multiply(
          a.GetInteger(),
          b.GetInteger(),
          BigInteger::Base(),
          true);
    end = clock();
    double timeTakenClassical = (double)(end - start) / CLOCKS_PER_SEC;

    cmp = BigIntegerComparator::CompareTo(rClassical, ans.GetInteger());
    if(cmp != 0)
    {
      cerr << "Classical algorithm failed." << endl;
      return cmp;
    }

    start = clock();
    ToomCookMultiplication tcm;
    vector<DataT> rToom = tcm.Multiply(
      a.GetInteger(),
      b.GetInteger(),
      BigInteger::Base());
    end = clock();
    double timeTakenToomCook = (double)(end - start) / CLOCKS_PER_SEC;

    cmp = BigIntegerComparator::CompareTo(rToom, ans.GetInteger());
    if(cmp != 0)
    {
      cerr << "Toom-Cook algorithm failed." << endl;
      return cmp;
    }

    // start = clock();
    // vector<DataT> rToomMem = ToomCookMultiplication::Multiply(
    //   a.GetInteger(),
    //   b.GetInteger(),
    //   BigInteger::Base());
    // end = clock();
    // double timeTakenToomCookMem = (double)(end - start) / CLOCKS_PER_SEC;

    // cmp = BigIntegerComparator::CompareTo(rToomMem, ans.GetInteger());
    // if(cmp != 0)
    // {
    //   cerr << "Toom-Cook-Memory-Optimized algorithm failed." << endl;
    //   return cmp;
    // }

    cerr.setf(ios::showpoint);
    // cerr << DATA << "," << timeTakenKarat << "," << timeTakenClassical << "," << timeTakenToomCook << "," << timeTakenToomCookMem << "," << rKarat.size() << endl;
    cerr << DATA << "," << timeTakenKarat << "," << timeTakenClassical << "," << timeTakenToomCook << "," << rKarat.size() << endl;
    // fprintf(timeFile, "%f,%f,%f,%f,%lu\n", timeTakenKarat, timeTakenClassical, timeTakenToomCook, timeTakenToomCookMem, rKarat.size());
    fprintf(timeFile, "%f,%f,%f,%lu\n", timeTakenKarat, timeTakenClassical, timeTakenToomCook, rKarat.size());
    fflush(timeFile);

    DATA++;
  
    if(cin.eof())
      break;
  }

  fclose(timeFile);

  return EXIT_SUCCESS;
}

