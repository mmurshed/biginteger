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
#include "../algorithms/karatsuba/KaratsubaMultiplication.h"

using namespace BigMath;

int main(int argc, char *argv[])
{
  if (argc < 3)
    return 1;

  if (!freopen(argv[1], "rt", stdin))
    return 1;

  FILE *timeFile = fopen(argv[2], "wt");
  if (!timeFile)
    return 1;

  char op;

  cerr << "Karat Cutoff,Time" << endl;

  BigInteger a, b;

  cin >> a >> op >> b;

  for (SizeT k = 4; k <= 256; k += 2)
  {
    // KaratsubaMultiplication::KARATSUBA_THRESHOLD = k;

    clock_t start = clock();

    vector<DataT> rKarat = KaratsubaMultiplication::Multiply(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());

    clock_t end = clock();

    double timeTakenKarat = (double)(end - start) / CLOCKS_PER_SEC;

    cerr.setf(ios::showpoint);
    cerr << k << "," << timeTakenKarat << endl;
    fprintf(timeFile, "%u,%f\n", k, timeTakenKarat);
  }

  fclose(timeFile);

  return EXIT_SUCCESS;
}
