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
#include <ctime>
#include <string>

using namespace std;

#include "biginteger/BigInteger.h"
#include "biginteger/ops/IO.h"
#include "biginteger/common/Parser.h"
#include "biginteger/common/Comparator.h"
#include "biginteger/algorithms/division/BurnikelZieglerDivision.h"
#include "biginteger/algorithms/division/FastDivision.h"
#include "biginteger/algorithms/division/KnuthDivision.h"

using namespace BigMath;

int main(int argc, char *argv[])
{
  if (!freopen("input.txt", "rt", stdin))
    return 1;

  if (!freopen("output.txt", "wt", stdout))
    return 1;

  ifstream ansFile("answer.txt");

  FILE *timeFile = fopen("result.csv", "wt");
  if (!timeFile)
    return 1;

  long DATA = 1;
  char op;

  fprintf(timeFile, "Results Digit,Knuth,Fast,BurnikelZiegler\n");
  fflush(timeFile);

  cerr << "Data,Results Digit,Knuth,Fast,BurnikelZiegler" << endl;

  while (!cin.eof())
  {
    BigInteger a, b;
    cin >> a >> op >> b;
    if (cin.eof())
      break;

    string line;
    getline(ansFile, line);
    BigInteger ansq = Parse(line.c_str());
    getline(ansFile, line);
    BigInteger ansr = Parse(line.c_str());

    clock_t start, end;
    vector<DataT> q, r;
    double timeTakenKnuth = 0.;
    double timeTakenFast = 0.;
    double timeTakenBZ = 0.;

    start = clock();
    tie(q, r) = KnuthDivision::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());
    end = clock();
    timeTakenKnuth = (double)(end - start) / CLOCKS_PER_SEC;

    if (Compare(q, ansq.GetInteger()) != 0 || Compare(r, ansr.GetInteger()) != 0)
      cerr << "Knuth algorithm failed." << endl;

    start = clock();
    tie(q, r) = FastDivision::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());
    end = clock();
    timeTakenFast = (double)(end - start) / CLOCKS_PER_SEC;

    if (Compare(q, ansq.GetInteger()) != 0 || Compare(r, ansr.GetInteger()) != 0)
      cerr << "Fast division algorithm failed." << endl;

    start = clock();
    tie(q, r) = BurnikelZieglerDivision::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());
    end = clock();
    timeTakenBZ = (double)(end - start) / CLOCKS_PER_SEC;

    if (Compare(q, ansq.GetInteger()) != 0 || Compare(r, ansr.GetInteger()) != 0)
      cerr << "Burnikel-Ziegler division algorithm failed." << endl;

    cerr.setf(ios::showpoint);
    cerr << DATA << "," << ansq.size() << "," << timeTakenKnuth << "," << timeTakenFast << "," << timeTakenBZ << endl;
    fprintf(timeFile, "%u,%f,%f,%f\n", a.size(), timeTakenKnuth, timeTakenFast, timeTakenBZ);
    fflush(timeFile);

    DATA++;
  }

  fclose(timeFile);

  return EXIT_SUCCESS;
}
