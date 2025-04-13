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
#include "../../ops/BigIntegerIO.h"
#include "../../ops/BigIntegerOperations.h"
#include "../../BigIntegerParser.h"
#include "../../BigIntegerComparator.h"
#include "../../algorithms/divison/KnuthDivision.h"
#include "../../algorithms/divison/NewtonRaphsonDivision.h"

using namespace BigMath;

int main(int argc, char *argv[])
{
  /*
  if(argc < 5)
  {
    cerr << "Usage: perftst [INPUT] [OUTPUT] [ANSWER] [CSVFILE]" << endl;
    return 1;
  }

  if(!freopen(argv[1],"rt",stdin))
    return 1;

  if(!freopen(argv[2],"wt",stdout))
    return 1;

  ifstream ansFile(argv[3]);

  FILE *timeFile = fopen(argv[4], "wt");
  if(!timeFile)
    return 1;
*/

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

  fprintf(timeFile, "Results Digit,Classic,Newton-Raphson,Newton-Raphson2\n");
  fflush(timeFile);

  cerr << "Data,Results Digit,Classic,Newton-Raphson,Newton-Raphson2" << endl;

  while (true)
  {
    BigInteger a, b;

    cin >> a >> op >> b;

    string line;
    getline(ansFile, line);
    BigInteger ansq = BigIntegerParser::Parse(line.c_str());
    getline(ansFile, line);
    BigInteger ansr = BigIntegerParser::Parse(line.c_str());

    clock_t start = clock();
    auto [q, r] = KnuthDivision::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());
    clock_t end = clock();
    double timeTakenClassic = (double)(end - start) / CLOCKS_PER_SEC;

    int cmpq = BigIntegerComparator::Compare(q, ansq.GetInteger());
    int cmpr = BigIntegerComparator::Compare(r, ansr.GetInteger());
    if (cmpq != 0 || cmpr != 0)
    {
      cerr << q << endl;
      cerr << r << endl;
      cerr << "Classic algorithm failed." << endl;
    }

    start = clock();
    tie(q, r) = NewtonRaphsonDivision::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());
    end = clock();
    double timeTakenNR = (double)(end - start) / CLOCKS_PER_SEC;

    cmpq = BigIntegerComparator::Compare(q, ansq.GetInteger());
    cmpr = BigIntegerComparator::Compare(r, ansr.GetInteger());
    if (cmpq != 0 || cmpr != 0)
    {
      cerr << q << endl;
      cerr << r << endl;
      cerr << "NewtonRaphsonDivision algorithm failed." << endl;
    }


    cerr.setf(ios::showpoint);
    cerr << DATA << "," << q.size() << "," << timeTakenClassic << "," << timeTakenNR << endl;
    fprintf(timeFile, "%lu,%f,%f\n", q.size(), timeTakenClassic, timeTakenNR);
    fflush(timeFile);

    DATA++;

    if (cin.eof())
      break;
  }

  fclose(timeFile);

  return EXIT_SUCCESS;
}
