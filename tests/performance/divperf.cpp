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
#include "../../algorithms/division/ClassicDivision.h"
#include "../../algorithms/division/KnuthDivision.h"
#include "../../algorithms/division/NewtonRaphsonDivision.h"
#include "../../algorithms/division/MontgomeryDivision.h"
#include "../../algorithms/division/BarrettDivision.h"

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

  fprintf(timeFile, "Results Digit,Classic,Knuth,Newton-Raphson,Montgomery,Barrett\n");
  fflush(timeFile);

  cerr << "Data,Results Digit,Classic,Knuth,Newton-Raphson,Montgomery,Barrett" << endl;

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
    auto [q, r] = ClassicDivision::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());
    clock_t end = clock();
    double timeTakenClassic = (double)(end - start) / CLOCKS_PER_SEC;

    int cmpq = BigIntegerComparator::Compare(q, ansq.GetInteger());
    int cmpr = BigIntegerComparator::Compare(r, ansr.GetInteger());
    if (cmpq != 0 || cmpr != 0)
    {
      cerr << "a: " << a << endl;
      cerr << "b: " << b << endl;
      cerr << "q: " << q << endl;
      cerr << "r: " << r << endl;
      cerr << "correct q: " << ansq << endl;
      cerr << "correct r: " << ansr << endl;
      cerr << "Classic algorithm failed." << endl;
    }

    start = clock();
    tie(q, r) = KnuthDivision::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());
    end = clock();
    double timeTakenKnuth = (double)(end - start) / CLOCKS_PER_SEC;

    cmpq = BigIntegerComparator::Compare(q, ansq.GetInteger());
    cmpr = BigIntegerComparator::Compare(r, ansr.GetInteger());
    if (cmpq != 0 || cmpr != 0)
    {
      cerr << "a: " << a << endl;
      cerr << "b: " << b << endl;
      cerr << "q: " << q << endl;
      cerr << "r: " << r << endl;
      cerr << "correct q: " << ansq << endl;
      cerr << "correct r: " << ansr << endl;
      cerr << "Knuth algorithm failed." << endl;
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
      cerr << "a: " << a << endl;
      cerr << "b: " << b << endl;
      cerr << "q: " << q << endl;
      cerr << "r: " << r << endl;
      cerr << "correct q: " << ansq << endl;
      cerr << "correct r: " << ansr << endl;
      cerr << "NewtonRaphsonDivision algorithm failed." << endl;
    }

    start = clock();
    tie(q, r) = MontgomeryDivision::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());
    end = clock();
    double timeTakenM = (double)(end - start) / CLOCKS_PER_SEC;

    cmpq = BigIntegerComparator::Compare(q, ansq.GetInteger());
    cmpr = BigIntegerComparator::Compare(r, ansr.GetInteger());
    if (cmpq != 0 || cmpr != 0)
    {
      cerr << "a: " << a << endl;
      cerr << "b: " << b << endl;
      cerr << "q: " << q << endl;
      cerr << "r: " << r << endl;
      cerr << "correct q: " << ansq << endl;
      cerr << "correct r: " << ansr << endl;
      cerr << "Montgomery algorithm failed." << endl;
    }


    start = clock();
    tie(q, r) = BarrettDivision::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());
    end = clock();
    double timeTakenB = (double)(end - start) / CLOCKS_PER_SEC;

    cmpq = BigIntegerComparator::Compare(q, ansq.GetInteger());
    cmpr = BigIntegerComparator::Compare(r, ansr.GetInteger());
    if (cmpq != 0 || cmpr != 0)
    {
      cerr << "a: " << a << endl;
      cerr << "b: " << b << endl;
      cerr << "q: " << q << endl;
      cerr << "r: " << r << endl;
      cerr << "correct q: " << ansq << endl;
      cerr << "correct r: " << ansr << endl;
      cerr << "Barrett algorithm failed." << endl;
    }

    cerr.setf(ios::showpoint);
    cerr << DATA << "," << ansq.size() << "," << timeTakenClassic << "," << timeTakenKnuth << "," << timeTakenNR << "," << timeTakenM << "," << timeTakenB << endl;
    fprintf(timeFile, "%lu,%f,%f,%f,%f,%f\n", ansq.size(), timeTakenClassic, timeTakenKnuth, timeTakenNR, timeTakenM, timeTakenB);
    fflush(timeFile);

    DATA++;

    if (cin.eof())
      break;
  }

  fclose(timeFile);

  return EXIT_SUCCESS;
}
