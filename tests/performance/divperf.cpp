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

#include "../../biginteger/BigInteger.h"
#include "../../biginteger/ops/IO.h"
#include "../../biginteger/ops/Operations.h"
#include "../../biginteger/common/Parser.h"
#include "../../biginteger/common/Comparator.h"
#include "../../biginteger/algorithms/division/ClassicDivision.h"
#include "../../biginteger/algorithms/division/KnuthDivision.h"
#include "../../biginteger/algorithms/division/NewtonRaphsonDivision.h"
#include "../../biginteger/algorithms/division/MontgomeryDivision.h"
#include "../../biginteger/algorithms/division/BarrettDivision.h"

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
    if (cin.eof())
      break;

    BigInteger a, b;

    cin >> a >> op >> b;

    string line;
    getline(ansFile, line);
    BigInteger ansq = Parse(line.c_str());
    getline(ansFile, line);
    BigInteger ansr = Parse(line.c_str());

    clock_t start, end;
    vector<DataT> q, r;
    double timeTakenClassic = 0., timeTakenKnuth = 0., timeTakenNR = 0., timeTakenM = 0., timeTakenB = 0.;
    int cmpq, cmpr;

    /*
    start = clock();
    tie(q, r) = ClassicDivision::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());
    end = clock();
    timeTakenClassic = (double)(end - start) / CLOCKS_PER_SEC;

    cmpq = Compare(q, ansq.GetInteger());
    cmpr = Compare(r, ansr.GetInteger());
    if (cmpq != 0 || cmpr != 0)
    {
      cerr << "Classic algorithm failed." << endl;
      cerr << "a: " << a << endl;
      cerr << "b: " << b << endl;
      if (cmpq != 0)
      {
        cerr << "q: " << q << endl;
        cerr << "correct q: " << ansq << endl;
      }
      if (cmpr != 0)
      {
        cerr << "r: " << r << endl;
        cerr << "correct r: " << ansr << endl;
      }
    }
    */

    start = clock();
    tie(q, r) = KnuthDivision::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());
    end = clock();
    timeTakenKnuth = (double)(end - start) / CLOCKS_PER_SEC;

    cmpq = Compare(q, ansq.GetInteger());
    cmpr = Compare(r, ansr.GetInteger());
    if (cmpq != 0 || cmpr != 0)
    {
      cerr << "Knuth algorithm failed." << endl;
      cerr << "a: " << a << endl;
      cerr << "b: " << b << endl;
      if (cmpq != 0)
      {
        cerr << "q: " << q << endl;
        cerr << "correct q: " << ansq << endl;
      }
      if (cmpr != 0)
      {
        cerr << "r: " << r << endl;
        cerr << "correct r: " << ansr << endl;
      }
    }

    /*
    start = clock();
    tie(q, r) = NewtonRaphsonDivision::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());
    end = clock();
    double timeTakenNR = (double)(end - start) / CLOCKS_PER_SEC;

    cmpq = Compare(q, ansq.GetInteger());
    cmpr = Compare(r, ansr.GetInteger());
    if (cmpq != 0 || cmpr != 0)
    {
      cerr << "Newton Raphson algorithm failed." << endl;
      cerr << "a: " << a << endl;
      cerr << "b: " << b << endl;
      if (cmpq != 0)
      {
        cerr << "q: " << q << endl;
        cerr << "correct q: " << ansq << endl;
      }
      if (cmpr != 0)
      {
        cerr << "r: " << r << endl;
        cerr << "correct r: " << ansr << endl;
      }
    }
    */

    start = clock();
    tie(q, r) = MontgomeryDivision::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());
    end = clock();
    timeTakenM = (double)(end - start) / CLOCKS_PER_SEC;

    cmpq = Compare(q, ansq.GetInteger());
    cmpr = Compare(r, ansr.GetInteger());
    if (cmpq != 0 || cmpr != 0)
    {
      cerr << "Montgomery algorithm failed." << endl;
      cerr << "a: " << a << endl;
      cerr << "b: " << b << endl;
      if (cmpq != 0)
      {
        cerr << "q: " << q << endl;
        cerr << "correct q: " << ansq << endl;
      }
      if (cmpr != 0)
      {
        cerr << "r: " << r << endl;
        cerr << "correct r: " << ansr << endl;
      }
    }

    start = clock();
    tie(q, r) = BarrettDivision::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());
    end = clock();
    timeTakenB = (double)(end - start) / CLOCKS_PER_SEC;

    cmpq = Compare(q, ansq.GetInteger());
    cmpr = Compare(r, ansr.GetInteger());
    if (cmpq != 0 || cmpr != 0)
    {
      cerr << "Barrett algorithm failed." << endl;
      cerr << "a: " << a << endl;
      cerr << "b: " << b << endl;
      if (cmpq != 0)
      {
        cerr << "q: " << q << endl;
        cerr << "correct q: " << ansq << endl;
      }
      if (cmpr != 0)
      {
        cerr << "r: " << r << endl;
        cerr << "correct r: " << ansr << endl;
      }
    }

    cerr.setf(ios::showpoint);
    cerr << DATA << "," << ansq.size() << "," << timeTakenClassic << "," << timeTakenKnuth << "," << timeTakenNR << "," << timeTakenM << "," << timeTakenB << endl;
    fprintf(timeFile, "%lu,%f,%f,%f,%f,%f\n", ansq.size(), timeTakenClassic, timeTakenKnuth, timeTakenNR, timeTakenM, timeTakenB);
    fflush(timeFile);

    DATA++;
  }

  fclose(timeFile);

  return EXIT_SUCCESS;
}
