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
#include "../../ops/BigIntegerOperations.h"
#include "../../BigIntegerParser.h"
#include "../../BigIntegerComparator.h"
#include "../../algorithms/toomcook/ToomCookMultiplication.h"
#include "../../algorithms/toomcook/ToomCookMultiplication2.h"
#include "../../algorithms/toomcook/ToomCookMultiplication2a.h"
#include "../../algorithms/stonehagestrassen/FFTMultiplication.h"

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

  if(!freopen("input.txt","rt",stdin))
  return 1;

  if(!freopen("output.txt","wt",stdout))
  return 1;

  ifstream ansFile("answer.txt");

  FILE *timeFile = fopen("result.csv", "wt");
  if(!timeFile)
    return 1;

  long DATA = 1;

  char op;

  fprintf(timeFile, "Results Digit,Toom-Cook2,FFT\n");
  fflush(timeFile);

  cerr << "Data,Results Digit,Toom-Cook2,FFT" << endl;

  while(true)
  {
    BigInteger a, b;
   
    cin >> a >> op >> b;

    string line;
    std::getline(ansFile, line);
    BigInteger ans = BigIntegerParser::Parse(line.c_str());
    
    clock_t start = clock();
    ToomCookMultiplication2 tcm2;
    vector<DataT> rToom2 = tcm2.Multiply(
      a.GetInteger(),
      b.GetInteger(),
      BigInteger::Base());
    clock_t end = clock();
    double timeTakenToomCook2 = (double)(end - start) / CLOCKS_PER_SEC;

    // cout << BigInteger(rToom2) << endl;

    int cmp = BigIntegerComparator::CompareTo(rToom2, ans.GetInteger());
    if(cmp != 0)
    {
      cerr << "Toom-Cook 2 algorithm failed." << endl;
    }

    start = clock();
    vector<DataT> fftm = FFTMultiplication::Multiply(
      a.GetInteger(),
      b.GetInteger(),
      BigInteger::Base());
    end = clock();
    double timeTakenFft = (double)(end - start) / CLOCKS_PER_SEC;

    cmp = BigIntegerComparator::CompareTo(fftm, ans.GetInteger());
    if(cmp != 0)
    {
      cerr << "FFT algorithm failed." << endl;
    }

    cerr.setf(ios::showpoint);
    cerr << DATA << "," << rToom2.size() << "," << timeTakenToomCook2 << "," << timeTakenFft << endl;
    fprintf(timeFile, "%lu,%f,%f\n", rToom2.size(), timeTakenToomCook2, timeTakenFft);
    fflush(timeFile);

    DATA++;
  
    if(cin.eof())
      break;
  }

  fclose(timeFile);

  return EXIT_SUCCESS;
}

