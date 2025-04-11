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

#include "../BigInteger.h"
#include "../BigIntegerIO.h"
#include "../BigIntegerParser.h"
#include "../BigIntegerComparator.h"
#include "../algorithms/toomcook/ToomCookMultiplication.h"
#include "../algorithms/toomcook/ToomCookMultiplication2.h"
#include "../algorithms/toomcook/ToomCookMultiplication2a.h"
#include "../algorithms/toomcookmemoptim/ToomCookMultiplicationMemOptimized.h"
#include "../algorithms/classic/ClassicMultiplication.h"
#include "../algorithms/stonehagestrassen/FFTMultiplication.h"

using namespace BigMath;

int main(int argc, char *argv[])
{
  /*
    if(argc < 3)
    {
      cerr << "Usage: toomtest [INPUT] [OUTPUT] [ANSWER]" << endl;
      return 1;
    }

    if(!freopen(argv[1],"rt",stdin))
      return 1;

    if(!freopen(argv[2],"wt",stdout))
      return 1;

    ifstream ansFile(argv[3]);
  */

  if (!freopen("input.txt", "rt", stdin))
    return 1;

  if (!freopen("output.txt", "wt", stdout))
    return 1;

  ifstream ansFile("answer.txt");

  long DATA = 1;

  char op;

  cerr << "Data,Results Digit,Time" << endl;

  while (true)
  {
    BigInteger a, b;
    cerr << "reading" << endl;

    clock_t start = clock();
    cin >> a >> op >> b;

    // vector<DataT> res1 = CommonAlgorithms::ShiftLeft(a.GetInteger(), 2);

    // cerr << a << endl;
    // cerr << res1 << endl;

    string line;
    getline(ansFile, line);
    BigInteger ans = BigIntegerParser::Parse(line.c_str());
    clock_t end = clock();
    double timeTaken = (double)(end - start) / CLOCKS_PER_SEC;

    cerr << "Reading time: " << timeTaken << endl;

    cerr << "started" << endl;

    start = clock();
    vector<DataT> result = FFTMultiplication::Multiply(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base());
    end = clock();
    timeTaken = (double)(end - start) / CLOCKS_PER_SEC;

    // cout << ans << endl;
    cout << result << endl;

    cerr.setf(ios::showpoint);
    cerr << DATA << "," << result.size() << "," << timeTaken << endl;

    int cmp = BigIntegerComparator::CompareTo(result, ans.GetInteger());
    if (cmp != 0)
    {
      cerr << "Multiplication algorithm failed." << endl;
      return cmp;
    }

    DATA++;

    if (cin.eof())
      break;
  }

  return EXIT_SUCCESS;
}
