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

#include "../biginteger/BigInteger.h"
#include "../biginteger/ops/IO.h"
#include "../biginteger/common/Parser.h"
#include "../biginteger/ops/Comparison.h"
#include "../biginteger/algorithms/multiplication/ToomCookMultiplication.h"
#include "../biginteger/algorithms/multiplication/ToomCookMultiplication.h"
#include "../biginteger/algorithms/multiplication/ToomCookBigIntegerMultiplication.h"
#include "../biginteger/algorithms/multiplication/toomcookmemoptim/ToomCookMultiplicationMemOptimized.h"
#include "../biginteger/algorithms/multiplication/ClassicMultiplication.h"
#include "../biginteger/algorithms/multiplication/FFTMultiplication.h"

#include "../biginteger/algorithms/division/KnuthDivision.h"
#include "../biginteger/algorithms/division/NewtonRaphsonDivision.h"
#include "../biginteger/algorithms/division/MontgomeryDivision.h"
#include "../biginteger/algorithms/division/FFTDivision.h"

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

    clock_t start = clock();
    cin >> a >> op >> b;

    string line;
    getline(ansFile, line);
    BigInteger ans = Parse(line.c_str());
    clock_t end = clock();
    double timeTaken = (double)(end - start) / CLOCKS_PER_SEC;

    cerr << "Reading time: " << timeTaken << endl;

    start = clock();
    auto result = FFTDivision::DivideAndRemainder(a.GetInteger(),b.GetInteger(),BigInteger::Base()).first;
    end = clock();
    timeTaken = (double)(end - start) / CLOCKS_PER_SEC;

    cout << result << endl;

    cerr.setf(ios::showpoint);
    cerr << DATA << "," << result.size() << "," << timeTaken << endl;

    if (Compare(result, ans.GetInteger()) != 0)
    {
      cerr << "Algorithm failed." << endl;
      cerr << "a: " << a << endl;
      cerr << "b: " << b << endl;
      cerr << "result: " << result << endl;
      cerr << "correct: " << ans << endl;
    }

    DATA++;

    if (cin.eof())
      break;
  }

  return EXIT_SUCCESS;
}
