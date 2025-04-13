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
#include "../ops/BigIntegerIO.h"
#include "../BigIntegerParser.h"
#include "../ops/BigIntegerComparison.h"
#include "../algorithms/multiplication/ToomCookMultiplication.h"
#include "../algorithms/multiplication/ToomCookMultiplication.h"
#include "../algorithms/multiplication/ToomCookBigIntegerMultiplication.h"
#include "../algorithms/multiplication/toomcookmemoptim/ToomCookMultiplicationMemOptimized.h"
#include "../algorithms/classic/ClassicMultiplication.h"
#include "../algorithms/multiplication/FFTMultiplication.h"

#include "../algorithms/classic/ClassicDivision.h"
#include "../algorithms/divison/NewtonRaphsonDivision.h"
#include "../algorithms/divison/MontgomeryDivision.h"

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
    BigInteger ans = BigIntegerParser::Parse(line.c_str());
    clock_t end = clock();
    double timeTaken = (double)(end - start) / CLOCKS_PER_SEC;

    cerr << "Reading time: " << timeTaken << endl;

    start = clock();
    BigInteger result = ClassicDivision::DivideAndRemainder(
        a.GetInteger(),
        b.GetInteger(),
        BigInteger::Base()).first;
    end = clock();
    timeTaken = (double)(end - start) / CLOCKS_PER_SEC;

    cout << result << endl;

    cerr.setf(ios::showpoint);
    cerr << DATA << "," << result.size() << "," << timeTaken << endl;

    if (result != ans)
    {
      cerr << "Algorithm failed." << endl;
    }

    DATA++;

    if (cin.eof())
      break;
  }

  return EXIT_SUCCESS;
}
