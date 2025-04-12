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

using namespace BigMath;

#define DEBUG 0

int main(int argc, char *argv[])
{
  /*
  if(argc < 3)
  {
    cerr << "Usage: stresstest [INPUT] [OUTPUT] [ANSWER] [CSVFILE]" << endl;
    return 1;
  } */

  string input("input.txt");
  string output("output.txt");
  string answer("answer.txt");

  if (!freopen(input.c_str(), "rt", stdin))
    return 1;
  if (!freopen(output.c_str(), "wt", stdout))
    return 1;

  ifstream ansFile(answer.c_str());

  FILE *timeFile = fopen("time.csv", "wt");
  if (!timeFile)
    return 1;

  long DATA = 1;

  char op;
  bool rp;
  bool un;

  bool first = true;

  clock_t startTotal = clock();

  while (true)
  {
    BigInteger m, n;
    BigInteger r;
    if (!first)
      cout << endl;

    first = false;
    clock_t start = clock();
    cin >> m >> op >> n;
    clock_t end = clock();

    double timeTaken = (double)(end - start) / CLOCKS_PER_SEC;
    cerr << "Read in : " << timeTaken << endl;
    cerr << m.size() << " digits " << op << " " << n.size() << " digits" << endl;

    un = false;

    start = clock();

    try
    {
      switch (op)
      {
      case '+':
        r = m + n;
        break;
      case '-':
        r = m - n;
        break;
      case '*':
        r = m * n;
        break;
      case '/':
        r = m / n;
        break;
      case '%':
        r = m % n;
        break;
      case '=':
        rp = m.CompareTo(n) == 0;
        un = true;
        break;
      case '>':
        rp = m.CompareTo(n) > 0;
        un = true;
        break;
      case '<':
        rp = m.CompareTo(n) < 0;
        un = true;
        break;
      default:
        break;
      }
    }
    catch (...)
    {
    }

    end = clock();

    if (un)
    {
      if (rp)
        cout << "true";
      else
        cout << "false";
    }
    else
      cout << r;

    string line;
    getline(ansFile, line);
    BigInteger ans = BigIntegerParser::Parse(line.c_str());

    timeTaken = (double)(end - start) / CLOCKS_PER_SEC;
    cerr << "Data : " << DATA << "  Time : ";

    cerr.setf(ios::showpoint);
    cerr << timeTaken << endl;
    fprintf(timeFile, "%d,%f\n", r.size(), timeTaken);

    int cmp = BigIntegerComparator::Compare(r.GetInteger(), ans.GetInteger());
    if (cmp != 0)
    {
      cerr << "Failed." << endl;
      return cmp;
    }

    DATA++;

    if (cin.eof())
      break;
  }

  clock_t endTotal = clock();
  cerr << "Total  Time : ";
  cerr.setf(ios::showpoint);
  cerr << (double)(endTotal - startTotal) / CLOCKS_PER_SEC << endl;

  fclose(timeFile);

  return EXIT_SUCCESS;
}
