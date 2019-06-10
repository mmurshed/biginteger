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
#include <string>
#include <stdexcept>

using namespace std;

#include "../../BigInteger.h"
#include "../../BigIntegerIO.h"
#include "../../BigIntegerOperations.h"
#include "../../BigIntegerParser.h"
#include "../../algorithms/toomcook/ToomCookMultiplication.h"

#include <easy/profiler.h>

using namespace BigMath;

int main(int argc, char *argv[])
{
  if(argc < 2)
  {
    cerr << "Usage: perfanalysis [INPUT]" << endl;
    return 1;
  }

  EASY_PROFILER_ENABLE;

  ifstream ansFile(argv[1]);
  
  BigInteger a;
  BigInteger b;
  ansFile >> a >> b;

  ToomCookMultiplication tcm;
  vector<DataT> mult = tcm.Multiply(a.GetInteger(), b.GetInteger(), BigInteger::Base());

  auto blocks_count = profiler::dumpBlocksToFile("test.prof");

  std::cout << "Blocks count: " << blocks_count << std::endl;

  return 0;
}
