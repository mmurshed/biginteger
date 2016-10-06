/***************************************************************************
                          main.cpp  -  description
                             -------------------
    begin                : Thu Aug 29 02:15:55 BDT 2002
    copyright            : (C) 2002 by Mahbub Murshed Suman
    email                : suman@bttb.net.bd
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <malloc.h>
#include <cmath>
#include <cstring>
#include <ctime>
#include <strstream>
#include <string>
#include <stdexcept>

using namespace std;

#include "6.7/BigInteger.h"

using namespace BigMath;

int main()
{
  cerr << BigIntegerAbout << endl;

  BigInteger m;

  while(true)
  {
    if(cin.eof()==true) break;
    cin >> m;
    cout << m << endl;
  }

  return 0;
}
