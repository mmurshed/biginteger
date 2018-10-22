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

#include "BigInteger.h"
#include "BigIntegerIO.h"
#include "BigIntegerOperations.h"
#include "BigIntegerParser.h"
#include "ToomCookAlgorithm.h"

using namespace BigMath;

int main(int argc, char *argv[])
{
  // BigInteger zero1;
  // bool z1 = zero1.IsZero(); // should be true

  // BigInteger zero2(5, false, 0);
  // bool z2 = zero2.IsZero(); // should be true

  // vector<DataT> a(zero2.GetInteger());
  // a[0] = 200;
  // BigInteger zero3(a, false);
  // bool z3 = zero3.IsZero(); // should be false

  // BigInteger zero1;
  // SizeT z1 = zero1.TrimZeros(); // should be 0

  // BigInteger zero2(5, false, 0);
  // vector<DataT> a(zero2.GetInteger());

  // SizeT z2 = zero2.TrimZeros(); // should be 5

  // a[0] = 200;
  // BigInteger zero3(a, false);
  // SizeT z3 = zero3.TrimZeros(); // should be 4

  // BigInteger zero1;
  // SizeT z1 = BigIntegerUtil::FindNonZeroByte(zero1.GetInteger()); // should be 0

  // BigInteger zero2(5, false, 0);
  // vector<DataT> a(zero2.GetInteger());

  // SizeT z2 = BigIntegerUtil::FindNonZeroByte(zero1.GetInteger()); // should be 0

  // a[0] = 200;
  // SizeT z3 = BigIntegerUtil::FindNonZeroByte(a); // should be 1

  BigInteger a = BigIntegerParser::Parse("82214287154264689604931794456724899414830547454821219178655835782565433365697983603518200983943510370081527699212628774843110186377172181980346514479605255547776236222412497742875656466706870594441529030048053593491548291546792328985506619217689815757748019133330759804793997550263816520811587878801354548277612602635107306295655161497370274892234924238408428179018411074675399864220904101233001017646549457125057137713246563887796510030876864103826492375814664295437708163761682514274698721074952682557079550789547397144242202986653944057693837758014199337539657543163113884161613550718887543070980228");
  BigInteger b = BigIntegerParser::Parse("18397345861553567652641514702113650630970754164581373527916193718645396260872396528569008683208407234874497922920031155051680425644240787486966840741976049543221725184622103569091906271138318001");

  // BigInteger& a = BigIntegerParser::Parse("14795936");
  // BigInteger& b = BigIntegerParser::Parse("33202659");

  // BigInteger& add = a + b;
  // string& str = BigIntegerParser::ToString(add);
  // bool cmp = (str == "4898444632449261293846336479321434");

  // BigInteger& sub = a - b;
  // str = sub.ToString();
  // cmp = (str == "123456789012345678901230000000");

  ToomCookAlgorithm tca;
  vector<DataT> mult = tca.MultiplyUnsigned(a.GetInteger(), b.GetInteger(), BigInteger::Base());
  string str = BigIntegerParser::ToString(mult);
  bool cmp = (str == "1512524675538088125755323413880095718933553181068000604535124927510223813860588304129125213254416954594759562548351581521617768525991990045478806879405992703532936635045982578128592563583651657606727166805858688694596312518668616892382307130543858966294346572318740504369765952779006011857593352964825981702999950051197593901099820461620900912177499417271653342109848514933946580997230216959370320602825654244224088706638006739885354067490608646841977189440249693195824836814177317188359572768945024088505009538891649761240113722703383770340279834193757371409808753708275666615086697618333260967175739065787086237370623068895449420368028848084973313932929803578710146152148001147857318031638690139077451413568853854704881347246344779305039423272276083001652268549982242789188050971589174247484228");
  // vector<DataT> re = tca.MultiplyRPart(b.GetInteger(), 4, 0, 6, BigInteger::Base());
  // string str = BigIntegerParser::ToString(re);
  // bool cmp = (str == "0");

  // BigInteger& div = a / b;
  // string& str = div.ToString();
  // bool cmp = (str == "44557413417469525");

  // BigInteger& rem = a % b;
  // string& str = rem.ToString();
  // bool cmp = (str == "1494561306735");

  return 0;
}
