#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "../biginteger/BigInteger.h"
#include "../biginteger/common/Parser.h"
#include "../biginteger/ops/Comparison.h"
#include "../biginteger/algorithms/multiplication/ClassicMultiplication.h"
#include "../biginteger/algorithms/multiplication/KaratsubaMultiplication.h"
#include "../biginteger/algorithms/multiplication/ToomCookMultiplication.h"
#include "../biginteger/algorithms/multiplication/NTTMultiplication.h"

using namespace BigMath;
using namespace std;

static string RandomDigits(int digits, mt19937 &gen)
{
  uniform_int_distribution<int> first(1, 9);
  uniform_int_distribution<int> rest(0, 9);
  string value;
  value.reserve(digits);
  value.push_back((char)('0' + first(gen)));
  for (int i = 1; i < digits; ++i)
    value.push_back((char)('0' + rest(gen)));
  return value;
}

static bool CheckSize(int digits)
{
  mt19937 gen(42 + digits);
  BigInteger a = Parse(RandomDigits(digits, gen).c_str());
  BigInteger b = Parse(RandomDigits(digits, gen).c_str());
  const vector<DataT> &av = a.GetInteger();
  const vector<DataT> &bv = b.GetInteger();

  vector<DataT> classic = ClassicMultiplication::Multiply(av, bv, BigInteger::Base());
  vector<DataT> kara = KaratsubaMultiplication::Multiply(av, bv, BigInteger::Base());
  vector<DataT> toom = ToomCookMultiplication::Multiply(av, bv, BigInteger::Base());
  vector<DataT> ntt = NTTMultiplication::Multiply(av, bv, BigInteger::Base());

  bool okKara = Compare(classic, kara) == 0;
  bool okToom = Compare(classic, toom) == 0;
  bool okNtt = Compare(classic, ntt) == 0;

  cout << digits
       << " digits: kara=" << okKara
       << " toom=" << okToom
       << " ntt=" << okNtt
       << " limbs=" << av.size()
       << endl;

  return okKara && okToom && okNtt;
}

static bool CheckVectors(vector<DataT> const &a, vector<DataT> const &b, const char *label)
{
  vector<DataT> classic = ClassicMultiplication::Multiply(a, b, BigInteger::Base());
  vector<DataT> kara = KaratsubaMultiplication::Multiply(a, b, BigInteger::Base());
  vector<DataT> toom = ToomCookMultiplication::Multiply(a, b, BigInteger::Base());
  vector<DataT> ntt = NTTMultiplication::Multiply(a, b, BigInteger::Base());

  bool okKara = Compare(classic, kara) == 0;
  bool okToom = Compare(classic, toom) == 0;
  bool okNtt = Compare(classic, ntt) == 0;

  cout << label
       << ": kara=" << okKara
       << " toom=" << okToom
       << " ntt=" << okNtt
       << " limbs=" << a.size() << "x" << b.size()
       << endl;

  return okKara && okToom && okNtt;
}

int main()
{
  int sizes[] = {1, 2, 9, 10, 50, 100, 1000, 10000, 100000};
  bool ok = true;
  for (int digits : sizes)
    ok = CheckSize(digits) && ok;
  ok = CheckVectors(vector<DataT>(257, 0xFFFFFFFFULL),
                    vector<DataT>(257, 0xFFFFFFFFULL),
                    "max-carry-257") &&
       ok;
  ok = CheckVectors(vector<DataT>(513, 0xFFFFFFFFULL),
                    vector<DataT>(17, 0xFFFFFFFFULL),
                    "unbalanced-max-carry") &&
       ok;
  return ok ? 0 : 1;
}
