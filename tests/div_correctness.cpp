#include <iostream>
#include <random>
#include <vector>

#include "../biginteger/BigInteger.h"
#include "../biginteger/common/Comparator.h"
#include "../biginteger/algorithms/Addition.h"
#include "../biginteger/algorithms/Division.h"
#include "../biginteger/algorithms/Multiplication.h"
#include "../biginteger/algorithms/division/BurnikelZieglerDivision.h"
#include "../biginteger/algorithms/division/ClassicDivision.h"
#include "../biginteger/algorithms/division/FastDivision.h"
#include "../biginteger/algorithms/division/NewtonDivision.h"
#include "../biginteger/algorithms/division/ReciprocalDivision.h"

using namespace BigMath;
using namespace std;

static vector<DataT> RandomNumber(SizeT limbs, mt19937_64 &gen)
{
  uniform_int_distribution<uint64_t> digit(0, 0xFFFFFFFFULL);
  vector<DataT> value(limbs);
  for (SizeT i = 0; i < limbs; ++i)
    value[i] = digit(gen);
  if (!value.empty() && value.back() == 0)
    value.back() = 1 + (digit(gen) & 0xFFFF);
  TrimZeros(value);
  return value.empty() ? vector<DataT>{0} : value;
}

static bool CheckDivision(vector<DataT> const &a, vector<DataT> const &b, const char *label)
{
  auto qr = DivideAndRemainder(a, b, BigInteger::Base());
  auto bz = BurnikelZieglerDivision::DivideAndRemainder(a, b, BigInteger::Base());
  vector<DataT> recomposed = Add(Multiply(qr.first, b, BigInteger::Base()), qr.second, BigInteger::Base());
  bool identity = Compare(recomposed, a) == 0;
  bool remainderRange = Compare(qr.second, b) < 0;
  bool bzMatch = Compare(qr.first, bz.first) == 0 && Compare(qr.second, bz.second) == 0;

  bool newtonMatch = true;
  if (b.size() > 1)
  {
    auto nt = NewtonDivision::DivideAndRemainder(a, b, BigInteger::Base());
    newtonMatch = Compare(qr.first, nt.first) == 0 && Compare(qr.second, nt.second) == 0;
  }

  cout << label
       << ": identity=" << identity
       << " remainder=" << remainderRange
       << " bz=" << bzMatch
       << " newton=" << newtonMatch
       << " limbs=" << a.size() << "x" << b.size()
       << endl;

  return identity && remainderRange && bzMatch && newtonMatch;
}

static bool CheckScalar(vector<DataT> const &a, DataT b, const char *label)
{
  auto qr = DivideAndRemainder(a, b, BigInteger::Base());
  vector<DataT> recomposed = Add(ClassicMultiplication::Multiply(qr.first, b, BigInteger::Base()),
                                 qr.second,
                                 BigInteger::Base());
  bool identity = Compare(recomposed, a) == 0;
  bool remainderRange = qr.second.size() == 1 && qr.second[0] < b;

  cout << label
       << ": identity=" << identity
       << " remainder=" << remainderRange
       << " limbs=" << a.size()
       << endl;

  return identity && remainderRange;
}

static bool CheckRepeatedDivisor(SizeT divisorLimbs, mt19937_64 &gen)
{
  vector<DataT> b = RandomNumber(divisorLimbs, gen);
  ReciprocalDivision::Divider divider(b, BigInteger::Base());

  bool ok = true;
  SizeT dividendSizes[] = {divisorLimbs, divisorLimbs + 1, divisorLimbs * 2, divisorLimbs * 8};
  for (SizeT dividendLimbs : dividendSizes)
  {
    vector<DataT> a = RandomNumber(dividendLimbs, gen);
    auto expected = DivideAndRemainder(a, b, BigInteger::Base());
    auto actual = divider.DivideAndRemainder(a);
    bool match = Compare(expected.first, actual.first) == 0 && Compare(expected.second, actual.second) == 0;
    cout << "reciprocal-reuse"
         << ": match=" << match
         << " limbs=" << a.size() << "x" << b.size()
         << endl;
    ok = match && ok;
  }

  return ok;
}

int main()
{
  mt19937_64 gen(12345);
  bool ok = true;

  ok = CheckDivision(vector<DataT>{5}, vector<DataT>{9}, "less-than") && ok;
  ok = CheckDivision(vector<DataT>{0xFFFFFFFFULL, 0xFFFFFFFFULL}, vector<DataT>{0xFFFFFFFFULL}, "scalar-max") && ok;
  ok = CheckDivision(vector<DataT>(64, 0xFFFFFFFFULL), vector<DataT>(32, 0xFFFFFFFFULL), "max-64x32") && ok;

  SizeT sizes[][2] = {
      {2, 1}, {4, 2}, {16, 8}, {64, 31}, {256, 97}, {1024, 513}, {4096, 2048}, {8192, 512}};

  for (auto &size : sizes)
  {
    vector<DataT> a = RandomNumber(size[0], gen);
    vector<DataT> b = RandomNumber(size[1], gen);
    ok = CheckDivision(a, b, "random") && ok;
  }

  ok = CheckScalar(RandomNumber(512, gen), 3, "scalar-3") && ok;
  ok = CheckScalar(RandomNumber(512, gen), 0xFFFFFFFFULL, "scalar-max-digit") && ok;
  ok = CheckRepeatedDivisor(128, gen) && ok;
  ok = CheckRepeatedDivisor(1024, gen) && ok;

  return ok ? 0 : 1;
}
