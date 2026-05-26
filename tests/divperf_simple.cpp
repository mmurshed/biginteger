#include <chrono>
#include <iostream>
#include <random>
#include <vector>

#include "biginteger/BigInteger.h"
#include "biginteger/common/Comparator.h"
#include "biginteger/algorithms/division/BurnikelZieglerDivision.h"
#include "biginteger/algorithms/division/FastDivision.h"

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

static void Run(SizeT dividendLimbs, SizeT divisorLimbs)
{
  mt19937_64 gen(98765 + dividendLimbs + divisorLimbs);
  vector<DataT> a = RandomNumber(dividendLimbs, gen);
  vector<DataT> b = RandomNumber(divisorLimbs, gen);

  auto start = chrono::high_resolution_clock::now();
  auto fast = FastDivision::DivideAndRemainder(a, b, BigInteger::Base());
  auto end = chrono::high_resolution_clock::now();
  double fastMs = chrono::duration<double, milli>(end - start).count();

  start = chrono::high_resolution_clock::now();
  auto bz = BurnikelZieglerDivision::DivideAndRemainder(a, b, BigInteger::Base());
  end = chrono::high_resolution_clock::now();
  double bzMs = chrono::duration<double, milli>(end - start).count();

  bool match = Compare(fast.first, bz.first) == 0 && Compare(fast.second, bz.second) == 0;

  cout << dividendLimbs << "x" << divisorLimbs
       << " fast=" << fastMs << "ms"
       << " bz=" << bzMs << "ms"
       << " speedup=" << (bzMs == 0 ? 0 : fastMs / bzMs)
       << " match=" << match
       << endl;

  if (!match)
    exit(1);
}

int main()
{
  Run(1024, 512);
  Run(4096, 2048);
  Run(8192, 512);
  Run(16384, 512);
  return 0;
}
