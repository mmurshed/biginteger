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

static void Run(SizeT dividendLimbs, SizeT divisorLimbs, int k)
{
  mt19937_64 gen(98765 + dividendLimbs + divisorLimbs);
  vector<DataT> a = RandomNumber(dividendLimbs, gen);
  vector<DataT> b = RandomNumber(divisorLimbs, gen);

  double fastMs = 0, bzMs = 0;
  pair<vector<DataT>, vector<DataT>> fast, bz;

  for (int run = 0; run < k; ++run)
  {
    auto start = chrono::high_resolution_clock::now();
    fast = FastDivision::DivideAndRemainder(a, b, BigInteger::Base());
    auto end = chrono::high_resolution_clock::now();
    fastMs += chrono::duration<double, milli>(end - start).count();

    start = chrono::high_resolution_clock::now();
    bz = BurnikelZieglerDivision::DivideAndRemainder(a, b, BigInteger::Base());
    end = chrono::high_resolution_clock::now();
    bzMs += chrono::duration<double, milli>(end - start).count();
  }

  fastMs /= k;
  bzMs /= k;

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

int main(int argc, char **argv)
{
  int k = (argc > 1) ? atoi(argv[1]) : 3;
  if (k < 1) k = 1;

  cout << "Runs per size (k): " << k << endl;
  Run(1024, 512, k);
  Run(4096, 2048, k);
  Run(8192, 512, k);
  Run(16384, 512, k);
  return 0;
}
