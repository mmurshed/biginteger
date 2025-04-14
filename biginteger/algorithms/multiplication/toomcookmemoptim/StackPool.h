/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef STACK_POOL
#define STACK_POOL

#include <vector>
#include <stack>
using namespace std;

#include "../../../common/Util.h"

namespace BigMath
{
  typedef pair<Int, Int> Range;

  template <class T>
  class StackPool
  {
  public:
    vector<T> data;
    vector<Range> ranges;

    StackPool(SizeT n = 0) : data(n, 0)
    {
    }

    void Resize(SizeT n)
    {
      data.resize(n);
    }

    void Push(Int first, Int second)
    {
      ranges.push_back({first, second});
    }

    void Push(vector<T> const &u, SizeT p)
    {
      SizeT start = ranges.empty() ? 0 : ranges.back().second + 1;
      SizeT end = start + p - 1;

      Copy(
          u, 0, u.size() - 1,
          data, start, end);

      ranges.push_back({start, end});
    }

    Range Pop()
    {
      Range range = ranges.back();
      ranges.pop_back();
      return range;
    }

    Range Top()
    {
      return ranges.back();
    }
  };
}

#endif
