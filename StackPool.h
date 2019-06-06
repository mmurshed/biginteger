/**
 * BigInteger Class
 * Version 8.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef STACK_POOL_H
#define STACK_POOL_H

#include <vector>
#include <stack>
using namespace std;

namespace BigMath
{
  typedef pair<Int, Int> Range;

  template<class T>
  class StackPool
  {
    public:
    vector<T> data;
    stack<Range> ranges;
    
    StackPool(SizeT n = 0) : data(n, 0)
    {}

    void Resize(SizeT n)
    {
      data.resize(n);
    }

    void Push(Int first, Int second)
    {
      ranges.push(make_pair(first, second));
    }

    void Push(vector<T> const& u, SizeT p)
    {
      SizeT start = ranges.empty() ? 0 : ranges.top().second + 1;
      SizeT end = start + p - 1;
      BigIntegerUtil::Copy(u, 0, u.size() - 1, data, start, end);
      ranges.push(make_pair(start, end));
    }

    Range Pop()
    {
      Range range = ranges.top();
      ranges.pop();
      return range;
    }
    Range Top()
    {
      return ranges.top();
    }
  };
}

#endif
