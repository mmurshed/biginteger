/**
 * BigInteger Class
 * Version 8.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef CLASSICAL_ALGORITHMS_H
#define CLASSICAL_ALGORITHMS_H

#include <vector>
#include <string>
using namespace std;

#include "BigIntegerUtil.h"

namespace BigMath
{
  // All operations are unsigned
  class ClassicalAlgorithms 
  {

    public:
    static Int Len(SizeT start, SizeT end)
    {
      return (Int) end - (Int) start + 1;
    }

    private:
    static void AddTo(
      vector<DataT> & a,
      ULong b,
      ULong base)
      {
        AddTo(a, 0, a.size() - 1, b, base);
      }

    static void AddTo(
      vector<DataT> & a, SizeT aStart, SizeT aEnd, 
      ULong b,
      ULong base)
    {
      ULong sum = b;
      ULong carry = 0;
      if(aEnd >= aStart)
      {
        sum += a[aStart];
        a[aStart] = (DataT)(sum % base);
      }
      else
      {
        a.push_back((DataT)(sum % base));
      }

      carry = sum / base;

      SizeT aPos = aStart + 1;
      while(carry > 0)
      {
        sum = carry;
        if(aEnd >= aPos)
        {
          sum += a[aPos];
          a[aPos] = (DataT)(sum % base);
        }
        else
          a.push_back((DataT)(sum % base));
        carry = sum / base;
        aPos++;
      }

      while(carry > 0)
      {
        if(aEnd >= aPos)
          a[aPos] = (DataT)(carry % base);
        else
          a.push_back((DataT)(carry % base));
        carry = carry / base;
        aPos++;
      }
    }    

    public:
    static void AddTo(
      vector<DataT>& a, SizeT aStart, SizeT aEnd, 
      vector<DataT> const& b, SizeT bStart, SizeT bEnd,
      ULong base)
    {
      Add(
          a, aStart, aEnd,
          b, bStart, bEnd,
          a, aStart,
          base);
    }

    static vector<DataT> Add(
      vector<DataT> const& a,
      vector<DataT> const& b,
      ULong base)
    {
      SizeT size = (SizeT)max(a.size(),  b.size()) + 1;
      vector<DataT> result(size);

      result[size - 1] = 0;
      
      Add(
          a, 0, (SizeT)a.size() - 1,
          b, 0, (SizeT)b.size() - 1,
          result, 0,
          base);
      
      return result;
    }

    // Implentation of addition by paper-pencil method
    // Runtime O(n), Space O(n)
    static void Add(
      vector<DataT> const& a, SizeT aStart, SizeT aEnd, 
      vector<DataT> const& b, SizeT bStart, SizeT bEnd, 
      vector<DataT>& result, SizeT rStart,
      ULong base)
    {
      // bool aZero = BigIntegerUtil::IsZero(a, aStart, aEnd);
      // bool bZero = BigIntegerUtil::IsZero(b, bStart, bEnd);
      // if(aZero && bZero)
      //   return;
      // else if(aZero)
      // {
      //   BigIntegerUtil::Copy(b, bStart, result, rStart);
      //   return;
      // }

      aEnd = min(aEnd, (SizeT) (a.size() - 1));
      bEnd = min(bEnd, (SizeT) (b.size() - 1));

      Int size = max(Len(aStart, aEnd),  Len(bStart, bEnd));
      Long carry = 0;

      for(Int i = 0; i < size; i++)
      {
        Long digitOps = 0;

        Int aPos = i + aStart;
        if(aPos <= aEnd && aPos < a.size())
          digitOps = a.at(aPos);

        digitOps += carry;
        
        Int bPos = i + bStart;
        if(bPos <= bEnd && bPos < b.size())
          digitOps += b.at(bPos);

        Int rPos = rStart + i;
        if(rPos < result.size()) 
          result.at(rPos) = (DataT)(digitOps % base);
        carry = digitOps / base;
      }
      Int rPos = rStart + size;
      if(carry > 0 && rPos < result.size())
        result.at(rPos) += carry;
    }

    static void SubtractFrom(
      vector<DataT>& a, SizeT aStart, SizeT aEnd, 
      vector<DataT> const& b, SizeT bStart, SizeT bEnd,
      ULong base)
    {
      Subtract(
           a, aStart, aEnd,
           b, bStart, bEnd,
           a, aStart,
           base);
    }

    static vector<DataT> Subtract(
      vector<DataT> const& a, 
      vector<DataT> const& b,
      ULong base)
    {
      SizeT size = (SizeT)max(a.size(),  b.size()) + 1;
      vector<DataT> result(size);

      Subtract(
           a, 0, (SizeT)a.size() - 1,
           b, 0, (SizeT)b.size() - 1,
           result, 0, base);
      
      return result;
    }

    // Implentation of subtraction by paper-pencil method
    // Assumption: a > b
    // Runtime O(n), Space O(n)
    static void Subtract(
      vector<DataT> const& a, SizeT aStart, SizeT aEnd, 
      vector<DataT> const& b, SizeT bStart, SizeT bEnd, 
      vector<DataT>& result, SizeT rStart,
      ULong base)
    {
      aEnd = min(aEnd, (SizeT) (a.size() - 1));
      bEnd = min(bEnd, (SizeT) (b.size() - 1));
      
      Int size = max(Len(aStart, aEnd),  Len(bStart, bEnd));
      Long carry = 0;

      for(Int i = 0; i < size; i++)
      {
        Long digitOps = 0;

        Int aPos = aStart + i;
        if(aPos <= aEnd && aPos < a.size())
          digitOps = a.at(aPos);

        digitOps -= carry;
        
        Int bPos = bStart + i;
        if(bPos <= bEnd && bPos < b.size())
          digitOps -= b.at(bPos);

        carry = 0;
        if(digitOps < 0)
        {
          digitOps += base;
          carry = 1;
        }

        Int rPos = rStart + i;
        if(rPos < result.size())
          result.at(rPos) = (DataT)digitOps;
      }
    }

    static void MultiplyTo(
      vector<DataT>& a,
      ULong b,
      ULong base)
    {
        MultiplyTo(
          a, 0, a.size() - 1,
          b,
          base);
    }    

    static SizeT MultiplyTo(
      vector<DataT>& a, SizeT aStart, SizeT aEnd,
      ULong b,
      ULong base)
    {
      if(b == 0 || // a times 0
        BigIntegerUtil::IsZero(a, aStart, aEnd)) // 0 times b
      {
        BigIntegerUtil::SetBit(a, aStart, aEnd, 0);
        return 0;
      }

      ULong carry = 0;

      for(Int j = aStart; j <= aEnd; j++)
      {
        ULong multiply = a.at(j);
        multiply *= b;
        multiply += carry;
        
        if(j < a.size())
        {
          a.at(j) = (DataT)(multiply % base);
        }
        else
        {
          a.push_back( (DataT)(multiply % base) );
        }
          
        carry = multiply / base;
      }

      Int j = aEnd + 1;
      while(carry > 0)
      {
        if(j < a.size())
        {
          a.at(j) = (DataT)(carry % base);
        }
        else
        {
          a.push_back( (DataT)(carry % base) );
        }
        
        carry = carry / base;
        j++;
      }

      return j;
    }    

    static vector<DataT> Multiply(
      vector<DataT> const& a,
      ULong b,
      ULong base)
      {
        vector<DataT> w(a.size());
        Multiply(
          a, 0, a.size() - 1,
          b,
          w, 0, w.size() - 1,
          base);
        return w;
      }

    static SizeT  Multiply(
      vector<DataT> const& a, SizeT aStart, SizeT aEnd,
      ULong b,
      vector<DataT>& w, SizeT wStart, SizeT wEnd,
      ULong base)
    {
      if(b == 0 || // a times 0
        BigIntegerUtil::IsZero(a, aStart, aEnd)) // 0 times b
      {
        BigIntegerUtil::SetBit(w, wStart, wEnd, 0);
        return 0;
      }

      SizeT len = Len(aStart, aEnd);
      SizeT wPos = wStart;

      ULong carry = 0;
      for(Int j = 0 ; j < len; j++)
      {
        ULong multiply = 0;
        SizeT aPos = aStart + j;
        if(aPos <= aEnd)
          multiply = a.at(aPos);
        multiply *= b;
        multiply += carry;
        
        wPos = wStart + j;
        if(wPos < w.size())
        {
          w.at(wPos) = (DataT)(multiply % base);
        }
        else
        {
          w.push_back( (DataT)(multiply % base) );
        }
          
        carry = multiply / base;
      }

      wPos = wStart + len - 1;
      while(carry > 0)
      {
        if(wPos < w.size())
        {
          w.at(wPos) = (DataT)(carry % base);
        }
        else
        {
          w.push_back( (DataT)(carry % base) );
        }
        
        carry = carry / base;
        wPos++;
      }

      return wPos;
    }

    public:
    static void MultiplyTo(
      vector<DataT> & a,
      vector<DataT> const& b,
      ULong base)
    {
      a = Multiply(a, b, base);
    }

    static vector<DataT> Multiply(
      vector<DataT> const& a,
      vector<DataT> const& b,
      ULong base,
      bool trim=false)
    {
      if(BigIntegerUtil::IsZero(a) || BigIntegerUtil::IsZero(b)) // 0 times
        return vector<DataT>();

      SizeT size = (SizeT)(a.size() + b.size() + 1);
      vector<DataT> result(size);

      Multiply(a, 0, (SizeT)a.size() - 1, b, 0, (SizeT)b.size() - 1, result, 0, base);

      if(trim)
      {
        BigIntegerUtil::TrimZeros(result);
      }

      return result;      
    }

    // Implentation of multiplication by paper-pencil method
    // Classical algorithm
    // Runtime O(n^2), Space O(n)
    static SizeT Multiply(
      vector<DataT> const& a, SizeT aStart, SizeT aEnd, 
      vector<DataT> const& b, SizeT bStart, SizeT bEnd, 
      vector<DataT>& result, SizeT rStart,
      ULong base)
    {
      aEnd = min(aEnd, (SizeT) (a.size() - 1));
      bEnd = min(bEnd, (SizeT) (b.size() - 1));

      SizeT k = rStart;
      SizeT lenA = aEnd - aStart + 1;
      
      for(SizeT i = bStart; i <= bEnd; i++)
      {
        ULong carry = 0;
        SizeT jStart = rStart + (i - bStart);
        for(SizeT j = aStart; j <= aEnd; j++)
        {
          SizeT k = jStart + (j - aStart);
          ULong multiply = a.at(j);
          multiply *= b.at(i);
          multiply += result.at(k);
          multiply += carry;
          
          result.at(k) = (DataT)(multiply % base);
          carry = multiply / base;
        }
        k = jStart + lenA;
        result.at(k) = (DataT)carry;
      }
      return k;
    }

    static void DivideTo(
      vector<DataT>& u,
      DataT d,
      ULong base)
      {
        Divide(u, 0, u.size() - 1, d, u, 0, u.size() - 1, base);
      }

    static void DivideTo(
      vector<DataT>& u, SizeT uStart, SizeT uEnd,
      DataT d,
      ULong base)
      {
        Divide(u, uStart, uEnd, d, u, uStart, uEnd, base);
      }

    static vector<DataT> Divide(
      vector<DataT> const& u,
      DataT d,
      ULong base)
      {
        SizeT n = (SizeT)u.size();
        // Divide (u_n−1 . . . u_1 u_0)_b by d.
        vector<DataT> w(n);

        Divide(u, 0, u.size() - 1, d, w, 0, w.size() - 1, base);

        return w;
      }
    
    static void Divide(
      vector<DataT> const& u, SizeT uStart, SizeT uEnd,
      DataT d,
      vector<DataT>& w, SizeT wStart, SizeT wEnd,
      ULong base)
      {
        SizeT n = Len(uStart, uEnd);
        // Divide (u_n−1 . . . u_1 u_0)_b by d.

        // Set r = 0, j = n − 1.
        ULong r = 0;
        for(Int j = n - 1; j >= 0; j--)
        {
          // Set val = r * b + u_j
          ULong val = r * base;
          SizeT uPos = uStart + j;
          if(uPos <= uEnd)
            val += u.at(uPos);
          // Set w_j = ⌊val / d⌋
          SizeT wPos = wStart + j;
          if(wPos <= wEnd)
            w.at(wPos) = (DataT)(val / d);
          else
          {
            w.push_back((DataT)(val / d));
          }
          // Set r = val mod d
          r = val % d;
        }
      }

    static void DivideTo(
      vector<DataT> & a,
      vector<DataT> const& b,
      ULong base)
    {
      vector< vector<DataT> > results = DivideAndRemainder(a, b, base);
      a = results[0];
    }

    // Implentation of division by paper-pencil method
    // Assumption: a > b
    // See: D.E.Knuth 4.3.1
    // Runtime O(n^2), Space O(n)
    static vector< vector<DataT> > DivideAndRemainder(
      vector<DataT> const& a,
      vector<DataT> const& b,
      ULong base)
    {
      // Given nonnegative integers u = (u_m+n−1 . . . u_1 u_0)b and v = (v_n−1 . . . v_1 v_0)_b, 
      // where v_n−1 != 0 and n > 1, we form the radix-b quotient ⌊u/v⌋ = (q_m q_m–1 . . . q_0)_b 
      // and the remainder u mod v = (r_n−1 . . . r_1 r_0)b.
      
      // u is indexed from 0 to m + n - 1
      // m + n = u.size()
      // n = v.size()
      // m = u.size() - v.size()
      SizeT n = (SizeT)b.size();
      SizeT m = (SizeT)(a.size() - b.size());

      // 1. Normalize
      // Set d = ⌊b / (v_n−1 + 1)⌋.
      DataT d = (DataT)(base / (b[n - 1] + 1));

      // Set (u_m+n u_m+n−1 . . . u_1 u_0)_b equal to (u_m+n−1 . . . u_1 u_0)_b times d
      // Notice the introduction of a new digit position u_m+n at the left of u_m+n−1
      // If d = 1, all we need to do in this step is to set u_m+n = 0.
      // On a binary computer it may be preferable to choose d to be a power of 2 
      // instead of using the value suggested here.
      // Any value of d that results in v_n−1 >= ⌊b/2⌋ will suffice.
      vector<DataT> u = Multiply(a, d, base);
      if(u.size() <= m + n)
        u.push_back(0);

      // Set (v_n−1 . . . v_1 v_0)_b equal to (v_n−1 . . . v_1 v_0)_b times d.
      vector<DataT> v = Multiply(b, d, base);

      vector<DataT> q (m + 1);

      // 2. Loop over j from m to 0
      // Set j = m. (The loop on j, steps 2 through 7, will be essentially a division of 
      // (u_j+n . . . u_j+1 u_j)_b by (v_n−1 . . . v_1 v_0)_b to get a single quotient digit q_j
      for(Int j = m; j >= 0; j--)
      {
        // 3. Calculate q_hat
        // Set val = u_j+n * b + u_j+n−1
        ULong val = u[j + n];
        val *= base;
        val += u[j + n - 1];
        // Set q_hat = ⌊val / v_n−1⌋
        ULong qhat = val / v[n - 1];
        // Set r_hat be the remainder, r_hat = val mod v_n−1.
        ULong rhat = val % v[n - 1];

        // Now test if q_hat >= b or q_hat * v_n−2 > b * r_hat + u_j+n−2;
        // The test on v_n−2 determines at high speed most of the cases in which 
        // the trial value q_hat is one too large, and it eliminates all cases where q_hat is two too large.
        while(qhat >= base || (n > 1 && qhat * v[n - 2] > rhat * base + u[j + n - 2]) )
        {
          // Decrease q_hat by 1
          qhat--;
          // Increase r_hat by v_n−1
          rhat += v[n - 1];
          // Repeat this test if r_hat < b.
          if(rhat >= base)
            break;
        }

        // 4. Multiply and subtract
        // Replace (u_j+n u_j+n−1 . . . u_j)_b by (u_j+n u_j+n−1 . . . u_j)_b - q_hat * (0 v_n-1 . . . v_1 v_0)_b
        // This computation consists of a simple multiplication by a one-place number, combined with a subtraction.
        ULong carry = 0;
        for(SizeT i = 0; i <= n; i++)
        {
          ULong mult = carry;
          if(i < v.size())
            mult += qhat * v[i];
          carry = mult / base;
          mult %= base;

          // The digits (u_j+n, u_j+n−1, . . . , u_j) should be kept positive
          // If the result of this step is actually negative, (u_j+n u_j+n−1 . . . u_j)_b should be left as the true 
          // value plus b^n+1, namely as the b’s complement of the true value, and a "borrow" to the left should be remembered.
          
          Long diff = (Long)u[i + j] - (Long)mult;
          if(diff < 0)
          {
            diff += base;
            carry++;
          }
          u[i+j]  = (DataT)diff;
        }

        // 5. Test remainder
        // Set q_j = q_hat
        q[j] = (DataT)qhat;
        
        // Result was negative
        // The probability that this step is necessary is very small, on the order of only 2/b.
        // Test data to activate this step should therefore be specifically contrived when debugging.
        if(carry > 0)
        {
          // 6. Add back
          // Decrease q_j by 1
          q[j]--;
          // add (0 v_n−1 . . . v_1 v_0)_b to (u_j+n u_j+n−1 . . . u_j+1 u_j)_b.
          for(SizeT i = 0; i <= n; i++)
          {
            ULong sum = 0;
            if(i < v.size())
              sum += v[i];

            sum += u[i + j];
            u[i + j] = (DataT)(sum % base);
            // A carry will occur to the left of u_j+n, and it should be ignored since it cancels with 
            // the borrow that occurred in step 4.
          }
        }
      } // 7. Loop on j

      // 8. Unnormalize
      // Now (q_m . . . q_1 q_0)b is the desired quotient, and 
      // the desired remainder may be obtained by dividing (u_n−1 . . . u_1 u_0)_b by d.
      vector<DataT> w(n);

      // Set r = 0, j = n − 1.
      ULong r = 0;
      for(Int j = n - 1; j >= 0; j--)
      {
        // Set val = r * b + u_j
        ULong val = r * base + u[j];
        // Set w_j = ⌊val / d⌋
        w[j] = (DataT)(val / d);
        // Set r = val mod d
        r = val % d;
      }

      vector< vector<DataT> > result(2);
      result[0] = q;
      result[1] = w;

      return result;
    }

    static vector<DataT> ConvertBase(
      vector<DataT> const& bigIntB1,
      ULong base1,
      ULong base2)
      {
        return ConvertBase(
          bigIntB1, 0, bigIntB1.size() - 1,
          base1,
          base2);
      }

    // Convert an integer from base1 to base2
    // Donald E. Knuth 4.4
    static vector<DataT> ConvertBase(
      vector<DataT> const& bigIntB1, SizeT start, SizeT end,
      ULong base1,
      ULong base2)
    {
      if(base1 == base2)
      {
        vector<DataT> bigIntB2(bigIntB1);
        return bigIntB2;
      }

      vector<DataT> bigIntB2(1, 0);

      for(Int i = (Int)end; i >= (Int)start; i--)
      {
        MultiplyTo(
          bigIntB2,
          base1,
          base2);
        AddTo(
          bigIntB2,
          bigIntB1[i],
          base2);
      }

      return bigIntB2;
    }

    // Returns an integer by shifting n places
    // Equivalent to a * B^n
    static vector<DataT> ShiftLeft(
      vector<DataT> const& bigInt,
      SizeT shift)
    {
      SizeT size = (SizeT)bigInt.size() + shift;
      vector<DataT> result(size);

      Int j = size - 1;

      // Copy
      for(Int i = (Int)bigInt.size() - 1; i >= 0; i--, j--)
      {
        result[j] = bigInt[i];
      }
      
      // Insert zeros
      for(;j >= 0; j--)
      {
        result[j] = 0;
      }

      return result;
    }
   };
}

#endif
