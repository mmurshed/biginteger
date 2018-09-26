/**
 * BigInteger Class
 * Version 8.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef CLASSICAL_ALGORITHMS_H
#define CLASSICAL_ALGORITHMS_H

#include <vector>
using namespace std;

#include "BigIntegerUtil.h"
#include "BigInteger.h"

namespace BigMath
{
  class ClassicalAlgorithms 
  {
    // Unsigned comparison
    public:
    // Compares this with `with' irrespective of sign
    // Returns
    // 0 if equal
    // +value if this > with
    // -value if this < with
    static Int UnsignedCompareTo(vector<DataT> const& a, vector<DataT> const& b)
    {
      // Case with zero
      bool aIsZero = BigIntegerUtil::IsZero(a);
      bool bIsZero = BigIntegerUtil::IsZero(b);

      if(aIsZero && bIsZero)
        return 0;
      else if(!aIsZero && bIsZero)
        return 1;
      else if(aIsZero && !bIsZero)
        return -1;

      // Different in size
      Long diff = a.size();
      diff -= b.size();
      if(diff != 0)
        return diff;

      // Both ints have same number of digits
      Int cmp = 0;
      for(SizeT i = a.size() - 1; i >= 0; i--)
      {
        diff = a[i];
        diff -= b[i];
        if(diff != 0)
          return diff;
      }

      return cmp;
    }

    public:
    static void AddToUnsigned(vector<DataT> & a, ULong const& b, ULong const& base)
    {
      ULong sum = b;
      ULong carry = 0;
      if(a.size() > 0)
      {
        sum += a[0];
        a[0] = sum % base;
      }
      else
      {
        a.push_back(sum % base);
      }

      carry = sum / base;

      SizeT i = 1;
      while(carry > 0 && i < a.size())
      {
        sum = a[i] + carry;
        a[i] = sum % base;
        carry = sum / base;
        i++;
      }

      while(carry > 0)
      {
        a.push_back(carry % base);
        carry = carry / base;
      }
    }

    static vector<DataT>& AddUnsigned(vector<DataT> const& a, vector<DataT> const& b, ULong const& base)
    {
      SizeT size = max(a.size(),  b.size()) + 1;
      vector<DataT>& result = *new vector<DataT>(size);

      Long carry = 0;

      for(SizeT i = 0; i < size; i++)
      {
        Long digitOps = 0;

        if(i < a.size())
          digitOps = a[i];

        digitOps += carry;
                
        if(i < b.size())
          digitOps += b[i];

        result[i] = digitOps % base;
        carry = digitOps / base;
      }
      result[size-1] = carry;
      
      BigIntegerUtil::TrimZeros(result);

      return result;
    }

    // Implentation of subtraction by paper-pencil method
    // Assumption: a > b
    static vector<DataT>& SubtractUnsigned(vector<DataT> const& a, vector<DataT> const& b, ULong const& base)
    {
      SizeT size = max(a.size(),  b.size()) + 1;
      vector<DataT>& result = *new vector<DataT>(size);

      Long carry = 0;

      for(SizeT i = 0; i < size; i++)
      {
        Long digitOps = 0;

        if(i < a.size())
          digitOps = a[i];

        digitOps -= carry;
        
        if(i < b.size())
          digitOps -= b[i];

        carry = 0;
        if(digitOps < 0)
        {
          digitOps += base;
          carry = 1;
        }

        result[i] = digitOps;
      }
      
      BigIntegerUtil::TrimZeros(result);

      return result;
    }
    
    static void MultiplyToUnsigned(vector<DataT>& a, ULong const& b, ULong const& base)
    {
      if(BigIntegerUtil::IsZero(a)) // 0 times b
        return;
      if(b == 0) // a times 0
      {
        a.clear();
        return;
      }

      ULong carry = 0;
      for(Int j = 0 ; j < a.size(); j++)
      {
        ULong multiply = a[j];
        multiply *= b;
        multiply += carry;
        
        a[j] = multiply % base;
        carry = multiply / base;
      }

      while(carry > 0)
      {
        a.push_back(carry % base);
        carry = carry / base;
      }
    }

    // Implentation of multiplication by paper-pencil method
    // Classical algorithm
    static vector<DataT>& MultiplyUnsigned(vector<DataT> const& a, vector<DataT> const& b, ULong const& base)
    {
      SizeT size = a.size() + b.size() + 1;
      vector<DataT>& result = *new vector<DataT>(size);
      
      for(SizeT i = 0; i < b.size(); i++)
      {
        ULong carry = 0;
        for(SizeT j = 0 ; j < a.size(); j++)
        {
          ULong multiply = a[j];
          multiply *= b[i];
          multiply += result[i + j];
          multiply += carry;
          
          result[i + j] = multiply % base;
          carry = multiply / base;
        }
        result[i + a.size()] = carry;
      }

      BigIntegerUtil::TrimZeros(result);

      return result;
    }

    // Implentation of division by paper-pencil method
    // Assumption: a > b
    // See: D.E.Knuth 4.3.1
    static vector< vector<DataT> >& DivideAndRemainderUnsigned(vector<DataT> const& a, vector<DataT> const& b, ULong const& base)
    {
      // Given nonnegative integers u = (u_m+n−1 . . . u_1 u_0)b and v = (v_n−1 . . . v_1 v_0)_b, 
      // where v_n−1 != 0 and n > 1, we form the radix-b quotient ⌊u/v⌋ = (q_m q_m–1 . . . q_0)_b 
      // and the remainder u mod v = (r_n−1 . . . r_1 r_0)b.
      
      // u is indexed from 0 to m + n - 1
      // m + n = u.size()
      // n = v.size()
      // m = u.size() - v.size()
      SizeT n = b.size();
      SizeT m = a.size() - b.size();

      // 1. Normalize
      // Set d = ⌊b / (v_n−1 + 1)⌋.
      vector<DataT> d(1);
      d[0] = base / (b[n-1] + 1);

      // Set (u_m+n u_m+n−1 . . . u_1 u_0)_b equal to (u_m+n−1 . . . u_1 u_0)_b times d
      // Notice the introduction of a new digit position u_m+n at the left of u_m+n−1
      // If d = 1, all we need to do in this step is to set u_m+n = 0.
      // On a binary computer it may be preferable to choose d to be a power of 2 
      // instead of using the value suggested here.
      // Any value of d that results in v_n−1 >= ⌊b/2⌋ will suffice.
      vector<DataT>& u = MultiplyUnsigned(a, d, base);
      if(u.size() <= m+n)
        u.push_back(0);

      // Set (v_n−1 . . . v_1 v_0)_b equal to (v_n−1 . . . v_1 v_0)_b times d.
      vector<DataT>& v = MultiplyUnsigned(b, d, base);

      vector<DataT>& q = *new vector<DataT>(m + 1);

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
          u[i+j]  = diff;
        }

        // 5. Test remainder
        // Set q_j = q_hat
        q[j] = qhat;
        
        // Result was negative
        // The probability that this step is necessary is very small, on the order of only 2/b.
        // Test data to activate this step should therefore be specifically contrived when debugging.
        if(carry > 0)
        {
          // 6. Add back
          // Decrease q_j by 1
          q[j]--;
          // add (0 v_n−1 . . . v_1 v_0)_b to (u_j+n u_j+n−1 . . . u_j+1 u_j)_b.
          carry = 0;
          for(SizeT i = 0; i <= n; i++)
          {
            ULong sum = 0;
            if(i < v.size())
              sum += v[i];

            sum += u[i + j];
            u[i + j] = sum % base;
            carry = sum / base;
            // A carry will occur to the left of u_j+n, and it should be ignored since it cancels with 
            // the borrow that occurred in step 4.
          }
        }
      } // 7. Loop on j

      // 8. Unnormalize
      // Now (q_m . . . q_1 q_0)b is the desired quotient, and 
      // the desired remainder may be obtained by dividing (u_n−1 . . . u_1 u_0)_b by d.
      vector<DataT>& w = *new vector<DataT>(n);

      // Set r = 0, j = n − 1.
      ULong r = 0;
      for(Int j = n - 1; j >= 0; j--)
      {
        // Set val = r * b + u_j
        ULong val = r * base + u[j];
        // Set w_j = ⌊val / d⌋
        w[j] = val / d[0];
        // Set r = val mod d
        r = val % d[0];
      }

      BigIntegerUtil::TrimZeros(q);
      BigIntegerUtil::TrimZeros(w);
      vector< vector<DataT> >& result = *new vector< vector<DataT> >(2);
      result[0] = q;
      result[1] = w;

      return result;
    }

    // Convert an integer from base1 to base2
    // Donald E. Knuth 4.4
    static vector<DataT>& ConvertBase(vector<DataT> const& bigIntB1, ULong base1, ULong base2)
    {
      if(base1 == base2)
      {
        return *new vector<DataT>(bigIntB1);
      }

      vector<DataT>& bigIntB2 = *new vector<DataT>();

      for(Int i = bigIntB1.size() - 1 ; i >= 0; i--)
      {
        MultiplyToUnsigned(bigIntB2, base1, base2);
        AddToUnsigned(bigIntB2, bigIntB1[i], base2);
      }

      return bigIntB2;
    }

    // Returns an integer by shifting n places
    // Equivalent to a * B^n
    static vector<DataT>& ShiftLeft(vector<DataT> const& bigInt, SizeT shift)
    {
      SizeT size = bigInt.size() + shift;
      vector<DataT>& result = *new vector<DataT>(size);

      Int j = size - 1;

      // Copy
      for(Int i = bigInt.size() - 1; i >= 0; i--, j--)
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