/**
 * BigInteger Class
 * Version 10.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef TOOMCOOK_MULTIPLICATION_2A
#define TOOMCOOK_MULTIPLICATION_2A

#include <vector>
using namespace std;

#include "../../BigInteger.h"
#include "../BigIntegerUtil.h"
#include "../classic/ClassicMultiplication.h"
#include "../../ops/BigIntegerAddition.h"
#include "../../ops/BigIntegerSubtraction.h"

namespace BigMath
{
  class ToomCookMultiplication2a
  {
    private:

    // Threshold for switching to the naive multiplication.
    const SizeT baseCaseThreshold = 256;      

    // For evaluation, we need to compute 5 points (e.g., points for Toom–3: 0, 1, -1, 2, ∞).
    const int numPoints = 5;

      
    public:
    vector<DataT> Multiply(
      vector<DataT> const& a,
      vector<DataT> const& b,
      BaseT base)
    {
      SizeT n = (SizeT)max(a.size(), b.size());

      if (n < baseCaseThreshold) {
          return ClassicMultiplication::Multiply(a, b, base);
      }

      // For a Toom–3 algorithm, we typically split the input into 3 parts.
      // Calculate the part size (rounding up if necessary).
      SizeT partSize = (n + 2) / 3;
            
      // Step 1: Evaluate the polynomial representation of A and B at the chosen points.
      vector<BigInteger> evalA = evaluatePoly(a, partSize);
      vector<BigInteger> evalB = evaluatePoly(b, partSize);
      
      // Step 2: Pointwise multiply the evaluations.
      // Use a container for the products at each evaluation point.
      vector<BigInteger> tempProducts(numPoints);
      for (int i = 0; i < numPoints; i++) {
          // The multiplication result at each point requires room for 2*segmentSize elements.
          // For pointwise multiplication, you can call the naive multiplication routine.
          tempProducts[i] = BigInteger(
            Multiply(evalA[i].GetInteger(), evalB[i].GetInteger(), base),
            evalA[i].IsNegative() ^ evalB[i].IsNegative());
      }
      
      // Step 3: Interpolate to combine the evaluated products back into the final result.
      BigInteger result = interpolate(tempProducts, partSize);

      return result.GetInteger();
    }

    private:
    // Evaluate the polynomial represented by 'poly' (a big number whose limbs
    // are stored in little–endian order) by splitting it into parts of length 'partSize'.
    //
    // Let a₀ = lower part, a₁ = middle part, a₂ = higher part (if available).
    // Then we compute:
    //   f(0)   = a₀
    //   f(1)   = a₀ + a₁ + a₂
    //   f(–1)  = a₀ – a₁ + a₂
    //   f(2)   = a₀ + 2*a₁ + 4*a₂
    //   f(∞)   = a₂
    //
    // The five evaluation results are stored in 'evalPoints', which is assumed
    // to have been preallocated with 5 vectors.
    vector<BigInteger> evaluatePoly(const vector<DataT>& poly, SizeT partSize)
    {
      vector<BigInteger> evalPoints(numPoints);
      // a0
      BigInteger a0(vector<DataT>(
          poly.begin(),
          poly.begin() + min(partSize, (SizeT)poly.size())
      ));

      // a1
      BigInteger a1;
      if (poly.size() > partSize) {
          a1 = BigInteger(vector<DataT>(
              poly.begin() + partSize,
              poly.begin() + min(2*partSize, (SizeT)poly.size())
          ));
      }

      // a2
      BigInteger a2;
      if (poly.size() > 2*partSize) {
          a2 = BigInteger(vector<DataT>(
              poly.begin() + 2*partSize,
              poly.end()
          ));
      }

      evalPoints[0] = a0;                     // f(0)
      evalPoints[1] = (a0 + a1 + a2);         // f(1)
      evalPoints[2] = (a0 + a2 - a1);         // f(-1)
      evalPoints[3] = a0 + (a1 * (DataT)2) + (a2 * (SizeT)4); // f(2)
      evalPoints[4] = a2;                     // f(∞)

      return evalPoints;
    }
    
    // Interpolates the evaluated products (contained in "products") to reconstruct
    // the final result vector after multiplication.
    // The five products correspond to the evaluations at:
    //    r0   = f(0)
    //    r1   = f(1)
    //    r-1  = f(–1)
    //    r2   = f(2)
    //    rInf = f(∞)
    //
    // Using the standard Toom-3 formulas, we set:
    //    c0 = r0,
    //    c4 = rInf,
    //    c2 = (r1 + r-1 − 2·r0 − 2·rInf) / 2,
    //    s  = (r1 − r-1) / 2,
    //    X  = (r2 − r0 − 4·c2 − 16·rInf) / 2,
    //    c3 = (X − s) / 3,
    //    c1 = s − c3.
    //
    // Then the final result is constructed as:
    //   result = c0 + (c1 << partSize) + (c2 << 2 * partSize)
    //          + (c3 << 3 * partSize) + (c4 << 4 * partSize)
    //
    // Note: The shifts are in limbs (each limb is 32 bits).
    BigInteger interpolate(vector<BigInteger> const& products, SizeT partSize)
    {
        // For clarity, assume the products are indexed as follows:
        //   products[0] = r0  (evaluation at 0)
        //   products[1] = r1  (evaluation at 1)
        //   products[2] = r-1 (evaluation at -1)
        //   products[3] = r2  (evaluation at 2)
        //   products[4] = rInf (evaluation at ∞)

        BigInteger r0 = products[0];
        BigInteger r1 = products[1];
        BigInteger rm1 = products[2];
        BigInteger r2 = products[3];
        BigInteger rInf = products[4];

        // c0 = r0 and c4 = rInf
        BigInteger c0 = r0;
        BigInteger c4 = rInf;

        // c2 = (r1 + r-1 - 2*c0 - 2*c4) / 2
        BigInteger c2 = ( (r1 + rm1) - ( (c0 * 2) + (c4 * 2) ) ) / (DataT)2;
        
        // s = (r1 - r-1) / 2
        BigInteger s = (r1 - rm1) / (DataT)2;

        // X = (r2 - c0 - 4*c2 - 16*c4) / 2
        BigInteger fourC2 = c2 * (DataT)4;
        BigInteger sixteenC4 = c4 * (DataT)16;
        BigInteger Xsum = c0 + fourC2 + sixteenC4;
        BigInteger X = (r2 - Xsum) / (DataT)2;

        // c3 = (X - s) / 3
        BigInteger c3 = (X - s) / (DataT)3; // (X - s) / 3

        // c1 = s - c3
        BigInteger c1 = s - c3; // s - c3

        // Combine coefficients into the final result:
        BigInteger res0 = c0;
        BigInteger res1 = c1 << partSize;
        BigInteger res2 = c2 << (2 * partSize);
        BigInteger res3 = c3 << (3 * partSize);
        BigInteger res4 = c4 << (4 * partSize);

        BigInteger res = res0 + res1 + res2 + res3 + res4;

        return res;
    }
   };
}

#endif
