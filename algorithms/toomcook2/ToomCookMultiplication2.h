/**
 * BigInteger Class
 * Version 10.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef TOOMCOOK_MULTIPLICATION_2
#define TOOMCOOK_MULTIPLICATION_2

#include <vector>
#include <stack>
#include <cassert>
using namespace std;

#include "../../BigInteger.h"
#include "../BigIntegerUtil.h"
#include "../classic/ClassicAddition.h"
#include "../classic/ClassicSubtraction.h"
#include "../classic/ClassicMultiplication.h"
#include "../classic/ClassicDivision.h"
#include "../classic/CommonAlgorithms.h"

namespace BigMath
{
  class ToomCookMultiplication2
  {
    private:
      // Pre-allocated buffer for temporary storage
      vector<DataT> tempBuffer;

    // Threshold for switching to the naive multiplication.
    static const SizeT baseCaseThreshold = 32;      
      
    public:
    vector<DataT> Multiply(
      vector<DataT> const& a,
      vector<DataT> const& b,
      ULong base,
      bool trim=true)
    {
      SizeT n = (SizeT)max(a.size(), b.size());

        if (n < baseCaseThreshold) {
            return ClassicMultiplication::Multiply(a, b, base);
        }      

      // For a Toom–3 algorithm, we typically split the input into 3 parts.
      // Calculate the part size (rounding up if necessary).
      SizeT partSize = (n + 2) / 3;
      
      // For evaluation, we need to compute 5 points (e.g., points for Toom–3: 0, 1, -1, 2, ∞).
      constexpr int numPoints = 5;
      // Determine the size required for each evaluation segment. It may be larger than
      // partSize if extra space is needed for carries in the arithmetic.
      SizeT segmentSize = partSize + 1; // adjust as needed
      
      // We now use the provided tempBuffer to simulate multiple temporary arrays.
      // We need numPoints segments for A and numPoints segments for B, i.e. 2*numPoints segments.
      SizeT requiredSize = 2 * numPoints * segmentSize;
      
      // Resize the buffer if needed
      if (tempBuffer.size() < requiredSize) {
        tempBuffer.resize(requiredSize, 0);
      }
      
      // Create containers for evaluated values. Instead of allocating new memory,
      // we assign each evaluation vector to a segment of the tempBuffer.
      // Here we use pointer arithmetic (or iterators) for clarity.
      vector< vector<DataT> > evalA(numPoints, vector<DataT>(segmentSize));
      vector< vector<DataT> > evalB(numPoints, vector<DataT>(segmentSize));
      
      // Optionally, you could have these evaluation vectors view memory directly within tempBuffer.
      // This example simply copies segments from the preallocated pool.
      for (int i = 0; i < numPoints; i++) {
          // For evalA, use segment [i*segmentSize, (i+1)*segmentSize)
          copy(tempBuffer.begin() + i * segmentSize,
                    tempBuffer.begin() + (i + 1) * segmentSize,
                    evalA[i].begin());
          // For evalB, use segment [ (numPoints + i)*segmentSize, (numPoints + i + 1)*segmentSize )
          copy(tempBuffer.begin() + (numPoints + i) * segmentSize,
                    tempBuffer.begin() + (numPoints + i + 1) * segmentSize,
                    evalB[i].begin());
      }
      
      // Step 1: Evaluate the polynomial representation of A and B at the chosen points.
      evaluatePoly(a, partSize, evalA, base);
      evaluatePoly(b, partSize, evalB, base);
      
      // Step 2: Pointwise multiply the evaluations.
      // Use a container for the products at each evaluation point.
      vector< vector<DataT> > tempProducts(numPoints);
      for (int i = 0; i < numPoints; i++) {
          // The multiplication result at each point requires room for 2*segmentSize elements.
          tempProducts[i].resize(2 * segmentSize);
          // For pointwise multiplication, you can call the naive multiplication routine.
          tempProducts[i] = Multiply(evalA[i], evalB[i], base);
      }
      
      // Step 3: Interpolate to combine the evaluated products back into the final result.
      vector<DataT> result = interpolate(tempProducts, partSize, base);

      if(trim)
      {
        BigIntegerUtil::TrimZeros(result);
      }
      
      return result;
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
    static void evaluatePoly(const vector<DataT>& poly,
                             SizeT partSize,
                             vector< vector<DataT> >& evalPoints,
                             ULong base)
    {
        // Split the input poly into three parts: a0, a1, and a2.
        vector<DataT> a0, a1, a2;

        // a0: limbs [0, partSize)
        a0.assign(poly.begin(), poly.begin() + (SizeT)min((SizeT)partSize, (SizeT)poly.size()));

        // a1: limbs [partSize, 2*partSize) if available; otherwise zero.
        if (poly.size() > partSize) {
            a1.assign(poly.begin() + partSize,
                      poly.begin() + (SizeT)min((SizeT)(2 * partSize), (SizeT) poly.size()));
        }
        else {
            // If poly does not have a middle part, represent a1 as zero.
            a1.assign(1, 0);
        }

        // a2: limbs [2*partSize, 3*partSize) if available; otherwise zero.
        if (poly.size() > 2 * partSize) {
            a2.assign(poly.begin() + 2 * partSize, poly.end());
        }
        else {
            a2.assign(1, 0);
        }

        // --- Evaluate at the five points ---

        // f(0) = a0.
        evalPoints[0] = a0;

        // f(1) = a0 + a1 + a2.
        {
            vector<DataT> tmp = ClassicAddition::Add(a1, a2, base);
            evalPoints[1] = ClassicAddition::Add(a0, tmp, base);
        }

        // f(–1) = a0 – a1 + a2 = (a0 + a2) – a1.
        {
            vector<DataT> tmp = ClassicAddition::Add(a0, a2, base);
            // NOTE: This implementation assumes tmp >= a1.
            evalPoints[2] = ClassicSubtraction::Subtract(tmp, a1, base);
        }

        // f(2) = a0 + 2*a1 + 4*a2.
        {
            vector<DataT> two_a1 = ClassicMultiplication::Multiply(a1, 2, base);
            vector<DataT> four_a2 = ClassicMultiplication::Multiply(a2, 4, base);
            vector<DataT> sum = ClassicAddition::Add(two_a1, four_a2, base);
            evalPoints[3] = ClassicAddition::Add(a0, sum, base);
        }

        // f(∞) = a2.
        evalPoints[4] = a2;
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
    static vector<DataT> interpolate(const vector< vector<DataT> >& products,
                            SizeT partSize,
                            ULong base)
    {
        // For clarity, assume the products are indexed as follows:
        //   products[0] = r0  (evaluation at 0)
        //   products[1] = r1  (evaluation at 1)
        //   products[2] = r-1 (evaluation at -1)
        //   products[3] = r2  (evaluation at 2)
        //   products[4] = rInf (evaluation at ∞)
        assert(products.size() == 5);

        const vector<DataT>& r0   = products[0];
        const vector<DataT>& r1   = products[1];
        const vector<DataT>& rm1  = products[2];
        const vector<DataT>& r2   = products[3];
        const vector<DataT>& rInf = products[4];

        // c0 = r0 and c4 = rInf
        vector<DataT> c0 = r0;
        vector<DataT> c4 = rInf;

        // c2 = (r1 + r-1 - 2*c0 - 2*c4) / 2
        vector<DataT> temp1 = ClassicAddition::Add(r1, rm1, base);               // r1 + r-1
        vector<DataT> twoC0 = ClassicMultiplication::Multiply(c0, 2, base);              // 2 * c0
        vector<DataT> twoC4 = ClassicMultiplication::Multiply(c4, 2, base);              // 2 * c4
        vector<DataT> temp2 = ClassicAddition::Add(twoC0, twoC4, base);            // 2*c0 + 2*c4
        vector<DataT> numeratorC2 = ClassicSubtraction::Subtract(temp1, temp2, base);       // r1 + r-1 - 2*c0 - 2*c4
        vector<DataT> c2 = ClassicDivision::Divide(numeratorC2, (DataT)2, base);         // Divide by 2

        // s = (r1 - r-1) / 2
        vector<DataT> numeratorS = ClassicSubtraction::Subtract(r1, rm1, base);
        vector<DataT> s = ClassicDivision::Divide(numeratorS, (DataT)2, base);

        // X = (r2 - c0 - 4*c2 - 16*c4) / 2
        vector<DataT> fourC2    = ClassicMultiplication::Multiply(c2, 4, base);
        vector<DataT> sixteenC4 = ClassicMultiplication::Multiply(c4, 16, base);
        vector<DataT> temp3 = ClassicAddition::Add(c0, fourC2, base);
        vector<DataT> temp4 = ClassicAddition::Add(temp3, sixteenC4, base);
        vector<DataT> numeratorX = ClassicSubtraction::Subtract(r2, temp4, base);
        vector<DataT> X = ClassicDivision::Divide(numeratorX, 2, base);

        // c3 = (X - s) / 3
        vector<DataT> numeratorC3 = ClassicSubtraction::Subtract(X, s, base);
        vector<DataT> c3 = ClassicDivision::Divide(numeratorC3, 3, base);

        // c1 = s - c3
        vector<DataT> c1 = ClassicSubtraction::Subtract(s, c3, base);

        // Combine coefficients into the final result:
        // result = c0 + (c1 << partSize) + (c2 << 2*partSize) + (c3 << 3*partSize) + (c4 << 4*partSize)
        vector<DataT> res0 = c0;
        vector<DataT> res1 = CommonAlgorithms::ShiftLeft(c1, partSize);
        vector<DataT> res2 = CommonAlgorithms::ShiftLeft(c2, 2 * partSize);
        vector<DataT> res3 = CommonAlgorithms::ShiftLeft(c3, 3 * partSize);
        vector<DataT> res4 = CommonAlgorithms::ShiftLeft(c4, 4 * partSize);

        vector<DataT> sum = ClassicAddition::Add(res0, res1, base);
        sum = ClassicAddition::Add(sum, res2, base);
        sum = ClassicAddition::Add(sum, res3, base);
        sum = ClassicAddition::Add(sum, res4, base);


        BigIntegerUtil::TrimZeros(sum);

        return sum;
    }
   };
}

#endif
