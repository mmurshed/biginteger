/**
 * BigInteger Class
 * Version 10.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef TOOMCOOK_MULTIPLICATION
#define TOOMCOOK_MULTIPLICATION

#include <vector>
#include <stack>
#include <cassert>
using namespace std;

#include "../../common/Util.h"
#include "../Addition.h"
#include "../Subtraction.h"
#include "../Shift.h"
#include "../multiplication/ClassicMultiplication.h"
#include "../multiplication/KaratsubaMultiplication.h"
#include "../division/ClassicDivision.h"

namespace BigMath
{
  class ToomCookMultiplication
  {
  private:
    // Threshold for switching to the naive multiplication.
    static const SizeT baseCaseThreshold = 256;

    struct SignedVector
    {
      vector<DataT> magnitude;
      int sign;

      SignedVector() : sign(1) {}
      SignedVector(vector<DataT> const &value, int valueSign = 1)
          : magnitude(value), sign(valueSign < 0 ? -1 : 1)
      {
        Normalize();
      }

      void Normalize()
      {
        TrimZeros(magnitude);
        if (IsZero(magnitude))
          sign = 1;
      }
    };

    // High-speed division helpers for interpolation
    static void DivideBy2Ptr(vector<DataT> &a, BaseT base)
    {
      if (a.empty()) return;
      if (base == Base2_32)
      {
        uint64_t carry = 0;
        for (int i = (int)a.size() - 1; i >= 0; --i)
        {
          uint64_t val = a[i] + (carry << 32);
          a[i] = (DataT)(val >> 1);
          carry = val & 1;
        }
      }
      else
      {
        ULong carry = 0;
        for (int i = (int)a.size() - 1; i >= 0; --i)
        {
          ULong val = a[i] + carry * base;
          a[i] = (DataT)(val / 2);
          carry = val % 2;
        }
      }
      TrimZeros(a);
    }

    static SignedVector Negate(SignedVector value)
    {
      if (!IsZero(value.magnitude))
        value.sign = -value.sign;
      return value;
    }

    static SignedVector AddSigned(SignedVector const &a, SignedVector const &b, BaseT base)
    {
      if (IsZero(a.magnitude))
        return b;
      if (IsZero(b.magnitude))
        return a;

      if (a.sign == b.sign)
        return SignedVector(Add(a.magnitude, b.magnitude, base), a.sign);

      Int cmp = Compare(a.magnitude, b.magnitude);
      if (cmp == 0)
        return SignedVector();
      if (cmp > 0)
        return SignedVector(Subtract(a.magnitude, b.magnitude, base), a.sign);
      return SignedVector(Subtract(b.magnitude, a.magnitude, base), b.sign);
    }

    static SignedVector SubtractSigned(SignedVector const &a, SignedVector const &b, BaseT base)
    {
      return AddSigned(a, Negate(b), base);
    }

    static SignedVector MultiplySmall(SignedVector const &a, DataT b, BaseT base)
    {
      if (b == 0 || IsZero(a.magnitude))
        return SignedVector();
      return SignedVector(ClassicMultiplication::Multiply(a.magnitude, b, base), a.sign);
    }

    static SignedVector DivideSmall(SignedVector value, DataT divisor, BaseT base)
    {
      if (divisor == 2)
        DivideBy2Ptr(value.magnitude, base);
      else if (divisor == 3)
        DivideBy3Ptr(value.magnitude, base);
      value.Normalize();
      return value;
    }

    static SignedVector ShiftSigned(SignedVector const &value, SizeT limbs)
    {
      if (IsZero(value.magnitude))
        return SignedVector();
      return SignedVector(ShiftLeft(value.magnitude, limbs), value.sign);
    }

    static void DivideBy3Ptr(vector<DataT> &a, BaseT base)
    {
      if (a.empty()) return;
      if (base == Base2_32)
      {
        uint64_t carry = 0;
        for (int i = (int)a.size() - 1; i >= 0; --i)
        {
          uint64_t val = a[i] + (carry << 32);
          a[i] = (DataT)(val / 3);
          carry = val % 3;
        }
      }
      else
      {
        ULong carry = 0;
        for (int i = (int)a.size() - 1; i >= 0; --i)
        {
          ULong val = a[i] + carry * base;
          a[i] = (DataT)(val / 3);
          carry = val % 3;
        }
      }
      TrimZeros(a);
    }

  private:
    static vector<DataT> MultiplyRecursive(
        vector<DataT> const &a,
        vector<DataT> const &b,
        vector<DataT> &tempBuffer,
        BaseT base)
    {
      if (IsZero(a) || IsZero(b)) // 0 times
        return vector<DataT>();

      // If b is a single digit, use the scalar multiplication
      if (b.size() == 1)
        return ClassicMultiplication::Multiply(a, b[0], base);
      // If a is a single digit, use the scalar multiplication
      if (a.size() == 1)
        return ClassicMultiplication::Multiply(b, a[0], base);

      SizeT n = (SizeT)max(a.size(), b.size());

      if (n < baseCaseThreshold)
      {
        return ClassicMultiplication::Multiply(a, b, base);
      }

      return KaratsubaMultiplication::Multiply(a, b, base);

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
      if (tempBuffer.size() < requiredSize)
      {
        tempBuffer.resize(requiredSize, 0);
      }

      // Create containers for evaluated values. Instead of allocating new memory,
      // we assign each evaluation vector to a segment of the tempBuffer.
      // Here we use pointer arithmetic (or iterators) for clarity.
      vector<vector<DataT>> evalA(numPoints, vector<DataT>(segmentSize));
      vector<vector<DataT>> evalB(numPoints, vector<DataT>(segmentSize));

      // Optionally, you could have these evaluation vectors view memory directly within tempBuffer.
      // This example simply copies segments from the preallocated pool.
      for (int i = 0; i < numPoints; i++)
      {
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
      int fm1SignA = evaluatePoly(a, partSize, evalA, base);
      int fm1SignB = evaluatePoly(b, partSize, evalB, base);

      // Step 2: Pointwise multiply the evaluations.
      // Use a container for the products at each evaluation point.
      vector<vector<DataT>> tempProducts(numPoints);
      vector<int> tempProductsSign(numPoints, 1);
      for (int i = 0; i < numPoints; i++)
      {
        // The multiplication result at each point requires room for 2*segmentSize elements.
        tempProducts[i].resize(2 * segmentSize);
        // For pointwise multiplication, you can call the naive multiplication routine.
        tempProducts[i] = MultiplyRecursive(evalA[i], evalB[i], tempBuffer, base);
      }

      int fm1Sign = fm1SignA * fm1SignB;

      // Step 3: Interpolate to combine the evaluated products back into the final result.
      vector<DataT> result = interpolate(tempProducts, fm1Sign, partSize, base);

      TrimZeros(result);

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
    static int evaluatePoly(const vector<DataT> &poly,
                            SizeT partSize,
                            vector<vector<DataT>> &evalPoints,
                            BaseT base)
    {
      // Split the input poly into three parts: a0, a1, and a2.
      vector<DataT> a0, a1, a2;

      // a0: limbs [0, partSize)
      a0.assign(poly.begin(), poly.begin() + (SizeT)min((SizeT)partSize, (SizeT)poly.size()));

      // a1: limbs [partSize, 2*partSize) if available; otherwise zero.
      if (poly.size() > partSize)
      {
        a1.assign(poly.begin() + partSize,
                  poly.begin() + (SizeT)min((SizeT)(2 * partSize), (SizeT)poly.size()));
      }
      else
      {
        // If poly does not have a middle part, represent a1 as zero.
        a1.assign(1, 0);
      }

      // a2: limbs [2*partSize, 3*partSize) if available; otherwise zero.
      if (poly.size() > 2 * partSize)
      {
        a2.assign(poly.begin() + 2 * partSize, poly.end());
      }
      else
      {
        a2.assign(1, 0);
      }

      // --- Evaluate at the five points ---

      // f(0) = a0.
      evalPoints[0] = a0;

      // f(1) = a0 + a1 + a2.
      {
        vector<DataT> tmp1 = Add(a1, a2, base);
        evalPoints[1] = Add(a0, tmp1, base);
      }

      // f(–1) = a0 – a1 + a2 = (a0 + a2) – a1.
      int fm1Sign = 1;
      {
        // cerr << a0 << " - " << a1 << " + " << a2;
        vector<DataT> tmp2 = Add(a0, a2, base);
        if (Compare(tmp2, a1) < 0)
        {
          fm1Sign = -1;
          evalPoints[2] = Subtract(a1, tmp2, base);
        }
        else
        {
          evalPoints[2] = Subtract(tmp2, a1, base);
        }
        // cerr << " = " << fm1Sign << " * " << evalPoints[2] << endl;
      }

      // f(2) = a0 + 2*a1 + 4*a2.
      {
        vector<DataT> two_a1 = ClassicMultiplication::Multiply(a1, 2, base);
        vector<DataT> four_a2 = ClassicMultiplication::Multiply(a2, 4, base);
        vector<DataT> sum = Add(two_a1, four_a2, base);
        evalPoints[3] = Add(a0, sum, base);
      }

      // f(∞) = a2.
      evalPoints[4] = a2;

      return fm1Sign;
    }

    // Interpolates the evaluated products (contained in "products") to reconstruct
    // the final result vector after multiplication.
    // The five products correspond to the evaluations at:
    //    r0   = f(0)
    //    r1   = f(1)
    //    rm1  = f(–1)
    //    r2   = f(2)
    //    rInf = f(∞)
    //
    // Using the standard Toom-3 formulas, we set:
    //    c0 = r0,
    //    c4 = rInf,
    //    c2 = (r1 + rm1 − 2·r0 − 2·rInf) / 2 = ((r1 + rm1) − 2(r0 + rInf)) / 2,
    //    s  = (r1 − rm1) / 2,
    //    X  = (r2 − r0 − 4·c2 − 16·rInf) / 2,
    //    c3 = (X − s) / 3,
    //    c1 = s − c3.
    //
    // Then the final result is constructed as:
    //   result = c0 + (c1 << partSize) + (c2 << 2 * partSize)
    //          + (c3 << 3 * partSize) + (c4 << 4 * partSize)
    //
    // Note: The shifts are in limbs (each limb is 32 bits).
    static vector<DataT> interpolate(const vector<vector<DataT>> &products,
                                     int rm1Sign,
                                     SizeT partSize,
                                     BaseT base)
    {
      // For clarity, assume the products are indexed as follows:
      //   products[0] = r0  (evaluation at 0)
      //   products[1] = r1  (evaluation at 1)
      //   products[2] = r-1 = rm1 (evaluation at -1)
      //   products[3] = r2  (evaluation at 2)
      //   products[4] = rInf (evaluation at ∞)
      assert(products.size() == 5);

      SignedVector v0(products[0]);
      SignedVector v1(products[1]);
      SignedVector v2(products[2], rm1Sign);
      SignedVector v3(products[3]);
      SignedVector v4(products[4]);

      v3 = DivideSmall(SubtractSigned(v3, v1, base), 3, base);
      v1 = DivideSmall(SubtractSigned(v1, v2, base), 2, base);
      v2 = SubtractSigned(v2, v0, base);
      v3 = AddSigned(DivideSmall(SubtractSigned(v2, v3, base), 2, base),
                     MultiplySmall(v4, 2, base),
                     base);
      v2 = SubtractSigned(AddSigned(v2, v1, base), v4, base);
      v1 = SubtractSigned(v1, v3, base);

      SignedVector result = v0;
      result = AddSigned(result, ShiftSigned(v1, partSize), base);
      result = AddSigned(result, ShiftSigned(v2, 2 * partSize), base);
      result = AddSigned(result, ShiftSigned(v3, 3 * partSize), base);
      result = AddSigned(result, ShiftSigned(v4, 4 * partSize), base);

      if (result.sign == -1)
        throw runtime_error("Negative result in multiplication.");

      return result.magnitude;
    }

  public:
    static vector<DataT> Multiply(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base)
    {
      // Pre-allocated buffer for temporary storage
      vector<DataT> tempBuffer;
      return MultiplyRecursive(a, b, tempBuffer, base);
    }
  };
}

#endif
