#ifndef QUOTIENT_SIZED_DIVISION
#define QUOTIENT_SIZED_DIVISION

#include <utility>
#include <vector>
using namespace std;

#include "../../common/Comparator.h"
#include "../../common/Util.h"
#include "../Addition.h"
#include "../Multiplication.h"
#include "../Subtraction.h"
#include "NewtonDivision.h"

namespace BigMath
{
  // Top-level dispatcher (defined in src/algorithms/Division.cpp). The tops
  // sub-division below routes through it so the recursion picks the right
  // strategy for its own size; its ratio is ~2, so it can never re-enter the
  // quotient-sized band (ratio < 4/3).
  pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
      vector<DataT> const &a,
      vector<DataT> const &b,
      BaseT base,
      bool computeRemainder);

  // Short-quotient division. When Δ = a.size() − b.size() is small relative
  // to the divisor, the quotient has only Δ+1 limbs and is determined to ±a
  // few units by the operands' TOPS: with t = Δ + GUARD,
  //   q_est = floor( (a >> B^(nb−t)) / (b >> B^(nb−t)) )
  // truncating b low limbs perturbs a/b by a relative ~B^(1−t), i.e. the
  // (Δ+1)-limb quotient by ≤ ~B^(Δ+2−t) = B^(2−GUARD) — sub-ulp. The exact
  // remainder then comes from ONE Δ×nb back-multiply plus ±few fixups.
  //
  // Cost: T(2Δ/Δ tops divide) + M(Δ, nb) + O(nb) — scales with the QUOTIENT,
  // not the divisor. This is the right shape for ratio ∈ (1, 4/3) at large b,
  // where full Newton pays an nb-sized reciprocal for a tiny quotient and
  // BZ's recursion blows up 7–128× on 2^k+1-family divisor sizes.
  class QuotientSizedDivision
  {
  public:
    static constexpr SizeT GUARD = 4;

    static pair<vector<DataT>, vector<DataT>> DivideAndRemainder(
        vector<DataT> const &a,
        vector<DataT> const &b,
        BaseT base,
        bool computeRemainder = true)
    {
      SizeT na = (SizeT)a.size();
      SizeT nb = (SizeT)b.size();
      // Caller contract (dispatch): a > b, and delta + GUARD < nb. The band
      // macros are tunable (-DBIGMATH_*), so enforce it instead of letting
      // `s` underflow below; fall back to Newton which handles any shape.
      if (na <= nb || (na - nb) + GUARD >= nb)
        return NewtonDivision::DivideAndRemainder(a, b, base, computeRemainder);
      SizeT delta = na - nb;
      SizeT t = delta + GUARD;
      SizeT s = nb - t;

      vector<DataT> a_top(a.begin() + s, a.end());
      vector<DataT> b_top(b.begin() + s, b.end());

      // Tops divide is ratio ~2. Post-#85-87 Newton beats BZ there from
      // ~6k limbs (measured: 12k limbs 13 vs 22 ms, 30k 40 vs 299 ms,
      // 65537 161 ms vs 5.3 s); below that the dispatcher's BZ/Fast pick
      // is near-tied.
      vector<DataT> q = (t >= 6144)
                            ? NewtonDivision::DivideAndRemainder(a_top, b_top, base, false).first
                            : BigMath::DivideAndRemainder(a_top, b_top, base, false).first;

      // Reconstruct: r = a − q·b, fixing q's ±few-unit estimate error.
      vector<DataT> qb = Multiply(q, b, base);

      static const vector<DataT> one{1};
      const int FIXUP_LIMIT = 8;

      int iters = 0;
      while (Compare(qb, a) > 0)
      {
        if (++iters > FIXUP_LIMIT)
          return NewtonDivision::DivideAndRemainder(a, b, base, computeRemainder);
        q = Subtract(q, one, base);
        qb = Subtract(qb, b, base);
      }
      vector<DataT> r = Subtract(a, qb, base);

      iters = 0;
      while (Compare(r, b) >= 0)
      {
        if (++iters > FIXUP_LIMIT)
          return NewtonDivision::DivideAndRemainder(a, b, base, computeRemainder);
        r = Subtract(r, b, base);
        q = Add(q, one, base);
      }

      TrimZerosToOne(q);
      if (!computeRemainder)
        return {q, vector<DataT>()};
      TrimZerosToOne(r);
      return {q, r};
    }
  };
}

#endif
