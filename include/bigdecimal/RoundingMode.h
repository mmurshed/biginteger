/**
 * Rounding modes for BigDecimal division and rescaling.
 *
 * Semantics mirror java.math.RoundingMode.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGDECIMAL_ROUNDING_MODE
#define BIGDECIMAL_ROUNDING_MODE

namespace BigMath
{
  enum class RoundingMode
  {
    UP,            // away from zero
    DOWN,          // toward zero (truncate)
    CEILING,       // toward +infinity
    FLOOR,         // toward -infinity
    HALF_UP,       // nearest; ties away from zero
    HALF_DOWN,     // nearest; ties toward zero
    HALF_EVEN,     // nearest; ties to even (banker's)
    UNNECESSARY    // throws if rounding would discard a non-zero remainder
  };
}

#endif
