/**
 * BigMath: Goldilocks-NTT result finalization — carry-propagate 16-bit
 * convolution coefficients and pack them into base-2^32 or base-2^64 limbs.
 *
 * Shared by NTTMultiplication and NTTSquare, which previously carried
 * duplicated (and, for the 2^64 case, divergent) private copies.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef NTT_FINALIZE
#define NTT_FINALIZE

#include <vector>

#include "../../common/Util.h"

namespace BigMath
{
  // 2-into-1: pack pairs of carry-propagated 16-bit coefficients into
  // base-2^32 limbs.
  inline std::vector<DataT> NttFinalizeBase2_32(const std::vector<ULong> &coeffs, SizeT coeffCount)
  {
    std::vector<DataT> result;
    result.reserve(coeffCount / 2 + 2);

    ULong carry = 0;
    ULong low = 0;
    bool hasLow = false;

    for (SizeT i = 0; i < coeffCount; ++i)
    {
      ULong total = coeffs[i] + carry;
      ULong digit = total & 0xFFFFULL;
      carry = total >> 16;

      if (hasLow)
      {
        result.push_back((DataT)(low | (digit << 16)));
        hasLow = false;
      }
      else
      {
        low = digit;
        hasLow = true;
      }
    }

    while (carry > 0)
    {
      ULong digit = carry & 0xFFFFULL;
      carry >>= 16;

      if (hasLow)
      {
        result.push_back((DataT)(low | (digit << 16)));
        hasLow = false;
      }
      else
      {
        low = digit;
        hasLow = true;
      }
    }

    if (hasLow)
      result.push_back((DataT)low);

    TrimZeros(result);
    return result;
  }

  // 4-into-1: pack quadruples of carry-propagated 16-bit coefficients into
  // base-2^64 limbs. Unrolled main loop + slot rotor for the tail.
  inline std::vector<DataT> NttFinalizeBase2_64(const std::vector<ULong> &coeffs, SizeT coeffCount)
  {
    std::vector<DataT> result;
    result.reserve(coeffCount / 4 + 2);

    ULong carry = 0;
    SizeT i = 0;
    for (; i + 3 < coeffCount; i += 4)
    {
      ULong total0 = coeffs[i] + carry;
      ULong d0 = total0 & 0xFFFFULL;
      carry = total0 >> 16;

      ULong total1 = coeffs[i + 1] + carry;
      ULong d1 = total1 & 0xFFFFULL;
      carry = total1 >> 16;

      ULong total2 = coeffs[i + 2] + carry;
      ULong d2 = total2 & 0xFFFFULL;
      carry = total2 >> 16;

      ULong total3 = coeffs[i + 3] + carry;
      ULong d3 = total3 & 0xFFFFULL;
      carry = total3 >> 16;

      result.push_back((DataT)(d0 | (d1 << 16) | (d2 << 32) | (d3 << 48)));
    }

    ULong limb_acc = 0;
    int slot = 0;
    for (; i < coeffCount; ++i)
    {
      ULong total = coeffs[i] + carry;
      ULong digit = total & 0xFFFFULL;
      carry = total >> 16;
      limb_acc |= digit << (slot * 16);
      ++slot;
    }
    while (carry > 0)
    {
      ULong digit = carry & 0xFFFFULL;
      carry >>= 16;

      limb_acc |= digit << (slot * 16);
      ++slot;
      if (slot == 4)
      {
        result.push_back((DataT)limb_acc);
        limb_acc = 0;
        slot = 0;
      }
    }

    if (slot != 0)
      result.push_back((DataT)limb_acc);

    TrimZeros(result);
    return result;
  }
}

#endif
