/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_CONSTANTS
#define BIGINTEGER_CONSTANTS

namespace BigMath
{
  typedef int64_t Long;
  typedef uint64_t ULong;

  typedef int32_t Int;
  typedef uint32_t UInt;

  typedef uint64_t DataT;
  typedef uint32_t SizeT;
  typedef int64_t BaseT;

    // Base 10
    const BaseT Base10 = 10;
    // Base 10^3
    const BaseT Base100 = 100;
    // Base 10^8
    const BaseT Base100M = 100000000lu;
    // Number of digits in `BASE'
    const SizeT Base100_Zeroes = 2;

    // Base 2
    const BaseT Base2 = 2;
    // Base 16
    const BaseT Base16 = 16;
    // Base 2^8
    const BaseT Base2_8 = 256;
    // Base 2^10
    const BaseT Base2_10 = 1024;
    // Base 2^16
    const BaseT Base2_16 = 65536ul;
    // Base 2^32
    const BaseT Base2_32 = 4294967296ul; // 2^32
    // Base 2^31
    const BaseT Base2_31 = 2147483648; // 2^31
}

#endif
