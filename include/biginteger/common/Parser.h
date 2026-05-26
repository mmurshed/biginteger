/**
 * BigMath: Decimal parsing and formatting (Parser + ToString).
 *
 * Architectural overview lives in STRING_CONVERSION.md. Brief:
 *
 *   Parsing:
 *     ParseUnsignedLinear      — 18-digit chunk accumulation, O(L²) but tight
 *     ParseUnsignedDivideConquer — D&C using cached Pow10, O(M(L) · log L)
 *
 *   Formatting:
 *     ToStringLinearAppend     — divmod-10^18 loop, O(L²)
 *     ToStringDivConquer       — D&C using Newton-Divider chain, O(M(L) · log L)
 *
 * Thresholds:
 *     DecimalDcThreshold       — parse linear→D&C cutoff (8192 digits)
 *     ToStringDcThreshold      — format linear→D&C cutoff (2048 digits)
 *
 * Caches (`Pow10`) and chain-builder helpers live in the .cpp with internal-linkage
 * `static thread_local` storage that is shared across all callers in the same TU.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER_PARSER
#define BIGINTEGER_PARSER

#include <string>
#include <vector>

#include "Util.h"
#include "../BigInteger.h"

namespace BigMath
{
  // Base 10^18 grouping. 10^18 < 2^60, fits in uint64 with room for 32-bit limb multiply.
  inline constexpr ULong Base10_18 = 1000000000000000000ULL;
  inline constexpr SizeT Base10_18_Zeroes = 18;
  inline constexpr SizeT DecimalDcThreshold = 8192;

#ifndef BIGMATH_TOSTR_DC_THRESHOLD
#define BIGMATH_TOSTR_DC_THRESHOLD 2048
#endif
  inline constexpr SizeT ToStringDcThreshold = BIGMATH_TOSTR_DC_THRESHOLD;

  // ─── parsing ────────────────────────────────────────────────────────────────
  Int FindRange(char const *num, Int &start, Int &end, bool &isNegative);
  ULong ParseChunk(char const *num, Int start, SizeT len);

  // Direct conversion of a 64-bit value into base-2^32 limbs.
  std::vector<DataT> Convert(ULong n);

  // Pow10: memoized recursive doubling. Even d uses Square(Pow10(d/2)) to save one FFT.
  std::vector<DataT> Pow10(SizeT digits);

  std::vector<DataT> ParseUnsignedLinear(char const *num, Int start, Int end);
  std::vector<DataT> ParseUnsignedDivideConquer(char const *num, Int start, Int end);
  std::vector<DataT> ParseUnsigned(char const *num, Int start, Int end);

  BigInteger Parse(char const *num, Int start, Int *char_processed = nullptr);
  BigInteger Parse(char const *num, Int *char_processed = nullptr);

  // ─── formatting ─────────────────────────────────────────────────────────────
  // Appends to `out`. padTo > 0 pads with leading zeros to exactly padTo digits.
  void ToStringLinearAppend(std::vector<DataT> r, SizeT padTo, std::string &out);

  std::string ToString(std::vector<DataT> const &bigInt, SizeT start, SizeT end, bool isNeg = false);
  std::string ToString(std::vector<DataT> const &bigInt, bool isNeg = false);
  std::string ToString(BigInteger const &bigInt);
}

#endif
