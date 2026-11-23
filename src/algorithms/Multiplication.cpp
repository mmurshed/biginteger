/**
 * BigMath: Multiplication dispatcher implementation.
 *
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#include "biginteger/algorithms/Multiplication.h"

#include <algorithm>
#include "biginteger/algorithms/Squaring.h"
#include "biginteger/algorithms/multiplication/FermatOuterMultiplication.h"
#include "biginteger/common/Parallel.h"

namespace BigMath
{
  const SizeT CLASSIC_MULTIPLICATION_THRESHOLD = BIGMATH_CLASSIC_MULTIPLICATION_THRESHOLD;
  const SizeT NTT_MULTIPLICATION_THRESHOLD = BIGMATH_NTT_MULTIPLICATION_THRESHOLD;
  const SizeT NTT_MIN_LIMB_THRESHOLD = BIGMATH_NTT_MIN_LIMB_THRESHOLD;
  const ULong NTT_AREA_THRESHOLD = BIGMATH_NTT_AREA_THRESHOLD;
  const ULong NTT_COLD_AREA_THRESHOLD = BIGMATH_NTT_COLD_AREA_THRESHOLD;
  const SizeT CLASSIC_MIN_LIMB_THRESHOLD = BIGMATH_CLASSIC_MIN_LIMB_THRESHOLD;
  const SizeT CLASSIC_SKEW_MIN_LIMB_THRESHOLD = BIGMATH_CLASSIC_SKEW_MIN_LIMB_THRESHOLD;
  const SizeT CLASSIC_SKEW_RATIO = BIGMATH_CLASSIC_SKEW_RATIO;

  namespace
  {
    enum class ProductPath { Recursive, ColdRecursive, Transform };

    ProductPath ChooseProductPath(SizeT shorter, SizeT longer, [[maybe_unused]] bool serialHoldback)
    {
      // Sliced Karatsuba stays competitive on skewed inputs with a small
      // shorter operand. The area test accounts for both operand lengths.
      const ULong area = (ULong)shorter * longer;
      if (shorter + longer < NTT_MULTIPLICATION_THRESHOLD ||
          shorter < NTT_MIN_LIMB_THRESHOLD || area < NTT_AREA_THRESHOLD)
        return ProductPath::Recursive;
#if !BIGMATH_USE_THREADS && BIGMATH_NTT_SERIAL_BALANCED_MIN_LIMBS > 0
      if (serialHoldback && longer < BIGMATH_NTT_SERIAL_BALANCED_MIN_LIMBS && longer - shorter <= 1)
        return ProductPath::Recursive;
#endif
      // A parked pool may make a small NTT slower than recursion. Prewarm
      // for the next product if these operations are arriving in a burst.
      return area >= NTT_COLD_AREA_THRESHOLD || ParallelWarmOrPrewarm()
          ? ProductPath::Transform : ProductPath::ColdRecursive;
    }

    ProductPath ChooseSquarePath(SizeT size)
    {
      if (size < NTT_SQUARE_THRESHOLD)
        return ProductPath::Recursive;
      return size >= NTT_SQUARE_COLD_THRESHOLD || ParallelWarmOrPrewarm()
          ? ProductPath::Transform : ProductPath::ColdRecursive;
    }
  }

  std::vector<DataT> Multiply(std::vector<DataT> const &a,
                              std::vector<DataT> const &b,
                              BaseT base)
  {
    return detail::MultiplyLimbs(a, b, base);
  }

  std::vector<DataT> detail::MultiplyLimbs(std::span<const DataT> a,
                                         std::span<const DataT> b, BaseT base)
  {
    if (IsZero(a) || IsZero(b))
      return std::vector<DataT>{0};

    if (b.size() == 1)
      return ClassicMultiplication::Multiply(a, b[0], base);
    if (a.size() == 1)
      return ClassicMultiplication::Multiply(b, a[0], base);

    // Below the measured small-square ceiling we can exploit symmetry before
    // general multiplication dispatch. CRT already recognizes aliases at scale.
    const bool sameOperand = a.data() == b.data() && a.size() == b.size();
    if (sameOperand && base == Base2_64 && a.size() <= BIGMATH_SMALL_SELF_SQUARE_MAX)
      return SquareLimbs(a, base);

    SizeT minSize = std::min(a.size(), b.size());
    SizeT maxSize = std::max(a.size(), b.size());

    if (ClassicProductShape(minSize, maxSize))
      return ClassicMultiplication::MultiplyLimbs(a, b, base);

    const auto path = ChooseProductPath(minSize, maxSize, base == Base2_64 && !sameOperand);
    if (path == ProductPath::Transform)
      return NTTMultiplication::MultiplyLimbs(a, b, base);
    if (path == ProductPath::ColdRecursive)
    {
      std::vector<DataT> product = (sameOperand && base == Base2_64)
          ? KaratsubaSquare::SquareLimbs(a, base)
          : KaratsubaMultiplication::MultiplyLimbs(a, b, base);
      ParallelNoteColdDone();
      return product;
    }

    return KaratsubaMultiplication::MultiplyLimbs(a, b, base);
  }

  void detail::MultiplyFermatPoint(std::span<const DataT> a, std::span<const DataT> b,
                                  std::vector<DataT> &product, std::vector<DataT> &scratch)
  {
    if (IsZero(a) || IsZero(b))
    {
      product.clear();
      return;
    }
    if (a.size() < b.size())
      std::swap(a, b);

    const SizeT la = (SizeT)a.size(), lb = (SizeT)b.size();
    const bool sameOperand = a.data() == b.data() && la == lb;
    const bool smallSquare = sameOperand && la <= BIGMATH_SMALL_SELF_SQUARE_MAX;
    if (lb == 1 || (!smallSquare && ClassicProductShape(lb, la)))
    {
      product.resize(la + lb);
      ClassicMultiplication::MultiplyInto(a, b, product.data(), Base2_64);
      TrimZeros(product);
      return;
    }

    const auto path = smallSquare ? ChooseSquarePath(la) : ChooseProductPath(lb, la, !sameOperand);
    if (path == ProductPath::Transform)
    {
#if BIGMATH_NTT_CRT
      if (smallSquare || la + lb >= BIGMATH_NTT_CRT_THRESHOLD)
      {
        // The Fermat ring limit guarantees that a point product fits one
        // standard CRT transform, even at the maximum allowed residue width.
        NttCrt::MultiplyWithinTransform(a, b, Base2_64, (ULong)(la + lb) * 2 - 1, product);
        return;
      }
#endif
      product = smallSquare ? NTTSquare::SquareLimbs(a, Base2_64)
                            : NTTMultiplication::MultiplyLimbs(a, b, Base2_64);
      return;
    }

    product.resize(la + lb);
    if (smallSquare || (sameOperand && path == ProductPath::ColdRecursive))
    {
      scratch.resize(KaratsubaSquare::ScratchLimbs(la));
      KaratsubaSquare::SquareInto(a.data(), la, product.data(), scratch.data(), Base2_64);
    }
    else
    {
      scratch.resize(KaratsubaMultiplication::ScratchLimbs(la, lb));
      KaratsubaMultiplication::MultiplyInto(a.data(), la, b.data(), lb,
                                           product.data(), scratch.data(), Base2_64);
    }
    TrimZeros(product);
    if (path == ProductPath::ColdRecursive)
      ParallelNoteColdDone();
  }

  std::vector<DataT> Multiply(std::vector<DataT> const &a,
                              DataT b,
                              BaseT base)
  {
    if (b == 0 || IsZero(a))
      return std::vector<DataT>{0};
    return ClassicMultiplication::Multiply(a, b, base);
  }
}
