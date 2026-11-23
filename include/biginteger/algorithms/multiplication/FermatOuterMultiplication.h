/**
 * BigMath: outer FFT over Z/(2^n + 1) for huge 64-bit-limb products.
 *
 * A zero-padded cyclic convolution multiplies the operand pieces. The
 * Fermat modulus makes twiddles power-of-two shifts and is large enough
 * to recover each integer coefficient without wraparound. Pointwise
 * products use the ordinary BigMath multiplier, followed by a fold.
 */
#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

#include "../../common/LimbKernels.h"
#include "../../common/Parallel.h"
#include "../../common/Util.h"

namespace BigMath::detail
{
  // Bounded 64-bit product for Fermat residues. Inputs can be the same view.
  void MultiplyFermatPoint(std::span<const DataT> a, std::span<const DataT> b,
                           std::vector<DataT> &product, std::vector<DataT> &scratch);
}

namespace BigMath::FermatOuter
{
  constexpr size_t MaxRingBits = size_t{1} << 28;
  constexpr unsigned MaxExponent = 14;

  // A ring element has L low limbs and a top limb. In canonical form the
  // top limb is zero or represents -1 (one with all L low limbs zero).
  inline void Normalize(std::span<DataT> x)
  {
    const DataT top = x.back();
    x.back() = 0;
    if (LimbKernel::SubBorrow(x.data(), x.data(), (SizeT)(x.size() - 1), top))
      x.back() = LimbKernel::AddCarry(x.data(), x.data(), (SizeT)(x.size() - 1), 1);
  }

  inline void AddMod(std::span<DataT> out, std::span<const DataT> a, std::span<const DataT> b)
  {
    const DataT top = a.back() + b.back();
    const SizeT n = (SizeT)(out.size() - 1);
    out.back() = top + LimbKernel::AddN(out.data(), a.data(), b.data(), n);
    Normalize(out);
  }

  inline void SubMod(std::span<DataT> out, std::span<const DataT> a, std::span<const DataT> b)
  {
    const Int top = (Int)a.back() - (Int)b.back();
    const SizeT n = (SizeT)(out.size() - 1);
    const Int high = top - (Int)LimbKernel::SubN(out.data(), a.data(), b.data(), n);
    out.back() = high >= 0 ? (DataT)high
                           : LimbKernel::AddCarry(out.data(), out.data(), n, (DataT)-high);
    Normalize(out);
  }

  inline void NegateMod(std::span<DataT> out, std::span<const DataT> x)
  {
    if (x.back())
    {
      std::fill(out.begin(), out.end(), 0);
      out.front() = 1;
      return;
    }
    DataT borrow = 0;
    for (size_t i = 0; i + 1 < x.size(); ++i)
      out[i] = LimbKernel::Sbb(0, x[i], borrow, borrow);
    out.back() = borrow ? LimbKernel::AddCarry(out.data(), out.data(), (SizeT)(out.size() - 1), 1) : 0;
  }

  // out and x are disjoint. Compute the two n-bit halves of the shifted
  // value on demand, since (low + high * 2^n) reduces to low - high.
  inline void MulPow2Mod(std::span<DataT> out, std::span<const DataT> x, size_t shift)
  {
    const size_t limbs = x.size() - 1;
    const size_t bits = limbs * 64;
    const bool negative = shift >= bits;
    if (negative) shift -= bits;
    if (shift == 0)
    {
      if (negative) NegateMod(out, x);
      else std::copy(x.begin(), x.end(), out.begin());
      return;
    }

    const size_t wordShift = shift / 64;
    const unsigned bitShift = (unsigned)(shift % 64);
    const auto shifted = [&](size_t i) {
      if (i < wordShift) return DataT{0};
      i -= wordShift;
      DataT word = i < x.size() ? x[i] << bitShift : 0;
      if (bitShift && i && i <= x.size()) word |= x[i - 1] >> (64 - bitShift);
      return word;
    };
    DataT borrow = 0;
    for (size_t i = 0; i < limbs; ++i)
    {
      const DataT low = shifted(i), high = shifted(i + limbs);
      out[i] = negative ? LimbKernel::Sbb(high, low, borrow, borrow)
                        : LimbKernel::Sbb(low, high, borrow, borrow);
    }
    out.back() = borrow ? LimbKernel::AddCarry(out.data(), out.data(), (SizeT)limbs, 1) : 0;
  }

  struct PointwiseScratch
  {
    std::vector<DataT> product, recursive;
  };

  inline void MulMod(std::span<DataT> out, std::span<const DataT> a,
                     std::span<const DataT> b, PointwiseScratch &scratch)
  {
    if (a.back() || b.back())
    {
      NegateMod(out, a.back() ? b : a);
      return;
    }
    a = a.first(LimbDetail::TrimmedLimbCount(a.first(a.size() - 1)));
    b = a.data() == b.data() ? a : b.first(LimbDetail::TrimmedLimbCount(b.first(b.size() - 1)));
    // Finish reading both inputs before folding into the possibly aliased output.
    auto &product = scratch.product;
    detail::MultiplyFermatPoint(a, b, product, scratch.recursive);
    const SizeT n = (SizeT)(out.size() - 1);
    const size_t low = std::min<size_t>(n, product.size());
    std::copy_n(product.data(), low, out.data());
    std::fill(out.begin() + low, out.end(), 0);
    if (product.size() > n && LimbKernel::SubInPlace(out.data(), n, product.data() + n,
                                                    (SizeT)(product.size() - n)))
      out.back() = LimbKernel::AddCarry(out.data(), out.data(), n, 1);
  }

  struct Parameters
  {
    unsigned exponent;
    size_t count;       // K
    size_t pieceLimbs;  // m
    size_t ringLimbs;   // n' / 64, excluding the top limb
  };

  inline size_t DivideUp(size_t value, size_t divisor)
  {
    return value / divisor + (value % divisor != 0);
  }

  inline std::optional<Parameters> MakeParameters(size_t a, size_t b, unsigned exponent,
                                                   size_t pieceLimbs = 0)
  {
    constexpr size_t limit = std::numeric_limits<Int>::max();
    if (!a || !b || a > limit || b > limit - a || exponent < 1 || exponent > MaxExponent)
      return std::nullopt;
    const size_t count = size_t{1} << exponent;
    if (!pieceLimbs) pieceLimbs = DivideUp(a + b, count - 1);
    if (pieceLimbs > (MaxRingBits - exponent - 1) / 128 ||
        DivideUp(a, pieceLimbs) + DivideUp(b, pieceLimbs) - 1 > count)
      return std::nullopt;

    // Both 64 and K/2 are powers of two, so their LCM is the larger one.
    const size_t alignment = std::max<size_t>(64, count / 2);
    const size_t bits = DivideUp(128 * pieceLimbs + exponent + 1, alignment) * alignment;
    if (bits > MaxRingBits) return std::nullopt;
    return Parameters{exponent, count, pieceLimbs, bits / 64};
  }

  inline std::optional<Parameters> SelectParameters(size_t a, size_t b)
  {
    for (unsigned exponent = 6; exponent <= MaxExponent; ++exponent)
      if (auto params = MakeParameters(a, b, exponent)) return params;
    return std::nullopt;
  }

  class Transform
  {
    Parameters params;
    size_t stride;
    size_t tasks;
    // Allocate outside the pool callback. One scratch element per task keeps
    // temporary storage off the worker stacks, including at full ring sizes.
    std::vector<DataT> scratch;

    template <class Body>
    void ForRanges(size_t count, const Body &run)
    {
      if (tasks == 1)
      {
        run(0, count, std::span<DataT>(scratch));
        return;
      }
      // ParallelFor's cutoff counts items. Each item here processes a whole
      // ring element, so use coarse ParallelDo over pre-split ranges.
      ParallelDo((Int)tasks, [&](Int start, Int end) {
        for (size_t i = start; i < (size_t)end; ++i)
          run(count * i / tasks, count * (i + 1) / tasks,
              std::span<DataT>(scratch).subspan(i * stride, stride));
      });
    }

    void Layer(std::span<DataT> plane, size_t width, bool inverse)
    {
      const size_t half = width / 2;
      const size_t order = params.ringLimbs * 128;
      const size_t step = order / width;
      ForRanges(params.count / 2, [&](size_t begin, size_t end, std::span<DataT> temp) {
        for (size_t i = begin; i < end; ++i)
        {
          const size_t j = i % half;
          const size_t base = i / half * width + j;
          auto x = plane.subspan(base * stride, stride);
          auto y = plane.subspan((base + half) * stride, stride);
          const size_t shift = step * j;
          if (inverse)
          {
            MulPow2Mod(temp, y, shift ? order - shift : 0);
            SubMod(y, x, temp);
            AddMod(x, x, temp);
          }
          else
          {
            SubMod(temp, x, y);
            AddMod(x, x, y);
            MulPow2Mod(y, temp, shift);
          }
        }
      });
    }

  public:
    explicit Transform(Parameters p)
        : params(p), stride(p.ringLimbs + 1),
          tasks((p.count / 2) * stride >= ParallelMinSize()
                    ? std::min<size_t>(ParallelNumThreads(), p.count / 2) : 1),
          scratch(tasks * stride)
    {
    }

    void Forward(std::span<DataT> plane)
    {
      assert(plane.size() == params.count * stride);
      for (size_t width = params.count; width > 1; width /= 2)
        Layer(plane, width, false);
    }

    void Inverse(std::span<DataT> plane)
    {
      assert(plane.size() == params.count * stride);
      for (size_t width = 2; width <= params.count; width *= 2)
        Layer(plane, width, true);
      const size_t scale = params.ringLimbs * 128 - params.exponent;
      ForRanges(params.count, [&](size_t begin, size_t end, std::span<DataT> temp) {
        for (size_t i = begin; i < end; ++i)
        {
          auto value = plane.subspan(i * stride, stride);
          MulPow2Mod(temp, value, scale);
          std::copy(temp.begin(), temp.end(), value.begin());
        }
      });
    }
  };

  inline std::vector<DataT> Pack(std::span<const DataT> value, const Parameters &p)
  {
    const size_t stride = p.ringLimbs + 1;
    std::vector<DataT> plane(p.count * stride);
    for (size_t offset = 0, i = 0; offset < value.size(); offset += p.pieceLimbs, ++i)
      std::copy_n(value.data() + offset, std::min(p.pieceLimbs, value.size() - offset),
                  plane.data() + i * stride);
    return plane;
  }

  inline std::vector<DataT> MultiplyFermat(const std::vector<DataT> &a, const std::vector<DataT> &b,
                                           const Parameters &p)
  {
    if (a.empty() || b.empty()) return {};
    const auto expected = MakeParameters(a.size(), b.size(), p.exponent, p.pieceLimbs);
    if (!p.pieceLimbs || !expected || expected->count != p.count || expected->ringLimbs != p.ringLimbs)
      throw std::invalid_argument("invalid Fermat outer product parameters");

    const size_t stride = p.ringLimbs + 1;
    auto left = Pack(a, p);
    Transform transform(p);
    transform.Forward(left);
    {
      const bool square = &a == &b;
      std::vector<DataT> right;
      if (!square)
      {
        right = Pack(b, p);
        transform.Forward(right);
      }
      PointwiseScratch scratch;
      // The inner dispatcher uses the same pool. Do not nest it in a
      // butterfly task. Its products fit below the CRT transform ceiling
      // because n' is capped at 2^28 bits.
      for (size_t i = 0; i < p.count; ++i)
      {
        auto x = std::span<DataT>(left).subspan(i * stride, stride);
        std::span<const DataT> y = square ? x : std::span<DataT>(right).subspan(i * stride, stride);
        MulMod(x, x, y, scratch);
      }
    } // Release the second plane before allocating the product.
    transform.Inverse(left);

    std::vector<DataT> product(a.size() + b.size());
    const size_t coefficients = DivideUp(a.size(), p.pieceLimbs) + DivideUp(b.size(), p.pieceLimbs) - 1;
    for (size_t i = 0; i < coefficients; ++i)
    {
      const size_t offset = i * p.pieceLimbs;
      const size_t room = product.size() - offset;
      auto term = std::span<const DataT>(left).subspan(i * stride, stride);
      const size_t used = std::min(room, p.ringLimbs);
      // True coefficients fit below 2^n'. Padding beyond the product end
      // is zero, even for a partial last piece.
      if (std::any_of(term.begin() + used, term.end(), [](DataT x) { return x != 0; }) ||
          LimbKernel::AddInPlace(product.data() + offset, (SizeT)room, term.data(), (SizeT)used))
        throw std::logic_error("Fermat outer product exceeds its coefficient bound");
    }
    TrimZeros(product);
    return product;
  }
}
