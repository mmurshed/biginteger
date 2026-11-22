/**
 * BigInteger Class
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef BIGINTEGER
#define BIGINTEGER

#include <algorithm>
#include <cstdint>
#include <span>
#include <vector>

#include "common/Util.h"
#include "common/Comparator.h"

namespace BigMath
{
#if defined(__has_cpp_attribute)
#if __has_cpp_attribute(clang::lifetimebound)
#define BIGMATH_LIFETIMEBOUND [[clang::lifetimebound]]
#endif
#endif
#ifndef BIGMATH_LIFETIMEBOUND
#define BIGMATH_LIFETIMEBOUND
#endif

// Limbs stored inside the BigInteger object itself. Values of up to this
// many limbs (256 bits with 64-bit limbs) are created, copied and destroyed
// without touching the heap; longer values live in a std::vector.
#ifndef BIGMATH_INLINE_LIMBS
#define BIGMATH_INLINE_LIMBS 4
#endif

  class BigInteger
  {
  public:
    static constexpr SizeT InlineLimbs = BIGMATH_INLINE_LIMBS;
    static_assert(InlineLimbs >= 2, "BIGMATH_INLINE_LIMBS must be at least 2");

    // A `std::vector<DataT> const &` view of the limbs for vector-based
    // algorithms. Refers to the stored vector when the value is on the heap;
    // copies the (at most InlineLimbs) limbs otherwise. The reference lives
    // only as long as the VectorRef: pass `x.AsVector()` straight to a call,
    // or keep the VectorRef in a local — never bind `x.AsVector().get()` to
    // a reference (clang warns via lifetimebound).
    class VectorRef
    {
    public:
      std::vector<DataT> const &get() const BIGMATH_LIFETIMEBOUND { return ptr ? *ptr : copy; }
      operator std::vector<DataT> const &() const BIGMATH_LIFETIMEBOUND { return get(); }

    private:
      friend class BigInteger;
      explicit VectorRef(std::vector<DataT> const *p) : ptr(p) {}
      explicit VectorRef(std::span<const DataT> limbs) : copy(limbs.begin(), limbs.end()), ptr(nullptr) {}
      std::vector<DataT> copy;
      std::vector<DataT> const *ptr;
    };

    // Data
  private:
    // Exactly one representation is live:
    //   inline — heap is empty; inl[0..inlSize) holds the limbs,
    //            1 <= inlSize <= InlineLimbs;
    //   heap   — heap holds the limbs, heap.size() > InlineLimbs; inl is {0}
    //            with inlSize 1, so a moved-from heap vector reads as zero.
    // Every constructor and Trim() keep the value inline exactly when it has
    // at most InlineLimbs limbs, so copying a small value never allocates.
    std::vector<DataT> heap;
    DataT inl[InlineLimbs];
    uint32_t inlSize;
    // True if the number is negative
    bool isNegative;

    void SetInlineZero()
    {
      inl[0] = 0;
      inlSize = 1;
    }

    // Store limbs[0..n), choosing the representation by length. Does not trim.
    void AssignLimbs(const DataT *limbs, SizeT n)
    {
      if (n == 0)
      {
        heap.clear();
        SetInlineZero();
      }
      else if (n <= InlineLimbs)
      {
        std::copy(limbs, limbs + n, inl);
        inlSize = (uint32_t)n;
        if (!heap.empty())
          std::vector<DataT>().swap(heap);
      }
      else
      {
        heap.assign(limbs, limbs + n);
        SetInlineZero();
      }
    }

    // Take ownership of v, or copy it inline when short. Does not trim.
    void AdoptVector(std::vector<DataT> &&v)
    {
      if (v.size() <= InlineLimbs)
      {
        AssignLimbs(v.data(), (SizeT)v.size());
        std::vector<DataT>().swap(v);
      }
      else
      {
        heap = std::move(v);
        SetInlineZero();
      }
    }

    // Trim leading zero limbs (keeping at least one), then move the value
    // inline if it now fits.
    SizeT Normalize()
    {
      if (heap.empty())
      {
        while (inlSize > 1 && inl[inlSize - 1] == 0)
          --inlSize;
        return inlSize;
      }
      TrimZerosToOne(heap);
      if (heap.size() <= InlineLimbs)
        AssignLimbs(heap.data(), (SizeT)heap.size());
      return size();
    }

    void ClearNegativeZero()
    {
      if (isNegative && Zero())
        isNegative = false;
    }

    // Set up the inline part for a value of n trimmed limbs: the limbs
    // themselves when they fit, zero otherwise (the caller then fills heap).
    void InitFromTrimmed(const DataT *limbs, SizeT n, bool negative)
    {
      if (n == 0)
      {
        SetInlineZero();
        isNegative = false;
        return;
      }
      if (n <= InlineLimbs)
      {
        std::copy(limbs, limbs + n, inl);
        inlSize = (uint32_t)n;
        isNegative = negative && !(n == 1 && limbs[0] == 0);
      }
      else
      {
        SetInlineZero();
        isNegative = negative;
      }
    }

    // Direct constructors for FromLimbs: limbs already trimmed to n >= 1.
    struct InlineTag {};
    struct HeapTag {};
    BigInteger(InlineTag, const DataT *limbs, SizeT n, bool negative)
        : inlSize((uint32_t)n), isNegative(negative && !(n == 1 && limbs[0] == 0))
    {
      std::copy(limbs, limbs + n, inl);
    }
    BigInteger(HeapTag, const DataT *limbs, SizeT n, bool negative)
        : heap(limbs, limbs + n), inlSize(1), isNegative(negative)
    {
      inl[0] = 0;
    }

    // Constructor, destructor, and assignment operator
  public:
    // A zero with `size` limbs (at least one). The limbs are not trimmed.
    explicit BigInteger(SizeT size = 0, bool negative = false) : isNegative(negative)
    {
      SetInlineZero();
      if (size > InlineLimbs)
        heap.assign(size, 0);
      else if (size > 1)
      {
        std::fill(inl, inl + size, 0);
        inlSize = (uint32_t)size;
      }
      ClearNegativeZero();
    }

    BigInteger(std::vector<DataT> const &aInt, bool negative)
    {
      SizeT n = (SizeT)aInt.size();
      while (n > 1 && aInt[n - 1] == 0)
        --n;
      InitFromTrimmed(aInt.data(), n, negative);
      if (n > InlineLimbs)
        heap.assign(aInt.data(), aInt.data() + n);
    }

    // Zero-copy adoption of a raw limb vector (little-endian limb order).
    // Precondition: every limb must be canonical for the current base —
    // < 2^32 when built with BIGMATH_LIMB_64=0. No validation is performed;
    // non-canonical limbs silently corrupt downstream arithmetic. A value of
    // at most InlineLimbs limbs is copied inline instead; aInt then keeps its
    // buffer (the caller's temporary frees it).
    BigInteger(std::vector<DataT> &&aInt, bool negative)
    {
      SizeT n = (SizeT)aInt.size();
      while (n > 1 && aInt[n - 1] == 0)
        --n;
      InitFromTrimmed(aInt.data(), n, negative);
      if (n > InlineLimbs)
      {
        aInt.resize(n);
        heap = std::move(aInt);
      }
    }

    // Filled with specified data
    BigInteger(SizeT size, bool negative, DataT fill) : isNegative(negative)
    {
      SetInlineZero();
      std::vector<DataT> limbs(size);
      SetBit(limbs, 0, size - 1, fill);
      TrimZerosToOne(limbs);
      AdoptVector(std::move(limbs));
      ClearNegativeZero();
    }

    // Build from little-endian limbs (trimmed; same canonical-limb
    // precondition as the vector constructors). Allocates only when more
    // than InlineLimbs limbs remain after trimming.
    static BigInteger FromLimbs(std::span<const DataT> limbs, bool negative)
    {
      SizeT n = (SizeT)limbs.size();
      while (n > 1 && limbs[n - 1] == 0)
        --n;
      if (n == 0)
        return BigInteger();
      if (n <= InlineLimbs)
        return BigInteger(InlineTag{}, limbs.data(), n, negative);
      return BigInteger(HeapTag{}, limbs.data(), n, negative);
    }

    BigInteger(BigInteger const &) = default;
    BigInteger &operator=(BigInteger const &) = default;
    ~BigInteger() = default;

    // Moves leave the source as zero.
    BigInteger(BigInteger &&o) noexcept
        : heap(std::move(o.heap)), inlSize(o.inlSize), isNegative(o.isNegative)
    {
      std::copy(o.inl, o.inl + InlineLimbs, inl);
      o.heap.clear();
      o.SetInlineZero();
      o.isNegative = false;
    }

    BigInteger &operator=(BigInteger &&o) noexcept
    {
      if (this != &o)
      {
        heap = std::move(o.heap);
        std::copy(o.inl, o.inl + InlineLimbs, inl);
        inlSize = o.inlSize;
        isNegative = o.isNegative;
        o.heap.clear();
        o.SetInlineZero();
        o.isNegative = false;
      }
      return *this;
    }

    // Accessors
  public:
    // The limbs, little-endian. No allocation; valid until *this changes.
    std::span<const DataT> Limbs() const
    {
      if (heap.empty())
        return std::span<const DataT>(inl, inlSize);
      return std::span<const DataT>(heap.data(), heap.size());
    }

    // The limbs as a `std::vector const &` for vector-based algorithms;
    // copies only values short enough to be stored inline.
    VectorRef AsVector() const
    {
      if (heap.empty())
        return VectorRef(Limbs());
      return VectorRef(&heap);
    }

    // A copy of the limbs. (Returned by value since small values are not
    // stored in a vector; use Limbs() or AsVector() to avoid the copy.)
    std::vector<DataT> GetInteger() const
    {
      auto l = Limbs();
      return std::vector<DataT>(l.begin(), l.end());
    }

    std::vector<DataT> ReleaseInteger()
    {
      std::vector<DataT> r = heap.empty() ? GetInteger() : std::move(heap);
      heap.clear();
      SetInlineZero();
      isNegative = false;
      return r;
    }

    // True when the value is stored inside the object (no heap block).
    bool IsInline() const
    {
      return heap.empty();
    }

    // Byte order for magnitude serialization. Named enum instead of a bool so
    // call sites read unambiguously next to FromByteArray's `negative` flag.
    enum class ByteOrder
    {
      BigEndian,
      LittleEndian
    };

    std::vector<uint8_t> ToByteArray(ByteOrder order = ByteOrder::BigEndian) const
    {
      if (Zero())
      {
        return {};
      }

      constexpr SizeT limbBytes = LimbBits / 8;
      std::vector<uint8_t> bytes;
      const auto limbs = Limbs();
      bytes.reserve(limbs.size() * limbBytes);

      for (DataT val : limbs)
      {
        for (SizeT b = 0; b < limbBytes; ++b)
        {
          bytes.push_back(static_cast<uint8_t>((val >> (b * 8)) & 0xFF));
        }
      }

      while (!bytes.empty() && bytes.back() == 0)
      {
        bytes.pop_back();
      }

      if (order == ByteOrder::BigEndian)
      {
        std::reverse(bytes.begin(), bytes.end());
      }

      return bytes;
    }

    static BigInteger FromByteArray(std::span<const uint8_t> bytes, bool negative, ByteOrder order = ByteOrder::BigEndian)
    {
      if (bytes.empty())
      {
        return BigInteger();
      }

      constexpr SizeT limbBytes = LimbBits / 8;
      const size_t numLimbs = (bytes.size() + limbBytes - 1) / limbBytes;
      std::vector<DataT> limbs(numLimbs, 0);

      auto getByte = [&](size_t idx) -> uint8_t {
        if (order == ByteOrder::BigEndian)
        {
          return bytes[bytes.size() - 1 - idx];
        }
        else
        {
          return bytes[idx];
        }
      };

      for (size_t i = 0; i < numLimbs; ++i)
      {
        DataT limbVal = 0;
        for (size_t b = 0; b < limbBytes; ++b)
        {
          size_t byteIdx = i * limbBytes + b;
          if (byteIdx < bytes.size())
          {
            limbVal |= (static_cast<DataT>(getByte(byteIdx)) << (b * 8));
          }
        }
        limbs[i] = limbVal;
      }

      return BigInteger(std::move(limbs), negative);
    }

    DataT operator[](const SizeT i) const
    {
      return heap.empty() ? inl[i] : heap[i];
    }

    // Properties
    SizeT size() const
    {
      return heap.empty() ? (SizeT)inlSize : (SizeT)heap.size();
    }

    static BaseT Base()
    {
      // Resolves to Base2_64 by default (BIGMATH_LIMB_64=1); Base2_32 under -DBIGMATH_LIMB_64=0.
      return CurrentBase;
    }

    bool IsNegative() const
    {
      return isNegative;
    }

    bool Zero() const
    {
      return IsZero(Limbs());
    }

  public:
    // Trims Leading Zeros
    SizeT Trim()
    {
      return Normalize();
    }

    BigInteger &SetSign(bool sign)
    {
      this->isNegative = sign && !Zero();
      return *this;
    }

    // Negation, returns a negated copy; *this is unchanged
    BigInteger operator-() const
    {
      BigInteger r(*this);
      if (!r.Zero())
        r.isNegative = !r.isNegative;
      return r;
    }

  public:
    // Compares this with `with'
    // Returns
    // 0 if equal
    // 1 if this>with
    // -1 if this<with
    Int CompareTo(BigInteger const &with) const
    {
      if (!isNegative && with.isNegative)
        return 1;
      if (isNegative && !with.isNegative)
        return -1;
      const Int magnitude = Compare(Limbs(), with.Limbs());
      return isNegative ? -magnitude : magnitude;
    }
  };
}

#endif
