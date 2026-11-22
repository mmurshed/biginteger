#pragma once

// Included within BigMath::NttCrt. Weighted cyclic products evaluate the
// coefficient blocks at powers of g^H. A short cyclic product supplies the
// final node at 1. Recover the top block first, then interpolate the others.
// Replacing the smallest standard CRT prime widens the centered coefficient
// range and supplies up to 54 distinct nodes at H = 2^25.
#if BIGMATH_NTT_TRUNC && BIGMATH_NTT_MFA && BIGMATH_NTT_MFA_FUSE
namespace Wide
{
  constexpr UInt W1 = P1, W2 = P3, W3 = 2113929217u;
  using FW1 = ModField<W1>;
  using FW2 = ModField<W2>;
  using FW3 = ModField<W3>;
  constexpr UInt GW1 = G1, GW2 = G3, GW3 = 5;
  constexpr ULong128 WideModulus = (ULong128)W1 * W2 * W3;
  constexpr ULong WideMaxBiasedTerms = (ULong)(((WideModulus - 1) / 2) >> 62);
  constexpr Int WideMaxResidues = 53;
  constexpr Int WideResidueLength = 1 << 25;
  static_assert(W1 < W3 && W2 < W3 && W3 < (1u << 31) && 2 * (ULong)W3 < (1ULL << 32));
  static_assert((ULong128)4 * (W3 - 1) * (W3 - 1) + W3 < ((ULong128)1 << 64));
  static_assert((W1 - 1) % WideResidueLength == 0 &&
                (W2 - 1) % WideResidueLength == 0 &&
                (W3 - 1) % WideResidueLength == 0);
  static_assert(WideMaxResidues < (W1 - 1) / WideResidueLength &&
                WideMaxResidues < (W2 - 1) / WideResidueLength &&
                WideMaxResidues < (W3 - 1) / WideResidueLength);
  static_assert(((ULong128)0xFFFFFFFFu + BiasedOperandView::Offset<FW3>()) * (W3 - 1) <
                ((ULong128)1 << 64));
  static_assert(WideMaxBiasedTerms == 836075521 && WideMaxBiasedTerms < (1ULL << 30));

  template <typename F>
  struct PowerTable
  {
    static constexpr Int Stride = 256;
    std::array<UInt, Stride> lo;
    std::vector<UInt> hi;
#if BIGMATH_NEON
    std::array<UInt, Stride> loShoup;
    std::vector<UInt> hiShoup;
#endif

    PowerTable(UInt w, Int H) : hi((H + Stride - 1) / Stride)
#if BIGMATH_NEON
                              , hiShoup(hi.size())
#endif
    {
      UInt power = 1;
      for (Int i = 0; i < Stride; ++i)
      {
        lo[i] = power;
#if BIGMATH_NEON
        loShoup[i] = ShoupCompanion<F>(power);
#endif
        power = F::Mul(power, w);
      }
      const UInt step = power;
      power = 1;
      for (size_t i = 0; i < hi.size(); ++i)
      {
        hi[i] = power;
#if BIGMATH_NEON
        hiShoup[i] = ShoupCompanion<F>(power);
#endif
        power = F::Mul(power, step);
      }
    }

    UInt At(SizeT i) const { return F::Mul(hi[i / Stride], lo[i % Stride]); }

#if BIGMATH_NEON
    // i is a multiple of four, so the lanes share the high factor.
    uint32x4_t Apply(SizeT i, uint32x4_t values) const
    {
      const SizeT low = i % Stride, high = i / Stride;
      values = NeonShoupMul<F>(vld1q_u32(lo.data() + low), vld1q_u32(loShoup.data() + low), values);
      return NeonShoupMul<F>(vdupq_n_u32(hi[high]), vdupq_n_u32(hiShoup[high]), values);
    }
#endif
  };

  template <typename F>
  void ApplyWeight(const PowerTable<F> &weight, SizeT first, Int count, UInt *values)
  {
    Int x = 0;
#if BIGMATH_NEON
    for (; x + 4 <= count; x += 4)
      vst1q_u32(values + x, weight.Apply(first + x, vld1q_u32(values + x)));
#endif
    for (; x < count; ++x)
      values[x] = F::Mul(values[x], weight.At(first + x));
  }

#if BIGMATH_NEON
  template <typename F>
  uint32x4_t LoadBiased(const BiasedOperandView &src, SizeT i)
  {
    uint32x4_t raw;
    if (src.base64)
      raw = vreinterpretq_u32_u64(vld1q_u64(src.v + i / 2));
    else
      raw = vmovn_high_u64(vmovn_u64(vld1q_u64(src.v + i)), vld1q_u64(src.v + i + 2));
    // All wide primes exceed 2^30 and three times the smallest exceeds 2^32.
    static_assert((ULong)3 * F::Prime > 0xFFFFFFFFULL);
    const auto p = vdupq_n_u32(F::Prime);
    raw = vminq_u32(raw, vsubq_u32(raw, p));
    raw = vminq_u32(raw, vsubq_u32(raw, p));
    return NeonAdd<F>(raw, vdupq_n_u32(BiasedOperandView::Offset<F>()));
  }

  // Four products and an already reduced partial sum fit in a 64-bit lane.
  inline void AccumulateFour(ULong *sums, uint32x4_t values, UInt factor)
  {
    const auto w = vdup_n_u32(factor);
    vst1q_u64(sums, vmlal_u32(vld1q_u64(sums), vget_low_u32(values), w));
    vst1q_u64(sums + 2, vmlal_u32(vld1q_u64(sums + 2), vget_high_u32(values), w));
  }
#endif

  template <typename F, typename Operand>
  void FoldTwist(const Operand &src, Int H, UInt s, const PowerTable<F> &w, UInt *plane)
  {
    std::vector<UInt> powers((src.coeffCount + H - 1) / H);
    UInt power = 1;
    for (UInt &p : powers)
    {
      p = power;
      power = F::Mul(power, s);
    }

    constexpr Int Tile = 512;
    ParallelFor(H / Tile + (H % Tile != 0), [&](Int first, Int end) {
      for (Int tile = first; tile < end; ++tile)
      {
        const SizeT i = (SizeT)tile * Tile;
        const Int n = (Int)std::min<SizeT>(Tile, H - i);
        ULong sums[Tile]{};
        for (ULong u = i, turn = 0; u < src.coeffCount; u += H, ++turn)
        {
          const Int count = (Int)std::min<ULong>(n, src.coeffCount - u);
          Int x = 0;
#if BIGMATH_NEON
          if constexpr (std::is_same_v<Operand, BiasedOperandView>)
            for (; x + 4 <= count; x += 4)
              AccumulateFour(sums + x, LoadBiased<F>(src, (SizeT)u + x), powers[turn]);
#endif
          for (; x < count; ++x)
            sums[x] += (ULong)src.template Coeff<F>((SizeT)u + x) * powers[turn];
          if (turn % 4 == 3)
            for (Int t = 0; t < n; ++t) sums[t] %= F::Prime;
        }
        for (Int x = 0; x < n; ++x) plane[i + x] = (UInt)(sums[x] % F::Prime);
        ApplyWeight(w, i, n, plane + i);
      }
    });
  }

  template <typename F>
  std::vector<UInt> LagrangeBasis(const std::vector<UInt> &nodes, Int count, Int j)
  {
    std::vector<UInt> poly(count);
    poly[0] = 1;
    UInt scale = 1;
    Int degree = 0;
    for (Int l = 0; l < count; ++l)
    {
      if (l == j) continue;
      for (Int d = ++degree; d > 0; --d)
        poly[d] = F::Sub(poly[d - 1], F::Mul(poly[d], nodes[l]));
      poly[0] = F::Sub(0, F::Mul(poly[0], nodes[l]));
      scale = F::Mul(scale, F::Sub(nodes[j], nodes[l]));
    }
    scale = F::Inv(scale);
    for (UInt &c : poly) c = F::Mul(c, scale);
    return poly;
  }

  // Reduce one group of at most four terms per lane. Keeping the widening
  // accumulators in registers avoids streaming 64-bit sums for each row.
  template <typename F, Int Terms>
  void InterpolationGroup(const UInt *const *rows, const UInt *factors, Int n, UInt *sums)
  {
    Int x = 0;
#if BIGMATH_NEON
    for (; x + 4 <= n; x += 4)
    {
      const auto prior = vld1q_u32(sums + x);
      auto low = vmovl_u32(vget_low_u32(prior));
      auto high = vmovl_u32(vget_high_u32(prior));
      for (Int j = 0; j < Terms; ++j)
      {
        const auto values = vld1q_u32(rows[j] + x);
        low = vmlal_n_u32(low, vget_low_u32(values), factors[j]);
        high = vmlal_n_u32(high, vget_high_u32(values), factors[j]);
      }
      sums[x] = (UInt)(vgetq_lane_u64(low, 0) % F::Prime);
      sums[x + 1] = (UInt)(vgetq_lane_u64(low, 1) % F::Prime);
      sums[x + 2] = (UInt)(vgetq_lane_u64(high, 0) % F::Prime);
      sums[x + 3] = (UInt)(vgetq_lane_u64(high, 1) % F::Prime);
    }
#endif
    for (; x < n; ++x)
    {
      ULong sum = sums[x];
      for (Int j = 0; j < Terms; ++j) sum += (ULong)rows[j][x] * factors[j];
      sums[x] = (UInt)(sum % F::Prime);
    }
  }

  template <typename F>
  void Interpolate(const std::vector<UInt> &nodes, Int H, Int k, Int L, Int tail, UInt *out)
  {
    const UInt topTail = LagrangeBasis<F>(nodes, k + 1, k)[k];
    std::vector<UInt> top(k), blockTop(k);
    std::vector<std::vector<UInt>> basis(k, std::vector<UInt>(k));
    for (Int j = 0; j < k; ++j)
    {
      top[j] = LagrangeBasis<F>(nodes, k + 1, j)[k];
      const auto column = LagrangeBasis<F>(nodes, k, j);
      const UInt power = F::Power(nodes[j], k);
      for (Int m = 0; m < k; ++m)
      {
        basis[m][j] = column[m];
        blockTop[m] = F::Sub(blockTop[m], F::Mul(column[m], power));
      }
    }

    constexpr Int Tile = 128;
    UInt *const topBlock = out + (SizeT)k * H;
    ParallelFor((tail + Tile - 1) / Tile, [&](Int first, Int end) {
      for (Int tile = first; tile < end; ++tile)
      {
        const Int i = tile * Tile, n = std::min(Tile, tail - i);
        ULong sums[Tile];
        for (Int x = 0; x < n; ++x) sums[x] = F::Mul(topTail, topBlock[i + x]);
        for (Int j = 0; j < k; ++j)
        {
          const UInt *const row = out + (SizeT)j * H;
          if (L == H)
          {
            Int x = 0;
#if BIGMATH_NEON
            for (; x + 4 <= n; x += 4)
              AccumulateFour(sums + x, vld1q_u32(row + i + x), top[j]);
#endif
            for (; x < n; ++x) sums[x] += (ULong)row[i + x] * top[j];
          }
          else
          {
            ULong folded[Tile]{};
            for (Int u = i; u < H; u += L)
              for (Int x = 0; x < n; ++x) folded[x] += row[u + x];
            for (Int x = 0; x < n; ++x)
              sums[x] += (ULong)top[j] * (folded[x] % F::Prime);
          }
          if (j % 4 == 3)
            for (Int x = 0; x < n; ++x) sums[x] %= F::Prime;
        }
        for (Int x = 0; x < n; ++x) topBlock[i + x] = (UInt)(sums[x] % F::Prime);
      }
    });

    ParallelFor(H / Tile + (H % Tile != 0), [&](Int first, Int end) {
      static thread_local std::vector<UInt> buffer;
      buffer.resize((SizeT)k * Tile);
      for (Int tile = first; tile < end; ++tile)
      {
        const Int i = tile * Tile, n = std::min(Tile, H - i);
        const Int topCount = std::clamp(tail - i, 0, n);
        std::array<const UInt *, WideMaxResidues> rows;
        for (Int j = 0; j < k; ++j) rows[j] = out + (SizeT)j * H + i;
        for (Int m = 0; m < k; ++m)
        {
          UInt *const sums = buffer.data() + (SizeT)m * Tile;
          for (Int x = 0; x < topCount; ++x) sums[x] = F::Mul(blockTop[m], topBlock[i + x]);
          std::fill(sums + topCount, sums + n, 0);
          const UInt *const factors = basis[m].data();
          Int j = 0;
          for (; j + 4 <= k; j += 4)
            InterpolationGroup<F, 4>(rows.data() + j, factors + j, n, sums);
          switch (k - j)
          {
          case 1: InterpolationGroup<F, 1>(rows.data() + j, factors + j, n, sums); break;
          case 2: InterpolationGroup<F, 2>(rows.data() + j, factors + j, n, sums); break;
          case 3: InterpolationGroup<F, 3>(rows.data() + j, factors + j, n, sums); break;
          }
        }
        // Commit only after every block has read the old values of this tile.
        for (Int m = 0; m < k; ++m)
          std::copy_n(buffer.data() + (SizeT)m * Tile, n, out + (SizeT)m * H + i);
      }
    });
  }

  template <typename F, UInt G, typename Operand>
  void WideProductPrime(const Operand &a, const Operand &b, bool square,
                        Int H, Int k, Int L, Int tail, UInt *out,
                        UInt *foldA, UInt *foldB, UInt *sa, UInt *sb)
  {
    const ResidueView left{foldA}, right{square ? foldA : foldB};
    FoldOperand<F>(a, L, foldA);
    if (!square) FoldOperand<F>(b, L, foldB);
    CyclicProductViews<F, G>(left, right, L, out, sa, sb);
    std::copy_n(out, tail, out + (SizeT)k * H);

    std::vector<UInt> nodes(k + 1, 1);
    const UInt step = F::Power(G, H), inverse = F::Inv(G);
    UInt w = 1, invW = 1, s = 1;
    for (Int j = 0; j < k; ++j)
    {
      w = F::Mul(w, G);
      invW = F::Mul(invW, inverse);
      nodes[j] = s = F::Mul(s, step);
      const PowerTable<F> weight(w, H);
      FoldTwist(a, H, s, weight, foldA);
      if (!square) FoldTwist(b, H, s, weight, foldB);
      UInt *const row = out + (SizeT)j * H;
      CyclicProductViews<F, G>(left, right, H, row, sa, sb);
      // Unweight immediately, so only one inverse table is live.
      const PowerTable<F> unweight(invW, H);
      ParallelFor(H / PowerTable<F>::Stride + (H % PowerTable<F>::Stride != 0), [&](Int first, Int end) {
        const Int from = first * PowerTable<F>::Stride, to = std::min(H, end * PowerTable<F>::Stride);
        ApplyWeight(unweight, (SizeT)from, to - from, row + from);
      });
    }
    Interpolate<F>(nodes, H, k, L, tail, out);
  }

  struct WideCrt
  {
    static constexpr ULong128 Modulus = WideModulus;
    struct Inverses { UInt w1w2, w1w3, w2w3; };

    static const Inverses &GetInverses()
    {
      static const Inverses inv{FW2::Inv(W1 % W2), FW3::Inv(W1 % W3), FW3::Inv(W2 % W3)};
      return inv;
    }

    static ULong128 Value(UInt r1, UInt r2, UInt r3, const Inverses &inv)
    {
      const UInt u2 = FW2::Mul(FW2::Sub(r2, r1 % W2), inv.w1w2);
      const UInt u3 = FW3::Mul(FW3::Sub(FW3::Mul(FW3::Sub(r3, r1 % W3), inv.w1w3), u2), inv.w2w3);
      return r1 + (ULong)W1 * u2 + (ULong128)W1 * W2 * u3;
    }
  };

  inline std::vector<DataT> MultiplyWide(const std::vector<DataT> &a, const std::vector<DataT> &b,
                                         BaseT base, Int H = WideResidueLength)
  {
    const PackedOperandView va(a, base), vb(b, base);
    const ULong shorter = std::min(va.coeffCount, vb.coeffCount);
    const ULong count = va.coeffCount + vb.coeffCount - (shorter != 0);
    if (H < 64 || H > WideResidueLength || !std::has_single_bit((UInt)H) ||
        count <= (ULong)H || count > (ULong)(WideMaxResidues + 1) * H ||
        shorter == 0 || shorter > WideMaxBiasedTerms)
      throw std::invalid_argument("NttCrt::Wide::MultiplyWide: product does not fit the wide residues");

    // W1/W2 can double H at full size. This halves their interpolation
    // dimensions and keeps the same node limit (W2 has 27 nodes at 2^26).
    const Int longH = H == WideResidueLength && count > 2 * (ULong)H && MfaFusedWindow(2 * H)
                          ? 2 * H : H;
    const bool square = &a == &b;
    auto &workspace = GetProductWorkspace();
    for (auto *plane : {&workspace.fa1, &workspace.fa2, &workspace.fa3})
    {
      if (plane->capacity() < count) std::vector<UInt>().swap(*plane);
      plane->resize(count);
    }
    const BiasedOperandView xa(va), xb(vb);
    {
      const bool single = ParallelNumThreads() == 1;
      auto &sa = single ? workspace.fb1 : workspace.transpose[0];
      auto &sb = single ? workspace.fb2 : workspace.transpose[1];
      sa.resize(longH);
      sb.resize(longH);
      std::vector<UInt> foldA(longH), foldB(square ? 0 : longH);
      auto prime = [&]<typename F, UInt G>(Int length, UInt *out) {
        const Int k = (Int)((count - 1) / length);
        const Int tail = (Int)(count - (ULong)k * length);
        const Int L = std::max(64, (Int)std::bit_ceil((UInt)tail));
        WideProductPrime<F, G>(xa, xb, square, length, k, L, tail, out,
                               foldA.data(), foldB.data(), sa.data(), sb.data());
      };
      prime.template operator()<FW1, GW1>(longH, workspace.fa1.data());
      prime.template operator()<FW2, GW2>(longH, workspace.fa2.data());
      prime.template operator()<FW3, GW3>(H, workspace.fa3.data());
      std::vector<UInt>().swap(sa);
      std::vector<UInt>().swap(sb);
    }
    return FinalizeBiasedProduct<WideCrt>(workspace.fa1, workspace.fa2, workspace.fa3,
                                          va, vb, base, (SizeT)(a.size() + b.size() + 2));
  }
}
#endif
