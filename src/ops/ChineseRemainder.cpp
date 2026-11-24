#include "biginteger/ops/ChineseRemainder.h"

#include <limits>
#include <stdexcept>

#include "biginteger/common/Builder.h"
#include "biginteger/ops/Operations.h"

namespace BigMath
{
  namespace
  {
    DataT WordRemainder(std::span<const DataT> magnitude, DataT modulus)
    {
      if (magnitude.size() == 1)
        return magnitude.front() % modulus;
      return ClassicDivision::RemainderLimbs(magnitude, modulus, CurrentBase);
    }

    DataT InvertWord(DataT value, DataT modulus)
    {
      const auto inverse = ExtendedGcd(BigIntegerBuilder::From(value), BigIntegerBuilder::From(modulus));
      if (inverse.gcd.size() != 1 || inverse.gcd[0] != 1)
        throw std::invalid_argument("Chinese remaindering requires pairwise-coprime moduli");
      const DataT residue = WordRemainder(inverse.x.Limbs(), modulus);
      return inverse.x.IsNegative() && residue != 0 ? modulus - residue : residue;
    }
  }

  ChineseRemainderPlan::ChineseRemainderPlan(std::span<const DataT> input)
      : moduli(input.begin(), input.end()), cofactors(input.size()), weights(input.size())
  {
    if (moduli.empty())
      throw std::invalid_argument("Chinese remaindering requires at least one modulus");

    std::size_t first = 0;
    DataT product = 1;
    for (std::size_t i = 0; i < moduli.size(); ++i)
    {
      const DataT modulus = moduli[i];
      if (modulus < 2)
        throw std::invalid_argument("Chinese remaindering requires moduli of at least two");
      const ULong128 combined = ULong128{product} * modulus;
      if (combined > std::numeric_limits<DataT>::max())
      {
        groups.push_back({first, i, product});
        first = i;
        product = modulus;
      }
      else
        product = static_cast<DataT>(combined);
    }
    groups.push_back({first, moduli.size(), product});

    // A subtree with k leaves has 2*k - 1 nodes. Store each left subtree
    // immediately after its parent, followed by the right subtree.
    products.resize(2 * groups.size() - 1);
    BuildProducts(0, 0, groups.size());
    BuildWeights(0, 0, groups.size(), BigIntegerBuilder::From(ULong{1}));
  }

  ChineseRemainderPlan &ChineseRemainderPlan::operator=(ChineseRemainderPlan other) noexcept
  {
    moduli.swap(other.moduli);
    groups.swap(other.groups);
    products.swap(other.products);
    cofactors.swap(other.cofactors);
    weights.swap(other.weights);
    return *this;
  }

  void ChineseRemainderPlan::BuildProducts(std::size_t node, std::size_t first, std::size_t count)
  {
    if (count == 1)
    {
      products[node] = BigIntegerBuilder::From(groups[first].product);
      return;
    }
    const auto leftCount = count / 2;
    const auto right = node + 2 * leftCount;
    BuildProducts(node + 1, first, leftCount);
    BuildProducts(right, first + leftCount, count - leftCount);
    products[node] = products[node + 1] * products[right];
  }

  void ChineseRemainderPlan::BuildWeights(std::size_t node, std::size_t first, std::size_t count,
                                         BigInteger const &outside)
  {
    if (count == 1)
    {
      const auto &group = groups[first];
      // outside is the product of all moduli outside this group, reduced
      // modulo the group product. Complete each cofactor inside the word.
      const DataT outer = WordRemainder(outside.Limbs(), group.product);
      for (auto i = group.first; i < group.end; ++i)
      {
        const DataT modulus = moduli[i];
        cofactors[i] = group.product / modulus;
        const DataT cofactor = (ULong128{outer} * cofactors[i]) % modulus;
        weights[i] = InvertWord(cofactor, modulus);
      }
      return;
    }

    const auto left = node + 1;
    const auto leftCount = count / 2;
    const auto right = node + 2 * leftCount;
    BuildWeights(left, first, leftCount, (outside * products[right]) % products[left]);
    BuildWeights(right, first + leftCount, count - leftCount, (outside * products[left]) % products[right]);
  }

  void ChineseRemainderPlan::ReduceAt(std::size_t node, std::size_t first, std::size_t count,
                                      std::span<const DataT> value, std::span<DataT> residues) const
  {
    if (count == 1)
    {
      const auto &group = groups[first];
      const DataT reduced = WordRemainder(value, group.product);
      if (group.end == group.first + 1)
      {
        residues[group.first] = reduced;
        return;
      }
      for (auto i = group.first; i < group.end; ++i)
        residues[i] = reduced % moduli[i];
      return;
    }

    const auto product = products[node].Limbs();
    std::vector<DataT> remainder;
    if (Compare(value, product) >= 0)
    {
      remainder = detail::DivideAndRemainderLimbs(value, product, CurrentBase, true).second;
      value = remainder;
    }
    const auto leftCount = count / 2;
    ReduceAt(node + 1, first, leftCount, value, residues);
    ReduceAt(node + 2 * leftCount, first + leftCount, count - leftCount, value, residues);
  }

  std::vector<DataT> ChineseRemainderPlan::Reduce(BigInteger const &value) const
  {
    std::vector<DataT> residues(moduli.size());
    if (residues.empty() || value.Zero())
      return residues;
    ReduceAt(0, 0, groups.size(), value.Limbs(), residues);
    if (value.IsNegative())
      for (std::size_t i = 0; i < residues.size(); ++i)
        if (residues[i] != 0)
          residues[i] = moduli[i] - residues[i];
    return residues;
  }

  BigInteger ChineseRemainderPlan::ReconstructAt(std::size_t node, std::size_t first, std::size_t count,
                                                 std::span<const DataT> residues) const
  {
    if (count == 1)
    {
      const auto &group = groups[first];
      // Each summand is below the group's word-size product and a group
      // contains at most 63 moduli, so the accumulator fits in 128 bits.
      ULong128 sum = 0;
      for (auto i = group.first; i < group.end; ++i)
      {
        const DataT weighted = (ULong128{residues[i]} * weights[i]) % moduli[i];
        sum += ULong128{weighted} * cofactors[i];
      }
      const auto reduced = sum < group.product ? sum : sum % group.product;
      return BigIntegerBuilder::From(static_cast<ULong>(reduced));
    }

    const auto left = node + 1;
    const auto leftCount = count / 2;
    const auto right = node + 2 * leftCount;
    const auto a = ReconstructAt(left, first, leftCount, residues);
    const auto b = ReconstructAt(right, first + leftCount, count - leftCount, residues);
    return a * products[right] + b * products[left];
  }

  BigInteger ChineseRemainderPlan::Reconstruct(std::span<const DataT> residues) const
  {
    if (residues.empty() || residues.size() != moduli.size())
      throw std::invalid_argument("Chinese remaindering requires one residue per modulus");
    auto result = ReconstructAt(0, 0, groups.size(), residues);
    if (result >= Modulus())
      return result % Modulus();
    return result;
  }
}
