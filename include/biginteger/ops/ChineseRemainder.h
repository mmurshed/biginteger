#ifndef BIGINTEGER_CHINESE_REMAINDER
#define BIGINTEGER_CHINESE_REMAINDER

#include <span>
#include <vector>

#include "../BigInteger.h"

namespace BigMath
{
  // Precomputed product tree for a nonempty list of pairwise-coprime word moduli.
  class ChineseRemainderPlan
  {
    struct WordGroup
    {
      std::size_t first;
      std::size_t end;
      DataT product;
    };

    std::vector<DataT> moduli;
    std::vector<WordGroup> groups;
    std::vector<BigInteger> products;
    std::vector<DataT> cofactors;
    std::vector<DataT> weights;

    void BuildProducts(std::size_t node, std::size_t first, std::size_t count);
    void BuildWeights(std::size_t node, std::size_t first, std::size_t count,
                      BigInteger const &outside);
    void ReduceAt(std::size_t node, std::size_t first, std::size_t count,
                  std::span<const DataT> value, std::span<DataT> residues) const;
    BigInteger ReconstructAt(std::size_t node, std::size_t first, std::size_t count,
                             std::span<const DataT> residues) const;

  public:
    // Throws std::invalid_argument if the list is empty, a modulus is less
    // than two, or any two moduli have a common factor.
    explicit ChineseRemainderPlan(std::span<const DataT> moduli);
    ChineseRemainderPlan(ChineseRemainderPlan const &) = default;
    ChineseRemainderPlan(ChineseRemainderPlan &&) noexcept = default;
    ChineseRemainderPlan &operator=(ChineseRemainderPlan other) noexcept;

    std::span<const DataT> Moduli() const { return moduli; }
    BigInteger const &Modulus() const { return products.at(0); }

    // Return the least nonnegative residue at each modulus, in input order.
    std::vector<DataT> Reduce(BigInteger const &value) const;

    // Return the unique representative in [0, Modulus()). Residues may
    // exceed their modulus. Throws if the number of residues is wrong.
    BigInteger Reconstruct(std::span<const DataT> residues) const;
  };
}

#endif
