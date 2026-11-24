#include "unit_test_framework.h"

#include <array>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>

#include "biginteger/common/Builder.h"
#include "biginteger/ops/Operations.h"

using namespace BigMath;

REGISTER_TEST(ChineseRemainder, SmallAndWordModuli)
{
  constexpr std::array<DataT, 3> small{3, 5, 7};
  const ChineseRemainderPlan plan(small);
  const auto word = [](ULong value) { return BigIntegerBuilder::From(value); };
  ASSERT_TRUE(plan.Modulus() == word(105));
  ASSERT_TRUE(plan.Reconstruct(std::array<DataT, 3>{2, 3, 2}) == word(23));
  ASSERT_TRUE(plan.Reconstruct(std::array<DataT, 3>{5, 8, 16}) == word(23));
  ASSERT_TRUE(plan.Reduce(-word(82)) == std::vector<DataT>({2, 3, 2}));
  ASSERT_TRUE(plan.Reduce(BigInteger(23)) == std::vector<DataT>({0, 0, 0}));

  constexpr DataT maxWord = std::numeric_limits<DataT>::max();
  const ChineseRemainderPlan singleton(std::array<DataT, 1>{maxWord});
  ASSERT_TRUE(singleton.Reconstruct(std::array<DataT, 1>{maxWord}).Zero());
  ASSERT_TRUE(singleton.Reduce(-word(2)) == std::vector<DataT>({maxWord - 2}));
  const ChineseRemainderPlan wide(std::array<DataT, 2>{maxWord, maxWord - 1});
  const auto input = wide.Modulus() - word(1);
  ASSERT_TRUE(wide.Reduce(input) == std::vector<DataT>({maxWord - 1, maxWord - 2}));
  ASSERT_TRUE(wide.Reconstruct(wide.Reduce(input)) == input);
  auto copy = plan;
  copy = wide;
  ASSERT_TRUE(copy.Reconstruct(copy.Reduce(input)) == input);
}

REGISTER_TEST(ChineseRemainder, ReusedPlan)
{
  std::mt19937_64 rng(0xF3415);
  std::vector<DataT> moduli;
  for (DataT candidate = 2; moduli.size() < 129; ++candidate)
    if (std::all_of(moduli.begin(), moduli.end(), [candidate](DataT m) {
          return std::gcd(candidate, m) == 1;
        }))
      moduli.push_back(candidate);
  const ChineseRemainderPlan plan(moduli);
  ASSERT_TRUE(std::equal(plan.Moduli().begin(), plan.Moduli().end(), moduli.begin()));

  auto product = BigIntegerBuilder::From(ULong{1});
  for (DataT m : moduli)
    product = product * m;
  ASSERT_TRUE(product == plan.Modulus());

  std::vector<uint8_t> bytes(product.ToByteArray().size() * 2 + 1);
  for (unsigned trial = 0; trial < 12; ++trial)
  {
    for (auto &byte : bytes)
      byte = static_cast<uint8_t>(rng());
    const auto input = BigInteger::FromByteArray(bytes, trial % 2 == 0);
    const auto residues = plan.Reduce(input);
    for (std::size_t i = 0; i < moduli.size(); ++i)
    {
      const auto remainder = input % moduli[i];
      const DataT magnitude = remainder[0];
      ASSERT_EQ(residues[i], remainder.IsNegative() ? moduli[i] - magnitude : magnitude);
    }
    auto expected = input % product;
    if (expected.IsNegative())
      expected = expected + product;
    ASSERT_TRUE(plan.Reconstruct(residues) == expected);
  }
}

REGISTER_TEST(ChineseRemainder, InvalidInputs)
{
  for (const auto &moduli : {std::vector<DataT>{}, {1}, {0, 5}, {3, 5, 3},
                             {std::numeric_limits<DataT>::max(), 6}})
  {
    bool threw = false;
    try { (void)ChineseRemainderPlan(moduli); }
    catch (std::invalid_argument const &) { threw = true; }
    ASSERT_TRUE(threw);
  }
  const ChineseRemainderPlan plan(std::array<DataT, 2>{3, 5});
  bool threw = false;
  try { (void)plan.Reconstruct(std::array<DataT, 1>{3}); }
  catch (std::invalid_argument const &) { threw = true; }
  ASSERT_TRUE(threw);
}
