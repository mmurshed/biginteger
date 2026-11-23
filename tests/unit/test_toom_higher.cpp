#include "unit_test_framework.h"

#include <array>
#include <limits>
#include <random>
#include <utility>
#include <vector>

#include "biginteger/algorithms/multiplication/ClassicMultiplication.h"
#include "biginteger/algorithms/multiplication/Toom5Multiplication.h"
#include "biginteger/algorithms/multiplication/ToomHigherMultiplication.h"

using namespace BigMath;

namespace
{
  using Split = ToomHigherMultiplication::Split;
  constexpr std::array splits{Split::Toom5, Split::Toom65, Split::Toom85};

  void Check(Split split, const std::vector<DataT> &a, const std::vector<DataT> &b, BaseT base = Base2_64)
  {
    auto expected = ClassicMultiplication::Multiply(a, b, base);
    TrimZerosToOne(expected);
    ASSERT_TRUE(ToomHigherMultiplication::Multiply(a, b, base, split) == expected);
    ASSERT_TRUE(ToomHigherMultiplication::Multiply(b, a, base, split) == expected);
  }
}

#if BIGMATH_LIMB_64

REGISTER_TEST(ToomHigher, DistinctBalancedAndUnequalProducts)
{
  std::mt19937_64 random(0x5685);
  constexpr std::array shapes{
      std::pair{2, 1}, std::pair{5, 5}, std::pair{7, 6}, std::pair{9, 8}, std::pair{20, 19},
      std::pair{80, 18}, std::pair{164, 159}, std::pair{513, 512}, std::pair{700, 600},
      std::pair{1153, 1025}, std::pair{2049, 2041}, std::pair{4097, 4089}};
  for (const auto &[la, lb] : shapes)
  {
    std::vector<DataT> a(la), b(lb);
    for (auto &limb : a)
      limb = random();
    for (auto &limb : b)
      limb = random();
    for (const Split split : splits)
      Check(split, a, b);
  }
}

REGISTER_TEST(ToomHigher, FullCarriesZeroAndNegativeEvaluations)
{
  constexpr DataT full = std::numeric_limits<DataT>::max();
  constexpr size_t width = 29;
  constexpr std::array tallerParts{5, 7, 9};
  constexpr std::array shorterParts{5, 6, 8};
  for (size_t i = 0; i < splits.size(); ++i)
  {
    std::vector<DataT> a(width * tallerParts[i], full), b(width * shorterParts[i], full);
    Check(splits[i], a, b);
    // Constant equal chunks give zero at -1 when the piece count is even.
    std::fill(a.end() - width, a.end(), 0);
    Check(splits[i], a, b);
    // Large odd chunks make the first operand's evaluation negative at -1.
    // The second operand has a positive evaluation at that point.
    for (size_t j = 0; j < a.size(); ++j)
      a[j] = (j / width) % 2 ? full : 1;
    for (size_t j = 0; j < b.size(); ++j)
      b[j] = (j / width) % 2 ? 1 : full;
    Check(splits[i], a, b);
    Check(splits[i], a, a);
    Check(splits[i], {}, a);
    Check(splits[i], std::vector<DataT>(a.size(), 0), a);
  }
}

#endif

REGISTER_TEST(ToomHigher, LegacyBaseAndDirectToom5Entry)
{
  const std::vector<DataT> a(513, 0xffffffffU), b(512, 0x80000001U);
  for (const Split split : splits)
    Check(split, a, b, Base2_32);
  ASSERT_TRUE(Toom5Multiplication::Multiply(a, b, Base2_32) ==
              ClassicMultiplication::Multiply(a, b, Base2_32));
#if BIGMATH_LIMB_64
  ASSERT_TRUE(Toom5Multiplication::Multiply(a, b, Base2_64) ==
              ClassicMultiplication::Multiply(a, b, Base2_64));
#endif
}
