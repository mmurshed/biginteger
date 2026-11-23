#include <algorithm>
#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <span>
#include <string_view>
#include <vector>

#include "biginteger/algorithms/Multiplication.h"
#include "biginteger/algorithms/multiplication/ToomHigherMultiplication.h"

using namespace BigMath;

namespace
{
  using Limbs = std::vector<DataT>;
  using Split = ToomHigherMultiplication::Split;
  using Product = Limbs (*)(const Limbs &, const Limbs &);

  struct Method
  {
    std::string_view name;
    Product multiply;
  };

  template <Split Candidate>
  Limbs Higher(const Limbs &a, const Limbs &b)
  {
    return ToomHigherMultiplication::Multiply(a, b, Base2_64, Candidate);
  }

  constexpr std::array methods{
      Method{"dispatch", [](const Limbs &a, const Limbs &b) { return Multiply(a, b, Base2_64); }},
      Method{"recursive", [](const Limbs &a, const Limbs &b) {
        return KaratsubaMultiplication::Multiply(a, b, Base2_64);
      }},
      Method{"ntt", [](const Limbs &a, const Limbs &b) { return NTTMultiplication::Multiply(a, b, Base2_64); }},
      Method{"toom5", Higher<Split::Toom5>},
      Method{"toom65", Higher<Split::Toom65>},
      Method{"toom85", Higher<Split::Toom85>},
  };

  struct Shape
  {
    std::string_view name;
    SizeT numerator;
    SizeT denominator;
  };

  constexpr std::array shapes{Shape{"balanced", 1, 1}, Shape{"7:6", 7, 6}, Shape{"9:8", 9, 8}};
  constexpr std::array<SizeT, 20> sweep{
      160, 192, 224, 256, 288, 320, 352, 384, 448, 512, 576, 640, 704, 768, 1024, 1536, 2048, 3072, 4096, 8192};
  constexpr std::array<SizeT, 20> holdout{
      173, 207, 239, 271, 303, 337, 367, 415, 479, 545, 607, 673, 735, 895, 1279, 1793, 2559, 3583, 6143, 10241};
  constexpr std::array<SizeT, 29> balancedBand{
      640, 648, 656, 664, 672, 680, 688, 696, 704, 712, 720, 728, 733, 734, 735, 736, 737,
      744, 752, 760, 768, 776, 784, 792, 800, 808, 816, 824, 832};
  constexpr std::array<SizeT, 25> balancedHoldout{
      644, 652, 660, 668, 676, 684, 692, 700, 708, 716, 724, 732, 740,
      748, 756, 764, 772, 780, 788, 796, 804, 812, 820, 828, 836};
  constexpr unsigned rounds = 11;

  Limbs RandomLimbs(SizeT length, std::mt19937_64 &random)
  {
    Limbs value(length);
    for (auto &limb : value)
      limb = random();
    value.back() |= DataT{1} << 63;
    return value;
  }

  double Median(std::array<double, rounds> values)
  {
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
  }

  bool RunCase(Shape shape, SizeT shorter, std::mt19937_64 &random, DataT &checksum)
  {
    const SizeT longer = (shorter * shape.numerator + shape.denominator - 1) / shape.denominator;
    const auto a = RandomLimbs(longer, random), b = RandomLimbs(shorter, random);
    const auto expected = ClassicMultiplication::Multiply(a, b, Base2_64);
    for (const auto &method : methods)
      if (method.multiply(a, b) != expected)
      {
        std::cerr << method.name << " disagrees with schoolbook at " << longer << " x " << shorter << '\n';
        return false;
      }

    const unsigned reps = std::clamp<unsigned>(8192 / longer, 2, 64);
    std::array<std::array<double, rounds>, methods.size()> samples{};
    std::array<size_t, methods.size()> order{};
    for (size_t i = 0; i < order.size(); ++i)
      order[i] = i;

    for (unsigned round = 0; round < rounds; ++round)
    {
      std::shuffle(order.begin(), order.end(), random);
      for (const size_t index : order)
      {
        const Product multiply = methods[index].multiply;
        // Establish steady-state bursts for each algorithm before timing.
        for (unsigned j = 0; j < 4; ++j)
          checksum += multiply(a, b)[j];
        const auto start = std::chrono::steady_clock::now();
        for (unsigned j = 0; j < reps; ++j)
          checksum += multiply(a, b)[j % expected.size()];
        const auto elapsed = std::chrono::steady_clock::now() - start;
        samples[index][round] = std::chrono::duration<double, std::nano>(elapsed).count() / reps;
      }
    }

    for (size_t i = 0; i < methods.size(); ++i)
    {
      std::array<double, rounds> speedup{}, recursiveSpeedup{};
      unsigned wins = 0;
      for (unsigned j = 0; j < rounds; ++j)
      {
        speedup[j] = samples[0][j] / samples[i][j];
        recursiveSpeedup[j] = samples[1][j] / samples[i][j];
        wins += speedup[j] > 1.0;
      }
      std::cout << shape.name << ',' << longer << ',' << shorter << ',' << methods[i].name << ','
                << rounds << ',' << reps << ',' << Median(samples[i]) << ',' << Median(speedup) << ','
                << Median(recursiveSpeedup) << ',' << wins << '\n';
    }
    return true;
  }
}

int main(int argc, char **argv)
{
  bool heldOut = false, narrow = false;
  for (int i = 1; i < argc; ++i)
  {
    const std::string_view arg = argv[i];
    if (arg == "--holdout")
      heldOut = true;
    else if (arg == "--balanced-band")
      narrow = true;
    else
    {
      std::cerr << "usage: " << argv[0] << " [--holdout] [--balanced-band]\n";
      return 2;
    }
  }
  std::span<const SizeT> lengths = heldOut ? holdout : sweep;
  auto selectedShapes = std::span<const Shape>(shapes);
  if (narrow)
  {
    lengths = heldOut ? std::span<const SizeT>(balancedHoldout) : std::span<const SizeT>(balancedBand);
    selectedShapes = selectedShapes.first(1);
  }
  std::mt19937_64 random((heldOut ? 0x91b85 : 0x1a5685) + unsigned(narrow));
  DataT checksum = 0;
  std::cout << std::fixed << std::setprecision(4)
            << "shape,longer_limbs,shorter_limbs,algorithm,rounds,reps,median_ns,"
               "speedup_vs_dispatch,speedup_vs_recursive,wins_vs_dispatch\n";
  for (Shape shape : selectedShapes)
    for (SizeT length : lengths)
      if (!RunCase(shape, length, random, checksum))
        return 1;
  std::cerr << "checksum=" << checksum << '\n';
}
