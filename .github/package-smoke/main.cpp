#include <biginteger/BigInteger.h>
#include <biginteger/algorithms/Division.h>
#include <biginteger/algorithms/Squaring.h>
#include <biginteger/common/Parser.h>
#include <biginteger/ops/Operations.h>

using BigMath::Parse;

int main()
{
    if (BigMath::DIV_BZ_MIN_B != BIGMATH_DIV_BZ_MIN_B ||
        BigMath::NTT_SQUARE_THRESHOLD != BIGMATH_NTT_SQUARE_THRESHOLD)
        return 1;

    const auto a = Parse("18446744073709551617");
    const auto b = Parse("3");
    const auto product = a * b;
    if (product != Parse("55340232221128654851") || product / b != a)
        return 1;

    if (BigMath::PowMod(Parse("7"), Parse("65537"), a) != Parse("8961139343211450995") ||
        BigMath::Gcd(product, a) != a)
        return 1;

    const auto [divisor, x, y] = BigMath::ExtendedGcd(a, b);
    const auto [root, remainder] = BigMath::SqrtRem(a);
    return divisor == Parse("1") && a * x + b * y == divisor &&
                   root == Parse("4294967296") && remainder == Parse("1")
               ? 0
               : 1;
}
