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
    return product == Parse("55340232221128654851") &&
                   product / b == a
               ? 0
               : 1;
}
