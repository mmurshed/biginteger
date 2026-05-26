# BigDecimal

Arbitrary-precision signed decimal numbers, built on top of `BigInteger`. Fixed-point Java-style semantics: every value is `unscaled · 10^(-scale)`.

---

## Number model

```
value  =  unscaled · 10^(-scale)

  unscaled : BigInteger     (signed)
  scale    : int            (positive → fractional digits; negative → trailing zeros)
```

Examples:

| value          | unscaled | scale |
|----------------|---------:|------:|
| `0`            |        0 |     0 |
| `123`          |      123 |     0 |
| `123.456`      |   123456 |     3 |
| `0.001`        |        1 |     3 |
| `1500`         |       15 |    −2 |
| `-3.14`        |     −314 |     2 |

The scale is preserved by every operation. `0 * 1.23` is `0.00` (scale 2), not `0`. `1.50 + 0.50` is `2.00` (scale 2). This matches `java.math.BigDecimal`.

## API

```cpp
#include "bigdecimal/BigDecimal.h"

using namespace BigMath;

BigDecimal a = BigDecimal::FromString("3.14159265358979323846");
BigDecimal b = BigDecimal::FromString("2.71828182845904523536");

BigDecimal sum  = a + b;                 // exact
BigDecimal diff = a - b;                 // exact
BigDecimal prod = a * b;                 // exact
BigDecimal quot = a.Divide(b, 30, RoundingMode::HALF_EVEN);  // explicit scale + mode

std::cout << prod.ToPlainString() << "\n";
```

| operation                                | exact | notes |
|------------------------------------------|-------|-------|
| `Add`, `Subtract`                        | yes   | result scale = max of input scales |
| `Multiply`                               | yes   | result scale = sum of input scales |
| `Divide(o, scale, RoundingMode)`         | no    | result has exactly `scale` fractional digits |
| `SetScale(newScale, RoundingMode)`       | up: yes; down: rounded | |
| `StripTrailingZeros()`                   | yes (value) | reduces scale while value unchanged |
| `Negate`, `Abs`                          | yes   | |
| `CompareTo`                              | —     | value-equal across scales |
| `Equals`                                 | —     | strict: matches both unscaled and scale |
| `ToPlainString`                          | —     | no scientific notation in output |

### Parse grammar

```
  [+-]? ( digits ('.' digits?)? | '.' digits ) ( [eE] [+-]? digits )?
```

Examples: `123`, `-123.456`, `1.5e2`, `2e-3`, `.5`, `+1.25`.

### Rounding modes

| mode | behavior |
|---|---|
| `UP`           | away from zero |
| `DOWN`         | toward zero (truncate) |
| `CEILING`      | toward +∞ |
| `FLOOR`        | toward −∞ |
| `HALF_UP`      | nearest; ties away from zero |
| `HALF_DOWN`    | nearest; ties toward zero |
| `HALF_EVEN`    | nearest; ties to even (banker's) |
| `UNNECESSARY`  | throws `std::runtime_error` if rounding would discard a non-zero remainder |

`Divide` and `SetScale` are the only operations that take a rounding mode.

## Performance

Every BigDecimal operation reduces to a small constant number of BigInteger operations plus rescaling by `10^k`. Rescaling uses `BigInteger::operator*` against a cached `10^k`, which dispatches through the same Classic / Karatsuba / NTT path as any other multiplication.

Measured on Apple M1 Max (`-O3 -march=native`):

| op (operand `int.frac`)  | time |
|---|---:|
| Add `20000.2000`        | 0.008 ms |
| Multiply `20000.2000`    | 1.6 ms |
| Divide `20000.2000`, scale 10 | 0.79 ms |
| Parse 10 000 digits      | 0.33 ms |
| ToString 10 000 digits   | 1.66 ms |

Add cost ≈ BigInteger add of rescaled unscaled values. Multiply cost ≈ BigInteger multiply at the combined unscaled size. Divide cost ≈ one BigInteger multiplication (to apply the rescale) plus one BigInteger divmod plus a small rounding adjustment.

The BigDecimal layer adds no quadratic overhead — every cost line above is dominated by the underlying BigInteger op at the same scale.

### Caching

`Pow10` (limb-form) is cached thread-locally inside `Parser.h`. BigDecimal additionally caches the `BigInteger` wrapper for `10^n` keyed by `n`, so repeated rescaling at the same scale skips both the vector build and the wrapper copy.

## See also

- [BASE.md](BASE.md) — number representation (BigDecimal builds on this)
- [MULTIPLICATION.md](MULTIPLICATION.md) — drives `Multiply` and `Pow10` cost
- [DIVISION.md](DIVISION.md) — drives `Divide`
- [STRING_CONVERSION.md](STRING_CONVERSION.md) — drives `FromString` / `ToPlainString`
