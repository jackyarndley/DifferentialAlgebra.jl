# Coefficient types

`TaylorPolynomial{T}` stores coefficients of a concrete real scalar type.
Choose `T` to match the precision and arithmetic needed by the calculation;
`variables(n; order)` defaults to `Float64`.

## Choose the scalar arithmetic

```@example exact
using DifferentialAlgebra
x, y = variables(Rational{BigInt}, 2; order=4)
p = (1 + x/3 + y/2)^2
coefficient(p, [1, 1])
```

Use `variables(BigFloat, n; order)` for higher floating precision and construct
all inputs inside the intended `setprecision` context. Fixed-width integer types
retain Julia's overflow behavior, including powers in affine substitutions;
`BigInt` or `Rational{BigInt}` avoids that overflow in exact arithmetic.

Arithmetic follows scalar promotion rules. Division and analytic functions may
change the coefficient type. `TaylorPolynomial(c)` uses `typeof(float(c))`;
`TaylorPolynomial{T}(c)` preserves an explicitly chosen type. Explicit polynomial
constants can also disambiguate custom number types with broad `Real` methods.

## Conversion and scalar queries

Converting to `Float32`, `Float64` or `BigFloat` requires a constant polynomial;
a nonconstant input raises `InexactError`. Use `constant_term(p)` to discard
higher terms explicitly. Conversion between polynomial coefficient types retains
stored terms, trimming only exact zeros introduced by conversion.

| Intent | Use |
|:--|:--|
| Inspect the expansion-center value | `constant_term(p)` |
| Select a larger coefficient type | `TaylorPolynomial{BigFloat}(p)` or `TaylorPolynomial{BigInt}(p)` |
| Round a scalar center | `round(constant_term(p))` |
| Construct rational coefficients | `variables(Rational{BigInt}, n; order=...)` |
| Divide Taylor expansions | `p/q` |

Polynomial `Int`, integer `floor`, `round`, `trunc`, `big`, `ldexp`, `//` and
`fma` are unsupported. Use the scalar operation explicitly or ordinary polynomial
scaling/division. `muladd` is supported: it filters the completed product and sum,
but promises no single fused rounding.

Comparisons and `sign` inspect constant parts; structural equality and full-term
predicates are described in [Getting started](getting-started.md#Compare-polynomials).

## Operation requirements

| Operation | Required scalar behavior |
|:--|:--|
| Construction | Concrete storage, `zero`, `one`, conversion |
| Polynomial arithmetic | `+`, `-`, `*`, `muladd`, `iszero`, promotion |
| Division/integration | Division and conversion of integer factors |
| Elementary/special functions | Corresponding scalar methods at the center |
| Map inversion/eigenpairs | Required linear algebra/factorization methods |

`T <: Real` alone does not supply these methods. Expansions also require regular
real centers: singular reciprocals, gamma poles and invalid real logarithms raise
`TaylorError`; invalid API parameters raise `ArgumentError`. NaN centers propagate.
Zero to a varying power is supported only around a positive exponent; negative
scalar bases require a constant exponent with a defined real scalar power.

Loading ForwardDiff 1 enables mixed dual/polynomial arithmetic, preserving partials
in coefficients. Loading IntervalArithmetic enables outward coefficient arithmetic
and Taylor models. Their [guarantee rules](interval-models.md) remain distinct
from the ordinary truncation contract.
