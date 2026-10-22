# Coefficient types

Converting a polynomial to `Float32`, `Float64`, or `BigFloat` requires a constant
polynomial. Nonconstant inputs throw `InexactError` rather than discard terms.
Use `constant_term(p)` explicitly when only the value at the expansion origin
is needed. Converting between polynomial coefficient types preserves the stored
terms, trimming only exact zeros introduced by conversion.

`TaylorPolynomial{T}` is parameterized by a real scalar type `T <: Real`. This includes the
`AbstractFloat` hierarchy and allows other real number types. A subtype bound
alone does not supply arithmetic: each operation requires the corresponding
scalar methods and conversions.

## Choose a type

Use `variables(T, n; order)` to choose the coefficient type explicitly.
`variables(n; order)` defaults to `Float64`. Construct constants with
`TaylorPolynomial{T}(c)` when their storage type should be explicit.

For example, rational coefficients can preserve exact polynomial arithmetic:

```@example exact
using DifferentialAlgebra
x, y = variables(Rational{Int}, 2; order=4)
p = (1 + x/3 + y/2)^2
coefficient(p, [1, 1])
```

Arithmetic follows Julia's promotion rules. Division or an elementary function
may produce a different coefficient type when its scalar result requires one.
Fixed-width integer coefficients have Julia's usual overflow behavior, including
the scalar powers in affine substitutions. Exact binomial factors do not prevent
those powers from overflowing. Use `BigInt` or `Rational{BigInt}` for unbounded
exact arithmetic; for example, scaling `x^19` by 10 requires a coefficient of
`big(10)^19`, which exceeds `typemax(Int64)`.
The convenience constructor `TaylorPolynomial(c)` uses `typeof(float(c))`; use `TaylorPolynomial{T}(c)`
to retain an exact scalar type.
Explicit polynomial constants also avoid competing mixed-scalar methods when
a custom number type defines broad operations on `Real` values.

## Scalar queries and unsupported operations

`sign(p)` and `signbit(p)` inspect the constant term, as branch comparisons do.
`iszero`, `isone` and `isinteger` inspect all stored coefficients. Operations that
need a scalar integer, including `Int(p)` and `floor(Int, p)`, are unsupported;
apply them to `constant_term(p)` explicitly. Polynomial `round` and `trunc` are
also unsupported.

Use `TaylorPolynomial{BigFloat}(p)` or `TaylorPolynomial{BigInt}(p)` to choose a
larger coefficient type; `big(p)` is unsupported. Use scalar multiplication in
place of `ldexp(p, n)`. Construct rational coefficients explicitly and use `/`
for Taylor-series division; `//` does not construct rational functions of
polynomials. `muladd` evaluates the polynomial product and sum before final
coefficient filtering. `fma` is unsupported: polynomial multiplication and
coefficient accumulation do not promise a single fused rounding.

Loading ForwardDiff enables mixed `Dual` and polynomial arithmetic. Dual
partials remain in each coefficient; only domain and branch checks inspect the
primal value.

## Scalar requirements

| Operation | Requirements on the coefficient type |
|:--|:--|
| Construction | `zero`, `one`, conversion and a concrete storage type |
| Polynomial arithmetic | Addition, subtraction, multiplication, `muladd`, `iszero` and promotion |
| Division and integration | Scalar division and conversion of integer factors |
| Elementary functions | The corresponding scalar function at the expansion point |
| Special functions | Suitable scalar methods, generally provided by SpecialFunctions.jl |
| Eigenpairs and map inversion | The required scalar linear algebra and factorization methods |

An analytic expansion also requires a regular center, such as a nonzero
constant for a reciprocal. Numerical precision and range are properties of the
chosen scalar type. An undefined real analytic expansion raises `TaylorError`,
including gamma poles, negative real logarithms, and singular time-series centers.
NaN constants propagate through polynomial elementary and special functions.
Invalid orders and other API parameters raise `ArgumentError`. A zero scalar base
raised to a polynomial or time series is identically zero near a positive exponent;
a varying exponent centered at zero is rejected. Negative scalar bases require a
constant exponent with a defined real scalar power.

For configurable-precision types, construct coefficients
and perform calculations within that type's precision context.

Loading IntervalArithmetic adds narrow coefficient hooks for exact-zero tests,
internal integer factors and interval accumulation. See
[intervals and Taylor models](interval-models.md) for guarantee flags, supported
operations and the distinction between interval coefficients and function remainders.
