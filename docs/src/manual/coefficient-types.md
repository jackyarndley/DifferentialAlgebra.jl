# Coefficient types

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
The convenience constructor `TaylorPolynomial(c)` uses `typeof(float(c))`; use `TaylorPolynomial{T}(c)`
to retain an exact scalar type.
Explicit polynomial constants also avoid competing mixed-scalar methods when
a custom number type defines broad operations on `Real` values.

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
chosen scalar type. For configurable-precision types, construct coefficients
and perform calculations within that type's precision context.

Loading IntervalArithmetic adds narrow coefficient hooks for exact-zero tests,
internal integer factors and interval accumulation. See
[intervals and Taylor models](interval-models.md) for guarantee flags, supported
operations and the distinction between interval coefficients and function remainders.
