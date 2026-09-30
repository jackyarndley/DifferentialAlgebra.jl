# Getting started

```@meta
CurrentModule = DifferentialAlgebra
```

## Create variables

```@example basics
using DifferentialAlgebra
x, y = variables(2; order=4)
p = 2 + 3x + x*y + y^2
```

[`variables`](@ref) initializes an algebra and returns its independent variables
in an ordinary Julia vector. Every monomial has total degree at most `order`;
multiplication and analytic functions discard higher-degree terms.

To expand about a nonzero point, add that point to the independent variables:
`sin(2 + x)` expands around 2. To retrieve a variable without reinitializing
the algebra, use [`variable`](@ref).

`DA(c)` creates a constant polynomial. In particular, `DA(1)` is constant one.

## Display polynomials

Polynomials print as expressions in `x1`, `x2`, and so on. These labels follow
the independent-variable indices, regardless of the Julia names assigned to them.

```@example basics
2 + 3x - y^2
```

Array displays use the same expressions. Julia's limited displays abbreviate
long polynomials with `…`; `print(p)` or `string(p)` includes every term.
Display is for reading, rather than a serialization format. Use
`DifferentialAlgebra.getMonomials(p)` to inspect coefficient and exponent data.

## Inspect and evaluate

```@example basics
constant_term(p), coefficient(p, [1, 1]), p([0.2, -0.1])
```

[`coefficient`](@ref) takes one nonnegative exponent per variable. It returns
the coefficient of the ordinary monomial, without factorial scaling.
`p(point)` is equivalent to [`evaluate(p, point)`](@ref evaluate).

## Differentiate and integrate

```@example basics
dx = differentiate(p, 1)
mixed = differentiate(p, [1, 1])
primitive = integrate(p, 2)
(constant_term(dx), constant_term(mixed), primitive)
```

Variable indices start at one. A vector of counts requests repeated partial
derivatives. Integration chooses a zero integration constant and truncates to
the current order.

[`gradient`](@ref), [`jacobian`](@ref) and [`hessian`](@ref) return arrays of
polynomials. Apply `constant_term` to obtain their values at the expansion point.

```@example basics
constant_term(gradient(p))
```

## Compare polynomials

Scalar comparisons use the constant part. To compare every coefficient, use
`DifferentialAlgebra.norm(p - q)`; `iszero(p)` also checks the full polynomial.

Calling `variables` again starts a new algebra and invalidates existing
polynomials. See [Configuration and storage](configuration.md) when managing
several calculations.
