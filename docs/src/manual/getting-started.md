# Getting started

```@meta
CurrentModule = DifferentialAlgebra
```

## Create variables

```@example basics
using DifferentialAlgebra
x, y = variables((:x, :y); order=4)
p = 2 + 3x + x*y + y^2
```

[`variables`](@ref) initializes an algebra and returns its independent variables
in an ordinary Julia vector. Every monomial has total degree at most `order`;
multiplication and analytic functions discard higher-degree terms.

Names may be symbols or strings in a tuple or vector. They must be distinct
identifiers and are displayed exactly as supplied, including Unicode names such
as `:δx` or `:q₁`. The number of variables is inferred from the names. You can
also write `variables(2; order=4, names=("x", "y"))` or choose a coefficient type
with `variables(T, (:x, :y); order=4)`.

For automatic names, use `variables(2; order=4)`. Assigning the returned values
to Julia bindings does not rename the polynomial coordinates. Names are copied
at initialization and affect display only; coefficients and derivatives still
use variable indices.

To expand about a nonzero point, add that point to the independent variables:
`sin(2 + x)` expands around 2. To retrieve a variable without reinitializing
the algebra, use [`variable`](@ref).

`TaylorPolynomial(c)` creates a constant polynomial. In particular, `TaylorPolynomial(1)` is constant one.

## Display polynomials

Polynomials print with Unicode superscripts for powers, such as `x²` and `y¹²`.
Automatic variable names use subscripts: `x₁`, `x₂`, …, `x₁₀`. With the names
chosen above, the following expression displays as `2.0 + 3.0 x - 1.0 y²`:

```@example basics
2 + 3x - y^2
```

Array displays use the same expressions. Julia's limited displays abbreviate
long polynomials with `…`; `print(p)` or `string(p)` includes every term.
Display is for reading, rather than a serialization format. Use
`DifferentialAlgebra.monomials(p)` to inspect coefficient and exponent data.

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
`DifferentialAlgebra.coefficient_norm(p - q)`; `iszero(p)` also checks the full polynomial.

Calling `variables` again starts a new algebra and invalidates existing
polynomials. See [Configuration and storage](configuration.md) when managing
several calculations.
