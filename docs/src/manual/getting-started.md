# Getting started

```@meta
CurrentModule = DifferentialAlgebra
```

A polynomial carries a value and its dependence on small input perturbations.
First choose the perturbations and maximum total degree; then evaluate the same
expression you would use with numbers.

## Create variables

```@example basics
using DifferentialAlgebra
δx, δy = variables((:δx, :δy); order=4)
p = sin(1 + δx) * exp(δy)
```

This expands around `(x, y) = (1, 0)`. [`variables`](@ref) returns a normal
Julia vector. Multiplication and analytic functions discard terms above total
degree four. To choose the coefficient type, use
`variables(BigFloat, (:δx, :δy); order=4)`; the default is `Float64`.

Use `variables(2; order=4)` for automatic names. Names are distinct identifiers,
copied for display only; coefficient and derivative queries use indices.
[`variable(1)`](@ref variable) retrieves a coordinate without reinitializing.
`TaylorPolynomial(1)` constructs the constant one, rather than a variable.

**Calling `variables` again invalidates existing polynomials**, even with the
same dimension and order. Finish a calculation or compile its numeric map before
starting another algebra. [Configuration](configuration.md) covers ownership,
working orders and concurrent calculations.

## Inspect and evaluate

```@example basics
(value_at_center = constant_term(p),
 mixed_coefficient = coefficient(p, [1, 1]),
 nearby_value = p([0.1, 0.2]))
```

[`coefficient`](@ref) takes one nonnegative exponent per variable and returns the
ordinary monomial coefficient: `c_α = ∂^α f / α!`. `p(point)` is equivalent to
[`evaluate(p, point)`](@ref evaluate); the input is a perturbation, so `[0.1, 0.2]`
corresponds to the physical point `(1.1, 0.2)` here.

Printed expressions use the chosen names and Unicode powers. `print(p)` shows
all terms; limited array displays may abbreviate them. Use [`monomials`](@ref)
for coefficient/exponent data rather than parsing display text.

## Differentiate and integrate

```@example basics
dx = differentiate(p, 1)
mixed = differentiate(p, [1, 1])
primitive = integrate(p, 2)
(gradient = constant_term(gradient(p)), hessian = constant_term(hessian(p)),
 mixed_derivative = constant_term(mixed))
```

Indices start at one; a vector of counts requests repeated partial derivatives.
`gradient`, `jacobian` and `hessian` return arrays of polynomials, which can be
evaluated away from the center. Apply `constant_term` for their center values.
Integration chooses a zero integration constant and obeys the working order.

## Reuse the calculation

```@example basics
map = CompiledMap([p, δx + δy])
map([0.1, 0.2])
```

A compiled map owns a shared evaluation tree and remains numerically evaluable
after algebra reinitialization. See [Polynomial maps](maps.md) for composition,
inversion and buffered evaluation.

## Compare polynomials

`==`, `<`, `<=`, `>`, `>=`, `isless`, `sign` and `signbit` inspect the constant
part so existing scalar branches can run on polynomials. They do not establish
equality of expansions. `isequal` and `hash` inspect the complete stored
coefficient sequence and dimension; avoid mutating dictionary keys.

`iszero`, `isone`, `isinteger` and `coefficient_norm` inspect all terms.
However, `coefficient_norm(p-q)` first applies the current arithmetic order and
coefficient tolerance. For an unfiltered comparison, extract matching scalar
coefficients with `coefficient` and compare them directly.
