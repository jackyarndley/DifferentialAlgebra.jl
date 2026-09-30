# DifferentialAlgebra.jl

DifferentialAlgebra.jl computes multivariate Taylor expansions using ordinary
Julia arithmetic. A polynomial carries its coefficients through a calculation,
so its output describes both a value and its local dependence on the inputs.

Use it to differentiate expressions, compose and invert polynomial maps, or
propagate an expansion through a numerical integrator. The coefficient type is
a parameter of `TaylorPolynomial{T}`; see [Coefficient types](manual/coefficient-types.md).

## Installation

Julia 1.10 or later is required.

```julia
using Pkg
Pkg.add(url = "https://github.com/jackyarndley/DifferentialAlgebra.jl")
```

## A first calculation

```@example introduction
using DifferentialAlgebra
x, y = variables(2; order=6)
p = sin(x) * exp(y)
(coefficient(p, [1, 1]), p([0.1, 0.2]))
```

Here `x` and `y` are independent perturbations about zero. The algebra retains
monomials through total degree six.

Start with [Getting started](manual/getting-started.md), then explore the
[examples](generated/sine.md) or the [API reference](api.md).
