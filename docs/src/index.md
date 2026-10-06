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
[tutorials grouped by topic](examples.md) or the [API reference](api.md).

## Astrodynamics examples

Explore [six-variable orbital-element uncertainty](generated/six_variable_ephemeris.md)
with linear and cubic maps, analytic polynomial moments and interval enclosures.
Follow [injection uncertainty under J2 gravity](generated/j2_uncertainty.md) through
inclined orbits, ground tracks and local-frame dispersion.

For ADS, compare the four error methods on an
[inclined orbital ribbon with polygon splitting](generated/orbit_polygon_ads.md),
or classify [B-plane flyby uncertainty](generated/flyby_uncertainty.md) using
guaranteed altitude bounds. Static-function certificates and numerical orbit
propagation have separate error contracts throughout these tutorials.

![Polygon uncertainty mapped into an inclined orbital ribbon](generated/orbit_polygon_ads-figure-2.png)
