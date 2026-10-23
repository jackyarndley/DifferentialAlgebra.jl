# DifferentialAlgebra.jl

[![Tests](https://github.com/jackyarndley/DifferentialAlgebra.jl/actions/workflows/test.yml/badge.svg)](https://github.com/jackyarndley/DifferentialAlgebra.jl/actions/workflows/test.yml)
[![Documentation](https://img.shields.io/badge/docs-online-blue.svg)](https://jackyarndley.github.io/DifferentialAlgebra.jl/)

Multivariate Taylor arithmetic for derivatives, local inverse maps, uncertainty
propagation and numerical integration. Write a calculation with polynomial inputs
to obtain its dependence on those inputs through a chosen total degree.

## Install

Requires Julia 1.12 or later:

```julia
using Pkg
Pkg.add(url = "https://github.com/jackyarndley/DifferentialAlgebra.jl")
```

## One calculation, many nearby inputs

```julia
using DifferentialAlgebra

δx, δy = variables((:δx, :δy); order = 6)
p = sin(1 + δx) * exp(δy)      # Expand around (x, y) = (1, 0)

p([0.1, 0.2])                 # Predict the value at (1.1, 0.2)
constant_term(gradient(p))      # Gradient at the expansion point
constant_term(hessian(p))       # Hessian at the expansion point
coefficient(p, [1, 1])         # Coefficient of δx*δy

map = CompiledMap([p, δx + δy])
map([0.1, 0.2])                # Reuse a shared evaluation tree
```

Coordinates are perturbations about the expansion point. Calling `variables`
again creates a new algebra and invalidates its existing polynomials; compiled
maps remain usable for numeric evaluation.

## Choose the calculation

| Need | Interface | Error contract |
|:--|:--|:--|
| Derivatives or a local map | `TaylorPolynomial`, `CompiledMap` | Truncated polynomial |
| A larger input region | `adaptive_map` | Heuristic with `GuardedTail`, `ExtrapolatedTail` or `LastTerms` |
| Original-function enclosures | `TaylorModel`, `IntervalBound` | Outward interval arithmetic plus an absolute remainder |
| Smooth ADS for optimization | `continuous_map(...; continuity=:c2)` | Smooth surrogate; derivative accuracy requires separate checks |
| ODE propagation | `taylor_expand`, SciML `TaylorMethod` | Independent time/uncertainty orders; numerical error estimate |

IntervalArithmetic, ForwardDiff and SciML integrations load only when requested.
Interval coefficients alone do not bound discarded polynomial terms; use Taylor
models when a function enclosure is required.

## Learn by solving a problem

- [Getting started](https://jackyarndley.github.io/DifferentialAlgebra.jl/manual/getting-started/)
  and [API reference](https://jackyarndley.github.io/DifferentialAlgebra.jl/api/).
- [Orbit maps and state transition matrices](examples/ode_integration.jl),
  [six-variable uncertainty](examples/six_variable_ephemeris.jl), and
  [angles-only orbit determination](examples/daiod.jl).
- [Certified interval fitting](examples/interval_models.jl) and
  [C2 surrogate optimization](examples/ads_optimization.jl).
- [All executable tutorials](docs/src/examples.md),
  [migration guide](docs/src/migration.md), and [contributing](docs/src/contributing.md).
