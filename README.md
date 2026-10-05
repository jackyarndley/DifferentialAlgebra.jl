# DifferentialAlgebra.jl

[![Tests](https://github.com/jackyarndley/DifferentialAlgebra.jl/actions/workflows/test.yml/badge.svg)](https://github.com/jackyarndley/DifferentialAlgebra.jl/actions/workflows/test.yml)
[![Documentation](https://img.shields.io/badge/docs-online-blue.svg)](https://jackyarndley.github.io/DifferentialAlgebra.jl/)

Multivariate Taylor polynomials in Julia for derivatives, polynomial maps and
numerical integration. Coefficients are parameterized by a real scalar type;
operations use that type's arithmetic and mathematical functions.

## Installation

Requires Julia 1.10 or later. Install directly from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/jackyarndley/DifferentialAlgebra.jl")
```

## Quick start

```julia
using DifferentialAlgebra

x, y = variables((:x, :y); order = 6)
p = sin(x) * exp(y)

p([0.1, 0.2])                  # Evaluate the polynomial
coefficient(p, [1, 1])          # Coefficient of x*y
differentiate(p, 1)             # Partial derivative with respect to x
constant_term(p)                # Value at the expansion point

map = CompiledMap([p, x + y])   # Reuse the evaluation tree
map([0.1, 0.2])
```

`variables(T, n; order)` selects a coefficient type. Creating a new algebra
invalidates existing polynomials; see the [user guide](https://jackyarndley.github.io/DifferentialAlgebra.jl/manual/getting-started/).

Loading optional IntervalArithmetic enables `enclose(p, box)` for certified
stored-polynomial ranges and native `taylor_models(box; order)` with absolute
function remainders. Interval coefficients alone do not certify truncation error.
See [intervals and Taylor models](docs/src/manual/interval-models.md).
`adaptive_map(f, box; estimator=IntervalBound(), order, atol)` reevaluates original
static maps and uses uniform interval error bounds for splitting.
`validated_adaptive_map` remains a convenience entry point. The
[interval fitting example](examples/interval_models.jl) plots enclosures,
partitions, and accepted error bounds. `splitter=:oriented` adds 2D convex polygon
ADS with automatic or supplied projection directions to all four error methods;
see the [polygon fitting plots](examples/polygon_ads.jl). GuardedTail,
ExtrapolatedTail, LastTerms and time solvers keep their heuristic contracts.

For smooth optimization objectives, `continuous_map(f, fit; continuity=:c2)`
optionally blends fresh overlapping fits with C0/C1/C2 continuity and supports
ForwardDiff gradients/Hessians. Interval sources retain function certificates
and report the new overlap error bounds. See the [continuity plots](examples/ads_continuity.jl),
[oriented optimization example](examples/ads_optimization.jl) and
[accuracy contract](docs/src/manual/continuous-ads.md).

- [Documentation](https://jackyarndley.github.io/DifferentialAlgebra.jl/)
- [Runnable examples](examples)
- [Contributing](docs/src/contributing.md)
