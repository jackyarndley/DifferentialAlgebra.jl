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

- [Documentation](https://jackyarndley.github.io/DifferentialAlgebra.jl/)
- [Runnable examples](examples)
- [Integration benchmarks](benchmark)
- [Contributing](docs/src/contributing.md)

## Acknowledgments

DifferentialAlgebra.jl grew out of the native Julia implementation developed for
[DACE.jl](https://github.com/UoA-AstroGroup/DACE.jl). It retains its MIT attribution.
