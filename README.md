# DifferentialAlgebra.jl

[![Tests](https://github.com/jackyarndley/DifferentialAlgebra.jl/actions/workflows/test.yml/badge.svg)](https://github.com/jackyarndley/DifferentialAlgebra.jl/actions/workflows/test.yml)
[![Documentation](https://github.com/jackyarndley/DifferentialAlgebra.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/jackyarndley/DifferentialAlgebra.jl/actions/workflows/documentation.yml)

Multivariate Taylor polynomials implemented in Julia, with `Float32`, `Float64`
and `BigFloat` coefficients. Includes arithmetic, elementary and special functions,
derivatives, substitutions, map inversion, symmetric eigenpairs and moments.
Requires Julia 1.10 or later.

## Install

Until the package is registered, install directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/jackyarndley/DifferentialAlgebra.jl")
```

## Quick start

```julia
using DifferentialAlgebra

DifferentialAlgebra.init(6, 2)             # Order 6, two independent variables
x, y = DifferentialAlgebra.identity()
p = sin(x) * exp(y)
DifferentialAlgebra.evaluate(p, [0.1, 0.2])
DifferentialAlgebra.getCoefficient(p, [1, 1])  # 1.0

x32 = DifferentialAlgebra.variable(1, Float32)
sin(x32)                                # DA{Float32}

setprecision(256) do
    xbig = DifferentialAlgebra.variable(1, BigFloat)
    DifferentialAlgebra.evaluate(exp(xbig), BigFloat[big"0.1", 0])
end
```

`DA(c)` creates a constant. `DA(i, c)` creates `c` times independent variable `i`.
Calling `init` starts a new global algebra and invalidates existing polynomials.

- [Engine and precision guide](docs/src/tutorials/native-julia.md)
- [Examples](examples): elementary functions, gradients, map inversion and orbit integration
- [API coverage and compatibility](docs/src/api-coverage.md)
- [Reproducible TaylorSeries.jl integration comparison](benchmark)
- [Development and documentation](docs/src/tutorials/setting-up-your-development-environment.md)

## Development

From this checkout, run `julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'`.
CI runs on Linux with Julia 1.10 and the latest stable Julia, and executes the
documented examples. Benchmarks and examples have separate dependency environments.

## Acknowledgments

This package grew out of the native Julia implementation developed for
[DACE.jl](https://github.com/UoA-AstroGroup/DACE.jl), retaining its mathematical
API, tests and MIT attribution. The polynomial engine uses no DACE binary,
CxxWrap or Eigen. DACE text and binary formats remain available for interoperability;
special functions use SpecialFunctions.jl and QuadGK.jl.
