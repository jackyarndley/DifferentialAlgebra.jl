# DifferentialAlgebra.jl

[DifferentialAlgebra.jl](https://github.com/jackyarndley/DifferentialAlgebra.jl) computes multivariate Taylor
polynomials using a Julia engine with `Float32`, `Float64` and `BigFloat`
coefficients. It requires Julia 1.10 or later.

## Getting started

Install from GitHub until the package is registered:

```julia
using Pkg
Pkg.add(url="https://github.com/jackyarndley/DifferentialAlgebra.jl")
```

For a local checkout, follow [Development](tutorials/setting-up-your-development-environment.md).

```julia
using DifferentialAlgebra
DifferentialAlgebra.init(6, 2)
x, y = DifferentialAlgebra.identity()
p = sin(x) * exp(y)
DifferentialAlgebra.evaluate(p, [0.1, 0.2])
DifferentialAlgebra.getCoefficient(p, [1, 1])
```

`DA(c)` creates a constant, including when `c` is an integer. `DA(i, c)` creates
`c` times independent variable `i`; index zero creates a constant. Use
`DifferentialAlgebra.identity(Float32)` or `DifferentialAlgebra.identity(BigFloat)` to choose coefficient types.

Polynomial evaluation uses `DifferentialAlgebra.evaluate`; `evalScalar` remains available.
See [Julia engine and precision](tutorials/native-julia.md) for details.
