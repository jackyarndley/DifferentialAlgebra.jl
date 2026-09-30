# Examples

From the repository root, install the example dependencies and develop this checkout:

```sh
julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=examples examples/sine.jl
julia --project=examples examples/gradient.jl
julia --project=examples examples/polynomial_inversion.jl
julia --project=examples examples/ode_integration.jl
```

Each script initializes its own algebra and checks its result. They cover
elementary functions, derivatives, inverse maps and orbit integration with
OrdinaryDiffEq's Vern9 solver. The documentation build also executes all four.

The [benchmark](../benchmark) compares the same orbit integration using
DifferentialAlgebra.jl and TaylorSeries.jl, including Float32 and BigFloat.
