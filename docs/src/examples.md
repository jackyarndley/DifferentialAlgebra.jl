# Examples

Each tutorial is a standalone Julia script. The documentation executes the
scripts with Literate.jl, showing their numerical results and Makie figures.

## Taylor arithmetic and calculus

- [Taylor expansions and truncation](generated/taylor_expansions.md)
- [Rational functions and convergence](generated/rational_functions.md)
- [Differentiation and integration](generated/calculus.md)
- [Multivariate derivatives and removable singularities](generated/multivariate_calculus.md)
- [Inverse maps](generated/inverse_maps.md)
- [Polynomial linear algebra](generated/linear_algebra.md)
- [Orthogonal polynomials](generated/orthogonal_polynomials.md)

## Implicit equations

- [Implicit Kepler maps](generated/implicit_maps.md)
- [Newton and fixed-point iteration](generated/implicit_solvers.md)
- [Pseudo-arclength continuation](generated/continuation.md)

## Numerical integration

- [Adaptive Runge–Kutta](generated/adaptive_integration.md)
- [Orbit integration and state transition matrices](generated/ode_integration.md)
- [Taylor integration in time](generated/taylor_integration.md)
- [Monte Carlo with a Taylor map](generated/damc_kepler.md)

## Automatic domain splitting

- [Comparing estimators and split directions](generated/ads_methods.md)
- [An analytic Kepler map](generated/ads_kepler.md)
- [Splitting during orbit propagation](generated/ads_flow.md)

## Storage

- [Serializing compiled maps](generated/serialization.md)

## Run locally

From the repository root, install the example environment:

```sh
julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=examples examples/ads_methods.jl
```

For interactive plotting, start Julia with `julia --project=examples` and run:

```julia
include("examples/ads_methods.jl")
```

Each plot ends with `fig`, so the interactive frontend or documentation
renderer displays it. A terminal script runs the calculations and prints
results without opening a viewer or writing image files.
