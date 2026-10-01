# Examples

The first five scripts reproduce the examples in
[DACE.jl](https://github.com/UoA-AstroGroup/DACE.jl/tree/7271c372ce86350e80a95a9c2c0d7e8739d09e83/examples)
using DifferentialAlgebra.jl's native Julia API. Each script is self-contained
and includes numerical checks. An additional Kepler example demonstrates
automatic domain splitting.

From the repository root:

```sh
julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=examples examples/sine.jl
julia --project=examples examples/polynomial_inversion.jl
julia --project=examples examples/tutorial1ex8.jl
julia --project=examples examples/ode_integration.jl
julia --project=examples examples/damc_kepler.jl
julia --project=examples examples/ads_kepler.jl
```

| Script | Calculation |
|:--|:--|
| `sine.jl` | Twentieth-order expansion and evaluation of sine |
| `polynomial_inversion.jl` | Tenth-order inverse of sine, compared with arcsine |
| `tutorial1ex8.jl` | First-order gradient of the sombrero function about (2, 3) |
| `ode_integration.jl` | Second-order Kepler flow and its state transition matrix |
| `damc_kepler.jl` | 10,000 Monte Carlo samples after 30 revolutions, compared with orders 2, 4 and 8 |
| `ads_kepler.jl` | Adaptive Taylor patches for an uncertain Kepler map, with independent error checks |

The Kepler examples use fixed random seeds and write PDF and PNG figures in the
working directory. CairoMakie renders without an interactive display. The
documentation build executes all six scripts and includes the resulting figures.

The [benchmarks](../benchmark) compare numerical integration with TaylorSeries.jl
and measure package loading and first-use latency.
