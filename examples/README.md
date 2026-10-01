# Examples

The first five scripts reproduce every example in
[DACE.jl](https://github.com/UoA-AstroGroup/DACE.jl/tree/7271c372ce86350e80a95a9c2c0d7e8739d09e83/examples)
using DifferentialAlgebra.jl's native Julia API. Each script is self-contained
and includes numerical checks and CairoMakie plots. An additional Kepler example
demonstrates automatic domain splitting (ADS). All six scripts become executable
documentation pages through Literate.jl.

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
| `sine.jl` | Twentieth-order expansion of sine, with approximation and error curves |
| `polynomial_inversion.jl` | Tenth-order inverse of sine, with convergence toward arcsine |
| `tutorial1ex8.jl` | Sombrero gradient about (2, 3), with contours and a tangent section |
| `ode_integration.jl` | Kepler trajectory, state transition matrix, and nonlinear flow-map errors |
| `damc_kepler.jl` | 10,000 Monte Carlo samples after 30 revolutions, with maps of orders 2, 4 and 8 |
| `ads_kepler.jl` | Kepler ADS partition, mapped uncertainty, and error heatmaps |

Every script writes PDF and PNG figures under `figures/` in the working directory.
The Monte Carlo and ADS examples use fixed random seeds. CairoMakie renders
without an interactive display. The documentation build executes all six scripts
and includes the figures with PDF download links. Generated figures are ignored
by Git.
