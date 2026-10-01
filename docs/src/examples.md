# Examples

These executable tutorials introduce Taylor arithmetic, inverse maps, derivatives,
and uncertainty propagation. Every page is generated from a self-contained Julia
script with Literate.jl. CairoMakie renders the figures during the documentation
build; each page includes a downloadable PDF.

| Example | What it demonstrates |
|:--|:--|
| [Sine function](generated/sine.md) | How polynomial order changes the approximation and its error |
| [Polynomial inversion](generated/polynomial_inversion.md) | Constructing a local inverse and understanding its range of accuracy |
| [Sombrero gradient](generated/tutorial1ex8.md) | Extracting derivatives and interpreting the tangent approximation |
| [Orbit integration](generated/ode_integration.md) | A Kepler flow map, its state transition matrix, and nonlinear corrections |
| [Kepler Monte Carlo](generated/damc_kepler.md) | Propagating 10,000 uncertain states with polynomial maps |
| [Kepler domain splitting](generated/ads_kepler.md) | Building an adaptive map, inspecting its partition, and validating its error |

The first five scripts adapt every example in
[DACE.jl](https://github.com/UoA-AstroGroup/DACE.jl/tree/7271c372ce86350e80a95a9c2c0d7e8739d09e83/examples)
to the DifferentialAlgebra API. Automatic domain splitting is an additional example.

## Run locally

From the repository root, install the example environment once:

```sh
julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```

Then run any script, for example:

```sh
julia --project=examples examples/ads_kepler.jl
```

PNG and PDF figures are written to `figures/` in the working directory. No
interactive display is required. The Monte Carlo and ADS examples use fixed
random seeds, and every example checks its results against an independent formula
or numerical calculation.
