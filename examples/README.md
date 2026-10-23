# Examples

Each Literate tutorial runs directly from the checkout, with numerical checks
and figures. Shared physical models live in `support/` so companion tutorials
use the same equations. Start with `implicit_maps.jl` for nonlinear sensitivity,
`ode_integration.jl` for a six-state flow map, or `interval_models.jl` for rigorous
enclosures. `daiod.jl` combines implicit maps, uncertainty and a measurement update.

```sh
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
julia --project=examples examples/implicit_maps.jl
```

For interactive figures, start `julia --project=examples` and run
`include("examples/implicit_maps.jl")`. Each plotting block returns its figure;
terminal execution prints reports. Use `CairoMakie.save("figure.png", fig)` to
export a figure.

The [catalog](../docs/src/examples.md) lists all 30 scripts by question and result.
The [rendered tutorials](https://jackyarndley.github.io/DifferentialAlgebra.jl/examples/)
include their executed reports and plots. Sampled discrepancies and numerical
integration checks are labeled separately from interval-certified bounds.

Examples carrying Apache-2.0 headers retain their [third-party notices](../THIRD_PARTY_NOTICES.md).
