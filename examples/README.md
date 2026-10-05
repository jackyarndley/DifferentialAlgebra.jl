# Examples

Run a standalone tutorial from the repository root:

```sh
julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=examples examples/ads_methods.jl
```

For interactive figures, start Julia with `julia --project=examples` and use
`include("examples/ads_methods.jl")`. Each figure is returned as the last
expression of its plotting block.

The [tutorial catalog](../docs/src/examples.md) and documentation sidebar group
the scripts into polynomial arithmetic/calculus, implicit equations, time
integration, orbit propagation, ADS methods/geometry, continuous ADS/optimization,
validated intervals, orbit determination, and storage.
`ads_methods.jl` compares the heuristic methods with interval-based splitting
and oriented geometry. `ads_kepler.jl` adds oriented fits and C2 blending;
`ads_flow.jl` compares three unvalidated flow-splitting strategies.
Interval bands cover complete cells; sampled error plots are numerical checks.
The [documentation](https://jackyarndley.github.io/DifferentialAlgebra.jl/examples/)
includes their executed numerical output and figures.
