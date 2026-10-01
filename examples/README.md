# Examples

Run a standalone tutorial from the repository root:

```sh
julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=examples examples/ads_methods.jl
```

For interactive figures, start Julia with `julia --project=examples` and use
`include("examples/ads_methods.jl")`. Each figure is returned as the last
expression of its plotting block.

The [tutorial catalog](../docs/src/examples.md) organizes the scripts by topic.
The [documentation](https://jackyarndley.github.io/DifferentialAlgebra.jl/examples/)
includes their executed numerical output and figures.
