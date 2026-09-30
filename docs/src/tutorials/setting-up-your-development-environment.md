# Development

Use Julia 1.10 or later. From the repository checkout:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.test()
```

The Julia polynomial engine requires no DACE binary artifact or compiler.
Special functions use SpecialFunctions.jl and QuadGK.jl.

## Documentation and examples

From the repository root:

```julia
using Pkg
Pkg.activate("docs")
Pkg.develop(path=".")
Pkg.instantiate()
include("docs/make.jl")
```

Open `docs/build/index.html` to read the result. The build executes the sine,
gradient, inversion and ODE examples. CI runs tests on Linux with Julia 1.10 and
the latest stable Julia release, and builds documentation with Julia 1.10.
Documentation deployment is not configured.

Examples and benchmarks have separate environments; see
[examples](https://github.com/jackyarndley/DifferentialAlgebra.jl/tree/main/examples)
and [benchmark](https://github.com/jackyarndley/DifferentialAlgebra.jl/tree/main/benchmark)
for setup commands. Manifests, generated documentation and local results are
excluded from version control.
