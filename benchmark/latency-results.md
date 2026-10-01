# Package loading and first-use latency

The package precompiles one small workload with ordinary vectors: arithmetic,
elementary functions, differentiation and compiled-map evaluation. It uses
Float64 coefficients at order four in two variables. Order and variable count
are runtime data, so this caches methods useful at other orders and dimensions
without creating a family of fixed-size coefficient types.

The measured workload uses order five and constructs
`exp(x/10)*sin(y/10) + sqrt(2+x-y)/(3+x+y)`, compiles it together with its first
partial derivative, and evaluates the map at vector coordinates. The package
workload does not load solvers, plotting libraries or ADS callbacks.

## Results

Measured on 2026-10-01 with Julia 1.13.1, Linux/WSL, an Intel Core i7-12700,
and one compute thread pinned to logical CPU 2. Both runs used identical source
(the hash is recorded in the raw files) and dependency versions. Only
DifferentialAlgebra's `precompile_workload` preference changed. Dependencies
were already installed and precompiled. Each entry is a median of three fresh
processes; the disabled batch preceded the enabled batch.

| Measurement | Workload disabled | Workload enabled |
|:--|--:|--:|
| Rebuild package cache | 1.481 s | 3.568 s |
| Package cache size (`.ji` + native image) | 1.102 MiB | 2.579 MiB |
| Import DifferentialAlgebra | 0.171 s | 1.796 s |
| First polynomial workload, ordinary coordinates | 1.472 s | 0.170 s |
| Import plus first polynomial workload | 1.652 s | 1.967 s |

The polynomial workload's first call is **8.64× faster**, but importing the
enabled cache was slower in this run. The combined import and first-call time
therefore increased; these measurements do not demonstrate an overall first-use
latency improvement.

The workload adds approximately 2.09 s of package-cache construction and
1.48 MiB of cache. Import timings depend on the environment and file-system
state. All timings exclude Julia process startup and first-time installation or
compilation of dependencies.

Every process checked that loading the package left its global algebra
uninitialized. The same first-use workload was run twice to
check its output; warm-call timings are recorded for completeness, not as a
microbenchmark. First-call allocation counts include compilation.

Application callbacks, other coefficient types and optional packages still
compile when first used. In particular, this workload does not eliminate the
compilation cost of a new ODE solve, an ADS callback or a plotting environment.

Raw measurements: [disabled](results-latency-before.toml),
[enabled](results-latency.toml).

## Reproduce

From the repository root, instantiate the benchmark environment and run:

```sh
julia --startup-file=no benchmark/setup.jl
julia --startup-file=no --project=benchmark benchmark/latency.jl latency.toml
```

The script rebuilds the package cache explicitly with `Base.compilecache`,
then measures imports and first use in separate processes. Imports are timed
at top level. Workload calls use `invokelatest` so compilation happens inside
the timer instead of during compilation of the measurement function.

For the disabled comparison, add the following preference to
`benchmark/LocalPreferences.toml`, preserving any existing settings:

```toml
[DifferentialAlgebra]
precompile_workload = false
```

Remove that setting to restore the default enabled workload, then rerun with
a different output filename. The script records the project preference,
source hash, platform, timings and cache-file sizes. The cache build measures
DifferentialAlgebra itself with cached dependencies, not the cost of an empty
Julia depot.
