# Benchmarks

Run from the repository root in the small benchmark environment. It includes
BenchmarkTools and the optional interval, AD and ODE interfaces used below.
Plotting, Documenter and GenericLinearAlgebra are not required.

```sh
julia --startup-file=no --project=benchmark -e 'using Pkg; Pkg.develop(path=".")'
julia --startup-file=no --project=benchmark -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=benchmark benchmark/compound_kernels.jl
julia --startup-file=no --project=benchmark benchmark/numerical_kernels.jl
julia --startup-file=no --project=benchmark benchmark/cleanup.jl
julia --startup-file=no --project=benchmark benchmark/unified_ads.jl
julia --startup-file=no --project=benchmark benchmark/interval_models.jl
julia --startup-file=no --project=benchmark benchmark/interval_kernels.jl
julia --startup-file=no --project=benchmark benchmark/polygon_ads.jl
julia --startup-file=no --project=benchmark benchmark/continuous_ads.jl
```

Use separate environments for Julia 1.10 and latest stable. The setup uses
`Pkg.develop`, including on Julia 1.10. Shared compatibility entries follow the
package's integration bounds; see [Contributing](../docs/src/contributing.md).

All scripts use `measurement.jl`. Each log starts with the loaded package's full
revision, Julia version, CPU, thread count and loaded dependency versions.
Workloads also print regression coefficients, interval bounds, patch counts or
other numerical diagnostics alongside their measurements.

The harness reports warmed minimum runtime, bytes, allocation count, evaluations
per sample, first-call elapsed time and compilation within that first call.
BenchmarkTools calibrates the sample batch for each workload. The defaults request
30 samples with a 0.5-second trial budget; `samples`, `seconds` and an explicit
`batch` can be supplied to `measure`. Use `batch=1` for nonrepeatable mutations.
Existing buffer workloads overwrite their output on every call.

First-call compilation is measured before tuning the warmed loop. Package loading
and harness compilation are excluded. Later workloads reuse earlier specializations,
so compilation rows describe their script order rather than isolated function costs.
Julia 1.10 reports first-call elapsed time and `missing` for the compilation field.

Run paired revisions sequentially on the same machine with matching dependencies
and the same script/harness source. Compare allocations and numerical results as
well as timings. Baselines with incorrect coefficients do less mathematical work;
their timings describe a correctness tradeoff.

Keep only the latest output per script in `measurements/`, with its revision in the
header. Historical logs remain in Git history. Release changes and API migration
instructions live in [NEWS](../NEWS.md) and the [migration guide](../docs/src/migration.md).

Current Julia 1.13.2 outputs (the header identifies the measured revision):

| Script | Output |
|:--|:--|
| Compound kernels | [Measurements](measurements/compound_kernels.txt) |
| Numerical kernels | [Measurements](measurements/numerical_kernels.txt) |
| Calculus, time series and solver buffers | [Measurements](measurements/cleanup.txt) |
| Shared ADS driver | [Measurements](measurements/unified_ads.txt) |
| Interval models | [Measurements](measurements/interval_models.txt) |
| Interval kernel comparison | [Measurements](measurements/interval_kernels.txt) |
| Polygon ADS | [Measurements](measurements/polygon_ads.txt) |
| Continuous ADS | [Measurements](measurements/continuous_ads.txt) |

The scripts retain separate algebra setup and numerical diagnostics. Run them in
fresh processes: preconstructing a shared suite across initializations would leave
earlier polynomial fixtures invalid. Performance comparisons use a controlled local
machine; portable CI checks numerical contracts and formatting.

The interval kernel experiment compares normalized-coordinate powers with direct
`IntervalArithmetic.pown`, and discarded products with the decorated fallback.
It also measures dense model multiplication and the uncertainty expression under
both multiplication-table budgets. Pass `--profile` to collect CPU samples and
allocation stacks. Run paired revisions sequentially using the same script; the
reference revision can omit the new private discarded-product comparison.
Full-precision bounds, decorations and guarantee flags accompany each workload.
See [kernel results](interval_kernels.md) for the measured scope and comparison.
