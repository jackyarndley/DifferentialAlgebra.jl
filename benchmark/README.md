# Benchmarks

Prepare the separate environment from the repository root:

```sh
julia --startup-file=no --project=benchmark -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=benchmark benchmark/review_kernels.jl
```

Run each script in a fresh process. Algebra initialization invalidates earlier
live polynomials, so preconstructing all workloads as one shared suite is unsafe.
Plotting, Documenter and GenericLinearAlgebra are not required.

## Choose a workload

| Script | Measures | Recorded output |
|:--|:--|:--|
| `compound_kernels.jl` | Compound polynomial arithmetic | [Log](measurements/compound_kernels.txt) |
| `numerical_kernels.jl` | Numerical kernel costs | [Log](measurements/numerical_kernels.txt) |
| `review_kernels.jl` | Hypot/Bessel, Duffing, eigenpairs, policy/model/geometry paths | [Log](measurements/review_kernels.txt) |
| `cleanup.jl` | Calculus, time series and solver buffers | [Log](measurements/cleanup.txt) |
| `unified_ads.jl` | Shared construction/refinement driver | [Log](measurements/unified_ads.txt) |
| `interval_models.jl` | Model arithmetic, ranges and six-coordinate uncertainty | [Log](measurements/interval_models.txt) |
| `interval_kernels.jl` | Coordinate powers, discarded products and both table budgets | [Log](measurements/interval_kernels.txt) |
| `interval_cleanup.jl` | Scalar model operations and live/snapshot queries | [Log](measurements/interval_cleanup.txt) |
| `polygon_ads.jl` | Box versus exact polygon geometry | [Log](measurements/polygon_ads.txt) |
| `continuous_ads.jl` | Blending, AD derivatives and enclosure queries | [Log](measurements/continuous_ads.txt) |

`interval_cleanup.jl --bigfloat` adds 256-bit endpoints.
`interval_kernels.jl --profile` collects CPU/allocation samples. Logs include
coefficients, bounds/metadata, patch counts or conservation checks alongside costs.
Their headers identify the measured revision; they are not claims about an
unmeasured current checkout.

## Interpret measurements

All scripts use `measurement.jl`, which records revision, Julia/dependency versions,
CPU and threads. It reports warmed minimum runtime, bytes, allocations, batch size,
first-call elapsed time and compilation. Defaults request 30 samples and a
0.5-second trial budget; `measure` accepts `samples`, `seconds` and `batch`.
Use `batch=1` for nonrepeatable mutations; existing buffer workloads overwrite
outputs on every call.

First-call measurements exclude package loading and harness compilation. Later
workloads reuse earlier specializations, so they describe script order rather
than isolated compilation costs. Run paired revisions sequentially with matching
hardware, dependencies and harness. Compare numerical results and allocations
before interpreting a timing difference; an inaccurate baseline may do less work.

Interval widths describe their stated target: a stored-polynomial range differs
from an original-function enclosure. Fewer polygon patches do not necessarily
mean less time because exact geometry costs more. Separate construction,
evaluation and derivative costs when choosing a method.

Keep one current log per script in `measurements/`, with older results in Git
history and paired comparisons in the PR description. Portable CI checks
numerical contracts and formatting; controlled local runs measure performance.
