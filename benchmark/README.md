# Benchmarks and verification

Run from the repository root. The benchmark environment contains only the package,
interval/AD support, and the ODE interfaces needed by these scripts; it does not
require the plotting or documentation environment.

```sh
julia --startup-file=no --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --startup-file=no --project=benchmark benchmark/cleanup.jl
julia --startup-file=no --project=benchmark benchmark/unified_ads.jl
julia --startup-file=no --project=benchmark benchmark/interval_models.jl
julia --startup-file=no --project=benchmark benchmark/polygon_ads.jl
julia --startup-file=no --project=benchmark benchmark/continuous_ads.jl
```

`measurement.jl` reports the minimum of 30 warmed batches of ten calls, allocated
bytes per call, first-call runtime, and compilation time inside that first call.
Package loading and compilation of the measurement wrapper are excluded. First
calls later in a process can reuse earlier specializations; zero compilation does
not mean that a path has no compilation cost. Scripts report both multiplication
budgets (zero and 32 MiB), numerical work, and interval widths where applicable.
Dynamic first-call entry keeps workload compilation inside the timer; a separate
specialized loop measures warmed calls without runtime dispatch. Julia 1.10
reports first-call elapsed time but returns `missing` for the compilation statistic.
Run paired revisions sequentially without competing Julia jobs. Timings are
machine dependent; compare allocations, patch counts, and numerical contracts too.

For a baseline, create a detached worktree with `git worktree add --detach
../da-baseline 17f74cc`. Copy `benchmark/Project.toml` and its locally generated
`Manifest.toml` into a temporary environment, then run `Pkg.develop(path="../da-baseline")`
there. Run the current scripts with that environment to measure identical workloads
and dependency versions against the old source. Local manifests are not committed.

## Cleanup from `17f74cc`

The cleanup retains the shared ADS driver and geometry/estimator adapters, both
multiplication kernels, BigFloat special functions, compiled snapshot ownership,
and validated rounding/remainder contracts. The new helpers share only checked
multi-index validation, a scalar dense Horner evaluator, and benchmark setup.
See the [API migration note](../docs/src/migration.md) for changed entry points.

Paired runs compare `17f74cc` with source revision `e93621a` on Windows, Julia
1.13.1, IntervalArithmetic 1.0.12, and OrdinaryDiffEqCore 4.19.0. Both environments
used the same dependency versions. The table uses the 32 MiB budget; both budgets
are retained in the raw results. Times are warmed minima in microseconds, and
bytes are allocated per call.

| Path | Before → after (μs) | Before → after (B) |
|:--|--:|--:|
| Six-variable Hessian, order four | 177.11 → 74.64 | 143,232 → 83,824 |
| `tan(::TimeSeries)`, order twenty | 1.51 → 0.81 | 1,120 → 672 |
| `tanh(::TimeSeries)` | 1.50 → 0.81 | 1,120 → 672 |
| Recorded `tan` ODE expansion | 1.98 → 1.53 | 3,056 → 2,592 |
| Recorded `tanh` ODE expansion | 1.96 → 1.50 | 3,056 → 2,592 |
| Automatic ordinary polygon frame/construction | 777.06 → 653.99 | 459,624 → 455,000 |
| Automatic validated polygon frame/construction | 1,031.19 → 909.92 | 481,304 → 473,560 |
| Mutating dense interpolation, two states/order twenty | 0.09 → 0.05 | 80 → 0 |

Both polygon cases use `exp(x+y)` on `[-1/2,1/2]²`, order three, tolerance `0.001`,
and four patches. Without multiplication tables, Hessian time was 90.71 → 39.77 μs;
validated polygon time was essentially unchanged (919.29 → 921.13 μs), with fewer
allocations. Very short timings approach timer resolution; allocation figures
provide a clearer comparison for dense interpolation.

First-call compiler costs below are milliseconds from each path's first invocation
(the zero-budget pass for paths inside the budget loop). The later table-enabled
pass reuses these specializations. Dense interpolation is measured after both passes. These costs
include workload/callback compilation, excluding package loading and compilation
of the measurement wrapper. They are individual first-call observations, not warmed
minima or complete package time-to-first-execution measurements.

| Path | Compilation before → after (ms) |
|:--|--:|
| Hessian | 32.84 → 90.89 |
| TimeSeries `tan` / `tanh` | 77.26 → 43.15 / 5.83 → 5.10 |
| Recorded `tan` / `tanh` expansion | 544.39 → 509.74 / 158.17 → 180.91 |
| Ordinary / validated polygons | 5,954.94 → 6,016.73 / 5,075.68 → 4,705.69 |
| Dense interpolation | 92.02 → 26.10 |

Raw paired results:

| Script | `17f74cc` | `e93621a` |
|:--|:--|:--|
| Affected paths | [before](measurements/17f74cc-cleanup.txt) | [after](measurements/e93621a-cleanup.txt) |
| Shared ADS | [before](measurements/17f74cc-unified_ads.txt) | [after](measurements/e93621a-unified_ads.txt) |
| Interval models | [before](measurements/17f74cc-interval_models.txt) | [after](measurements/e93621a-interval_models.txt) |
| Polygons | [before](measurements/17f74cc-polygon_ads.txt) | [after](measurements/e93621a-polygon_ads.txt) |
| Continuous ADS | [before](measurements/17f74cc-continuous_ads.txt) | [after](measurements/e93621a-continuous_ads.txt) |

All four broader scripts completed with both multiplication budgets. After
excluding timing, allocation, and snapshot-storage fields, their printed patch
counts, estimates, rigorous bounds, and enclosure widths matched exactly.

The isolated baseline passed **9,602/9,602** tests on Julia 1.13.1. The cleanup
passed **9,984/9,984** on **both Julia 1.10.12 and 1.13.1**, with bounds checking
enabled. Aqua's only exclusions are the five exact generated outer-constructor
methods documented in `test/quality.jl`, for empty-tuple type inference on Julia
1.10. Literate's catalog/figure checks passed 10/10; all **30** standalone
tutorials executed, the complete Documenter build passed, and its rendered
navigation/link/figure checks passed **196/196**. Runic (including benchmarks),
the benchmark helper on Julia 1.10, and `git diff --check` passed.

## Verification commands

```sh
julia --startup-file=no --project -e 'using Pkg; Pkg.test(; julia_args=["--startup-file=no", "--check-bounds=yes"])'
julia +1.10 --startup-file=no --project -e 'using Pkg; Pkg.test(; julia_args=["--startup-file=no", "--check-bounds=yes"])'
julia --startup-file=no --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --startup-file=no --project=docs docs/test_literate.jl
julia --startup-file=no --project=docs docs/make.jl
julia --startup-file=no --project=@runic -e 'using Runic; exit(Runic.main(ARGS))' -- --check --diff --docstrings src test ext examples benchmark docs/make.jl docs/literate.jl docs/test_literate.jl docs/test_navigation.jl
git diff --check
```

The Ubuntu CI matrix remains Julia 1.10 and latest stable. Local runs on Windows
do not substitute for that CI job. Namespace isolation uses SafeTestsets, while
context-mutating files remain sequential. Independent coefficient/moment oracles,
interval guarantees and metadata, precision, aliasing, both multiplication tables,
loading order, catalog registration, links, and generated figures remain tested.
