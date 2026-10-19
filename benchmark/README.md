# Benchmarks and verification

Run from the repository root. The benchmark environment contains only the package,
interval/AD support, and the ODE interfaces needed by these scripts; it does not
require the plotting or documentation environment.

```sh
julia --startup-file=no --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --startup-file=no --project=benchmark benchmark/numerical_kernels.jl
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

For the numerical-kernel baseline, create a detached worktree at `720b46c` and
copy the root project's local `Manifest.toml` into it. Run
`benchmark/numerical_kernels.jl` first with `--project=<baseline-worktree>` and
then with `--project=.`. This script needs only the root package dependencies;
use the same Julia version and manifest for both runs. The broader ADS scripts
use the benchmark environment above. Local manifests are not committed.

## Numerical cleanup from `720b46c`

The paired script measures analytic Horner/recurrence paths, quotient construction,
repeated differentiation, Hessians, coefficient norms and symmetric eigenpairs.
It prints the formerly lost coefficients alongside timings, so an incorrectly
shortened baseline result is visible. Global tolerance is configured outside the
measured calculations. See the [migration note](../docs/src/migration.md).

Paired runs compare `720b46c` with `1a66559` on Windows, an Intel i7-12700,
Julia 1.13.1 and SpecialFunctions 2.9.0, with identical dependency manifests.
BigFloat measurements use the default 256-bit precision. The numerical source is
unchanged by the results-only commit that records this report.

The following paths use the zero-byte multiplication budget. Times are warmed
minima and bytes are per call; baseline rows with discarded coefficients do less
mathematical work and should not be interpreted as equivalent correct results.

| Path | Before → after (μs) | Before → after (B) |
|:--|--:|--:|
| `exp`, three variables/order six | 0.91 → 0.98 | 1,888 → 1,888 |
| `log` | 0.93 → 1.29 | 2,064 → 2,064 |
| Polynomial quotient | 1.27 → 1.99 | 2,976 → 2,976 |
| `erf` | 0.97 → 0.95 | 2,000 → 2,000 |
| Hessian, three variables/order six | 1.93 → 1.94 | 10,784 → 10,784 |
| Scaled quotient, tolerance `1e-6` | 0.12 → 0.16 | 672 → 704 |
| Scaled order-ten `exp`, tolerance `1e-4` | 0.25 → 0.29 | 608 → 608 |
| Scaled order-two `erf`, tolerance `1e-4` | 0.11 → 0.09 | 496 → 496 |
| Repeated derivative, stored four/working two | 0.05 → 0.08 | 384 → 576 |
| Hessian, stored four/working two | 0.08 → 0.08 | 528 → 528 |

Norms and eigenpairs:

| Path | Before → after (μs) | Before → after (B) |
|:--|--:|--:|
| Float32 coefficient norm | 0.12 → 0.03 | 0 → 0 |
| Float32 degree norms | 0.14 → 0.12 | 96 → 160 |
| Float64 coefficient norm | 0.15 → 0.04 | 0 → 0 |
| Float64 degree norms | 0.14 → 0.14 | 112 → 192 |
| BigFloat coefficient norm | 8.86 → 2.15 | 7,208 → 7,112 |
| BigFloat degree norms | 8.81 → 11.77 | 8,080 → 8,520 |
| Symmetric eigenpair lifting | 5.29 → 5.57 | 19,136 → 18,720 |
| Diagonal eigenpairs | 0.30 → 0.63 | 1,840 → 2,896 |
| Empty integer eigenpairs | 0.01 → 0.01 | 80 → 80 |

Stable degree norms use a second scaled accumulation buffer. Repeated derivatives
retain one additional intermediate result in the stored-order regression. Diagonal
eigenpairs now consult the scalar solver for result types, adding work relative
to the former copy-only path. The empty path derives that working type directly,
so it remains usable for BigFloat even when an empty scalar solve is unsupported.
Submicrosecond timings approach timer resolution; allocation counts are useful
alongside them. No solver or polygon allocation tuning is included.

The 32 MiB multiplication budget is also measured in the raw results. The quotient
constant/linear terms are now `1`/`-0.01` with either budget. The tenth exponential
coefficient is about `1.85679756368096e11`; the quadratic error-function coefficient
is `-0.005079293868746624`. Both derivative regressions now return coefficient `12`
of `x²`. Previously the fallback quotient/exponential, both error-function paths,
and both derivative paths returned zero for these requested coefficients. Empty
integer eigenpairs now use Float64 coefficient containers.

Compilation below is measured inside each first call, excluding package loading
and the measurement wrapper. Later calls can reuse earlier specializations; these
are individual observations rather than complete package time-to-first-execution.

| First invocation | Compilation before → after (ms) |
|:--|--:|
| `exp` | 2.34 → 2.41 |
| `log` | 51.00 → 58.79 |
| `erf` | 40.48 → 50.21 |
| Hessian | 94.05 → 70.51 |
| Float32 coefficient norm | 88.50 → 174.18 |
| BigFloat coefficient norm | 139.33 → 183.25 |
| BigFloat degree norms | 49.38 → 48.29 |
| Eigenpair lifting | 1225.39 → 1194.58 |
| Empty integer eigenpairs | 556.10 → 637.36 |

All first-call runtimes and both budgets are retained in the
[baseline measurements](measurements/720b46c-numerical_kernels.txt) and
[cleanup measurements](measurements/1a66559-numerical_kernels.txt).
The only new numerical helper normalizes analytic arguments without filtering;
existing multiplication/scaling/calculus kernels receive local order/filter
controls. Coefficient norms reuse LinearAlgebra's stable implementation.

The isolated baseline passed **9,984/9,984** tests on Julia 1.13.1. The final
cleanup passed **9,991/9,991** on **both Julia 1.10.12 and 1.13.1**, with bounds
checking and Aqua enabled. Aqua retains only the five exact constructor-method
exceptions documented in `test/quality.jl`. Concurrent shutdown probes initially
exceeded Aqua's 30-second limit; staggered probes passed with that default unchanged.
Literate's catalog/rendering checks passed **10/10**, all **30** tutorials executed,
and Documenter's rendered catalog/link/asset checks passed **3,553/3,553** for flat
URLs and **3,269/3,269** for directory URLs. The asset audit excludes only the two
version-selector files created by deployment. Runic 1.11.1 and `git diff --check`
also passed. The benchmark script completed on Julia 1.10.12 as well; its
compilation statistic is `missing`, as documented above.

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
