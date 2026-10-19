# Benchmarks and verification

Run from the repository root. The benchmark environment contains the package,
interval/AD support, and the ODE interfaces needed by the scripts. It does not
require the plotting or documentation environment.

```sh
julia --startup-file=no --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --startup-file=no --project=benchmark benchmark/compound_kernels.jl
julia --startup-file=no --project=benchmark benchmark/numerical_kernels.jl
julia --startup-file=no --project=benchmark benchmark/cleanup.jl
julia --startup-file=no --project=benchmark benchmark/unified_ads.jl
julia --startup-file=no --project=benchmark benchmark/interval_models.jl
julia --startup-file=no --project=benchmark benchmark/polygon_ads.jl
julia --startup-file=no --project=benchmark benchmark/continuous_ads.jl
```

`measurement.jl` reports the minimum of 30 warmed batches of ten calls, allocated
bytes per call, first-call runtime, and compilation inside that first call.
Package loading and compilation of the measurement wrapper are excluded. Later
calls can reuse earlier specializations; zero compilation does not mean that a
path has no compilation cost. Julia 1.10 reports first-call elapsed time and
`missing` for the separate compilation statistic. Run paired revisions
sequentially without competing Julia jobs. Submicrosecond minima are sensitive
to timer resolution and CPU frequency; allocations and numerical results also
matter.

For the compound-kernel comparison, create a detached worktree at `24dc1a1`.
Use the current `compound_kernels.jl` script for both revisions, changing only
the benchmark environment's DifferentialAlgebra path with `Pkg.develop`.
Keep dependency versions and source hashes identical within each Julia pair.
The script prints the loaded package revision, multiplication budget, retained
regression coefficients, interval bounds, and metadata. Configuration occurs
outside the timed calls. It never loads GenericLinearAlgebra.

## Compound kernels from `24dc1a1` to `1efc220`

Measurements used Windows, an Intel i7-12700 (`alderlake`), one Julia thread,
IntervalArithmetic 1.0.12, SpecialFunctions 2.9.0, and default 256-bit BigFloat
precision. Compound rows use one variable/order four; interval rows use
`Interval{Float64}`, two variables/order four; diagonal matrices are 2×2.
Both Julia versions used `--startup-file=no --check-bounds=yes`.
Dependency versions and source hashes matched within each pair. The measured
numerical source is unchanged by the results-only commit recording this report.
See the [API migration note](../docs/src/migration.md).

Julia 1.13.2, zero-byte multiplication budget:

| Path | Before → after (μs) | Before → after (B) |
|:--|--:|--:|
| Integer cube | 0.26 → 0.18 | 432 → 352 |
| Negative integer power | 0.36 → 0.42 | 1,088 → 1,088 |
| `log2` | 0.23 → 0.51 | 672 → 704 |
| `hypot` | 0.34 → 0.99 | 896 → 1,152 |
| `hypot` with a small scalar | 1.20 → 1.10 | 1,312 → 1,216 |
| `atan` with a small scalar | 0.50 → 0.50 | 992 → 1,376 |
| `gamma` | 0.55 → 1.04 | 1,056 → 1,072 |
| Polynomial exponent | 0.44 → 1.00 | 1,152 → 1,184 |
| Stored order four/working order two power | 0.22 → 0.26 | 432 → 352 |
| Point-interval coefficient norm | 0.29 → 2.96 | 0 → 1,376 |
| Positive-width interval degree norm | 0.90 → 4.17 | 352 → 1,544 |
| Float64 diagonal eigenpairs | 1.73 → 0.42 | 2,912 → 1,312 |

Correct interval norms cost additional storage and arithmetic. The shared
interval accumulator groups coefficients by degree or into one group for the
full norm, using guaranteed point scales and IA's integer powers/roots. The
ordinary-number norm implementations are unchanged. Compound operations reuse
private entry points to the existing multiplication, scaling and analytic
kernels, with local filtering flags; no global tolerance or cutoff changes.

The positive-width degree norm previously returned about `[4, 8.004]` and lost
its guarantee flag. It now encloses `[5, sqrt(41)]` to outward-rounding accuracy,
with `com` decoration and its guarantee intact. The point degree-norm case now
returns the exact point interval `[5, 5]`. Zero-containing intervals previously threw; the
new full/degree norm calls take 0.35/0.47 μs and allocate 912/1,040 B.
BigFloat/Rational{BigInt} diagonal solves previously threw without an optional
solver; they now take 2.19/6.72 μs and allocate 2,552/5,096 B.

With either multiplication budget, the `log2` linear coefficient changes from
zero to `0.01298425536800067`, and the cube's quadratic coefficient changes
from `0.162` to `0.243`. Incorrect baseline results do less mathematical work;
the timings do not compare equivalent correct calculations. At 32 MiB, the
integer cube takes 0.29 → 0.12 μs, `log2` 0.13 → 0.08 μs, and `hypot`
0.46 → 0.23 μs; their allocations are 432 → 352, 224 → 256, and 704 → 704 B.

Julia 1.10.12, zero-byte multiplication budget:

| Path | Before → after (μs) | Before → after (B) |
|:--|--:|--:|
| Integer cube | 0.32 → 0.30 | 432 → 352 |
| `hypot` | 0.27 → 0.49 | 896 → 1,152 |
| Stored order four/working order two power | 0.28 → 0.14 | 432 → 352 |
| Positive-width interval degree norm | 0.54 → 2.41 | 352 → 1,832 |
| Float64 diagonal eigenpairs | 1.14 → 0.36 | 2,848 → 1,200 |

Compilation inside first invocations on Julia 1.13.2:

| Path | Before → after (ms) |
|:--|--:|
| Integer cube | 12.11 → 2.98 |
| `log2` | 6.07 → 3.77 |
| `gamma` | 2.95 → 3.97 |
| Point-interval coefficient norm | 8.40 → 7.61 |
| Float64 diagonal eigenpairs | 4.75 → 6.69 |
| Point-interval degree norms | 3.63 → 4.05 |

Julia 1.10 first-call elapsed times for the cube and diagonal eigenpairs are
3.39 → 3.80 ms and 4.39 → 4.66 ms. These observations include reuse from package
precompilation, input setup, and earlier paths in the same script; they are not
whole-package compilation totals. Full measurements retain first-call times and both budgets:

- [Julia 1.13.2 baseline](measurements/24dc1a1-compound_kernels-julia1.13.2.txt)
- [Julia 1.13.2 follow-up](measurements/1efc220-compound_kernels-julia1.13.2.txt)
- [Julia 1.10.12 baseline](measurements/24dc1a1-compound_kernels-julia1.10.12.txt)
- [Julia 1.10.12 follow-up](measurements/1efc220-compound_kernels-julia1.10.12.txt)

Earlier numerical-kernel measurements remain revision-tagged in
[720b46c](measurements/720b46c-numerical_kernels.txt) and
[1a66559](measurements/1a66559-numerical_kernels.txt).

## Verification commands

```sh
julia --startup-file=no --project -e 'using Pkg; Pkg.test(; julia_args=["--startup-file=no", "--check-bounds=yes"])'
julia +1.10 --startup-file=no -e 'using Pkg; Pkg.activate(; temp=true); Pkg.develop(path=pwd()); Pkg.test("DifferentialAlgebra"; julia_args=["--startup-file=no", "--check-bounds=yes"])'
git diff --check
```

Use separate environments for each Julia version to avoid rewriting a shared
manifest. See [Contributing](../docs/src/contributing.md) for Runic, documentation,
and Literate commands. The Ubuntu CI matrix remains Julia 1.10 and latest stable;
local Windows checks do not substitute for that CI job. Context-mutating tests
remain sequential, and Aqua keeps its five exact constructor-method exceptions.
No solver or polygon allocation tuning is included in this patch.

## Verification at `1efc220`

`Pkg.test()` passed **10,464/10,464** checks across 31 sequential SafeTestsets on
**Julia 1.10.12 and 1.13.2**, with bounds checking and Aqua enabled. Fresh-process
checks covered optional extension loading and BigFloat/Rational{BigInt}
polynomial diagonal eigenpairs without GenericLinearAlgebra. The suite also
checked snapshot ownership.
Runic 1.11.1 passed on both Julia versions; all added Julia lines fit 92 columns.
Literate's catalog/rendering checks passed **10/10**, and the audit of the existing
rendered documentation passed **3,553/3,553** catalog/link/asset checks using the
Julia 1.13.1 documentation environment. A fresh full documentation build was not
run for this follow-up. No failing checks remain. IntervalArithmetic's existing
invalid-input and rounding-probe warnings still appear in the test output.
