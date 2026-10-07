# Interval models and ADS

## Unified ADS refactor from `c4b256b`

The checkout was clean `main` at `c4b256b4d397c4fba916568905b0e64b0a0840c6`.
There were no repository or parent `AGENTS.md` instructions. An isolated worktree
of that revision passed all **8,853** tests. The initial test attempt reported
one continuity failure; neither its isolated focused rerun nor the isolated
complete baseline reproduced it. No reproducible baseline failure remains.

`src/ads_driver.jl` now owns traversal, partition trees, shape checks, limits,
status transitions, refinement and result assembly. The former ordinary box,
interval box and polygon construction loops were removed. Geometry supplies
splits and intersections; ordinary polynomial fitting and validated model fitting
supply local candidates. Checkpoint and accepted-step flow monitoring retain
their existing early rejection and restart behavior through the same driver.
Continuity still builds overlapping fits separately, using fresh local fits.

The interval extension entry point now includes coherent model, local ADS,
polygon and continuity files. Retained polynomial arithmetic still uses the
original kernels. Live and compiled models share a monomial bounder with integer
powers; outputs of a patch share owned coordinate/exponent snapshots and power
calculations. Public copies and accessors remain independent, including BigFloat
data. Interval boxes retain trees for closed point and intersecting-subbox queries.

`adaptive_map` is the canonical constructor and refiner for all supported static
estimators/geometries. Refinement now saves/defaults to the selected estimator
and requested retained order, rather than the actual degree. Validated box
refinement is supported and reevaluates the original callback on every child.
`validated_adaptive_map` is only a wrapper: its defaults remain `order=3`,
`atol=1e-6`, `splitter=:width`; canonical defaults remain `order=5`, `atol=1e-8`,
`splitter=:tail`. Fresh ADS construction defaults and numerical certification
contracts are unchanged. `continuous_map` also inherits the source's saved
retained order and estimator. The interval split heuristic now excludes satisfied
outputs and exact affine variation, using nonlinear retained tails and coefficient uncertainty
with a relative-width fallback. Acceptance still uses the complete rigorous bound.
The isolated baseline split both requested large-affine-plus-exponential examples
along coordinate one; both now split first along coordinate two.
Fixed subnormal ordinary coordinates now retain their exact center and zero radius.
Ordinary construction and refinement also clone BigFloat endpoint storage rather
than sharing it with caller bounds or a source map; copies/accessors stay independent.

### Paired measurements

`unified_ads.jl` measures warmed minima over 30 batches of 10 calls (100 batches
for queries), plus allocations and snapshot size. The same script was run against
the isolated baseline and this checkout on Windows, Julia 1.13.1 and
IntervalArithmetic 1.0.12. Both multiplication budgets completed. Timings below
use the 32 MiB budget and the paired runs after all validation jobs finished;
they are machine dependent. Allocation/storage figures are often more informative.
Raw results for both budgets are retained in
[ADS before](measurements/ads-before.txt), [ADS after](measurements/ads-after.txt),
[polygons before](measurements/polygon-before.txt) and
[polygons after](measurements/polygon-after.txt).

The two-output map is `[exp(x+y/2), sin(x*y)]` on `[-1/2,1/2]^2`, order three,
`atol=1e-5`; queries use `(1/8,-1/4)`.

| Operation | Before | After | Work / error |
|:--|--:|--:|:--|
| Ordinary construction | 275.31 μs / 906,796 B | 268.10 μs / 907,548 B | 59 patches; estimate `9.212827e-6` |
| Ordinary reusable query | 0.04 μs / 0 B | 0.04 μs / 0 B | Same stored map |
| Validated construction | 3.593 ms / 4,428,679 B | 3.703 ms / 4,352,821 B | 62 patches; bound `8.977223e-6` |
| Validated point query | 41.54 μs / 53,296 B | 9.59 μs / 10,368 B | Widths `[7.076203e-6, 5.016032e-6]` |
| Validated full enclosure | 224.26 μs / 325,248 B | 143.73 μs / 153,328 B | Widths `[1.647929, 0.495384]` |

All patch counts and errors/widths in this table are unchanged. The ordinary
snapshot grows from 71,432 to 72,512 bytes to save estimator/order settings;
the validated snapshot changes from 89,360 to 87,648 bytes, including its tree.
The modest validated construction slowdown is reported alongside the query gains.

The six-variable nonlinear case is
`exp(x1+x2/2)*cos(x3-x4)+log(2+x5*x6)` on `[-1/50,1/50]^6`, order three,
`atol=1e-6`. The eight-output case adds `k*x1`, `k=1:8`, to that expression.

| Operation | Before | After |
|:--|--:|--:|
| Six-variable ordinary construction | 172.96 μs / 296,532 B | 184.01 μs / 295,828 B |
| Six-variable ordinary reusable query | 0.07 μs / 0 B | 0.07 μs / 0 B |
| Six-variable validated construction | 104.25 μs / 171,630 B | 105.75 μs / 173,758 B |
| Six-variable validated point query | 4.32 μs / 6,592 B | 4.22 μs / 6,512 B |
| Six-variable validated full enclosure | 5.71 μs / 7,792 B | 5.73 μs / 7,712 B |
| Eight-output construction | 175.13 μs / 340,528 B | 159.13 μs / 263,206 B |
| Eight-output point query | 24.54 μs / 37,744 B | 20.13 μs / 6,656 B |
| Eight-output full enclosure | 36.01 μs / 47,008 B | 26.94 μs / 7,856 B |
| Eight-output snapshot size | 36,432 B | 16,368 B |

All these cases use one patch before and after. The ordinary six-variable estimate
is unchanged at `5.276273096438478e-7`; its snapshot grows from 118,400 to 118,432
bytes. GuardedTail uses working degree five here, whereas IntervalBound uses three;
these timings compare revisions within each method, not identical work across methods.
The validated scalar full width is
`0.062307791912239896`; every point width is `7.919122362309849e-7`.
The scalar bound is `4.083095687002589e-7`; the maximum eight-output bound is
`4.083095687280145e-7`, with full widths from `0.10230779191223993` through
`0.3823077919122402`. These values are unchanged. Single-output storage grows
from 5,240 to 5,336 bytes; metadata sharing primarily benefits multiple outputs.

Point scaling uses `x^2` on `[-1,1]`, order one, `splitter=:width`, and tolerance
`1/n^2`. Queries are at `1/7`, away from shared faces.

| Patches | Before query / bytes | After query / bytes | Point width / error bound |
|--:|--:|--:|:--|
| 8 | 4.05 μs / 3,792 B | 0.67 μs / 992 B | `0.015625` / `0.015625` |
| 32 | 14.50 μs / 13,008 B | 0.68 μs / 992 B | `0.0009765625` / `0.0009765625` |
| 128 | 56.47 μs / 49,872 B | 0.68 μs / 992 B | `6.103516e-5` / `6.103516e-5` |
| 512 | 222.75 μs / 197,328 B | 0.69 μs / 992 B | `3.814697e-6` / `3.814697e-6` |

The point widths/errors and patch counts are unchanged. Queries on shared faces
deliberately visit both closed branches. Full enclosures must still visit all
intersecting leaves.

`polygon_ads.jl` was also run on both revisions after validation finished.
It reports warmed minima of ten individual calls, with the same 32 MiB budget:

| Case / geometry | Before → after patches | Construction before → after | Full-enclosure before → after | Unchanged full / point width / bound |
|:--|--:|:--|:--|:--|
| `(x+y)^2`, axis boxes, order one | 64 → 64 | 0.750 → 0.925 ms; 933,199 → 1,087,165 B | 97.6 → 33.9 μs; 91,808 → 45,360 B | `4.09375` / `0.15625` / `0.0625` |
| `(x+y)^2`, oriented polygons, order one | 8 → 24 | 0.633 → 1.981 ms; 816,280 → 2,282,439 B | 69.5 → 262.1 μs; 90,784 → 332,216 B | `4.0625` / `0.0625` / `0.0625` |
| `exp(x+y)`, axis boxes, order three | 14 → 14 | 0.300 → 0.338 ms; 389,736 → 423,583 B | 34.3 → 21.6 μs; 43,040 → 33,392 B | `2.405255` / `0.000672592` / `0.000823975` |
| `exp(x+y)`, oriented polygons, order three | 4 → 4 | 0.362 → 0.392 ms; 486,232 → 482,440 B | 33.3 → 41.0 μs; 49,216 → 58,256 B | `2.365317` / `0.000442429` / `0.000442429` |

At order one no nonlinear retained coefficients can identify a direction for a
scalar remainder. The documented width fallback changes the quadratic partition
and makes this example more expensive. No trial-child solves or guessed directional
certificate were added. The axis-box cases remain 64 quadratic and 14 exponential
patches, with their previous widths/errors. Different estimators and geometries
are not expected to produce identical partitions.

The existing `interval_models.jl` and `continuous_ads.jl` benchmarks were also
run sequentially on both revisions with no other Julia jobs running. They report
warmed minima of ten individual calls. Raw results for both budgets are in
[models before](measurements/models-before.txt),
[models after](measurements/models-after.txt),
[continuity before](measurements/continuity-before.txt) and
[continuity after](measurements/continuity-after.txt). The following tables use
the 32 MiB budget.

The standalone model benchmark uses the same six-variable expression and box
as above. Its ADS construction reuses the already initialized interval algebra;
compare revisions within this table rather than its construction times with
the fresh-algebra benchmark above.

| Operation | Before | After |
|:--|--:|--:|
| Live model construction | 56.6 μs / 97,705 B | 65.6 μs / 97,705 B |
| Live model enclosure | 4.1 μs / 5,280 B | 5.2 μs / 5,280 B |
| Validated ADS construction | 73.2 μs / 129,351 B | 84.8 μs / 131,431 B |
| Validated ADS full enclosure | 7.2 μs / 7,792 B | 6.0 μs / 7,712 B |
| Ordinary polynomial construction | 0.8 μs / 5,536 B | 1.0 μs / 5,536 B |
| Stored-polynomial enclosure | 10.1 μs / 6,000 B | 16.4 μs / 6,000 B |

Both revisions use one ADS patch. The model and ADS full width remains
`0.062307791912239896`, the live remainder width `7.919122357670971e-7`, and
the rigorous ADS error bound `4.083095687002589e-7`. The stored-polynomial width
remains `0.062307000000003665`; it does not include generating-function
truncation error. Reusable ordinary multiplication and evaluation allocate zero
bytes on both revisions. Their individual-call timings reach timer resolution;
the batched measurements above provide more useful timing evidence.

Continuity uses fresh order-four C2 overlap fits over order-three source maps,
with source tolerance `1e-5`. The two-dimensional function is
`exp(x+y)+(x-y)^2` on `[-1/4,1/4]^2`; the six-variable function is
`exp(sum(x)/6)*cos(x1*x2)+sum(xi^2)` on `[-1/100,1/100]^6`.
The smooth query measures the surrogate; the enclosure measures the original
function with its overlap certificates.

| Case | Source patches before → after | Source construction before → after | Overlap construction before → after |
|:--|--:|:--|:--|
| Two-dimensional boxes | 29 → 29 | 0.809 → 0.952 ms; 996,728 → 1,062,623 B | 1.293 → 1.266 ms; 1,535,191 → 1,538,015 B |
| Two-dimensional oriented polygons | 20 → 36 | 1.874 → 3.590 ms; 2,368,648 → 4,220,717 B | 0.960 → 1.779 ms; 1,189,671 → 2,160,133 B |
| Six-variable boxes | 1 → 1 | 150.4 → 157.3 μs; 216,647 → 218,823 B | 620.1 → 633.7 μs; 787,376 → 788,560 B |

| Case | Source point before → after | Smooth point before → after | Original full enclosure before → after |
|:--|:--|:--|:--|
| Two-dimensional boxes | 18.6 → 3.1 μs; 20,864 → 4,928 B | 56.4 → 55.0 μs; 396,664 → 408,168 B | 176.4 → 199.6 μs; 202,752 → 207,024 B |
| Two-dimensional oriented polygons | 100.6 → 51.7 μs; 153,184 → 80,744 B | 34.6 → 83.4 μs; 207,760 → 623,600 B | 122.5 → 213.7 μs; 139,552 → 255,536 B |
| Six-variable boxes | 21.0 → 21.7 μs; 8,320 → 8,240 B | 21.7 → 20.6 μs; 37,360 → 37,424 B | 41.9 → 46.6 μs; 18,432 → 18,496 B |

For two-dimensional boxes, full width `1.0620090517687055`, point width
`1.2333677876785742e-6` and uniform blend error `1.1741914258388644e-6` are
unchanged. Oriented polygons have a slightly tighter full width,
`1.1172706627671407` → `1.1065284752671407`; point width
`1.6854236228880382e-7` and blend error `6.43411277441072e-7` are unchanged.
Their changed source partition increases overlap fitting and smooth evaluation
cost. Gradient/Hessian minima also rise from 43.4/44.8 μs to 185.1/198.4 μs;
full timings and allocations for every case are in the raw records.
The six-variable full width `0.02069200590077569`, point width
`1.0219003421241268e-10` and blend error `5.1092549499346795e-11` are unchanged.

### Verification and scope

The final full suites pass **9,602/9,602** on Julia 1.10.12 (4m17.8s) and
Julia 1.13.1 (5m57.4s),
including ordinary Float32/Float64 allocation guarantees, independent rational
coefficient/product oracles, incoming remainders, guarantee/decorations and domain
rejections. The shared ADS contracts exercise all four estimators on boxes and
polygons; the targeted vector/scalar split regressions exercise Float64/BigFloat
and both multiplication budgets. Interval fixtures are explicitly included by
each consuming test file.
After the last scoring cleanup, all **1,622** focused ordinary/interval ADS,
polygon, checkpoint-flow and online-flow checks also passed on Julia 1.13.1.
The subsequent BigFloat ownership follow-up passes all **590** ADS contracts,
including **90** construction, refinement, copy and accessor ownership checks.
Both complete suites also pass **9,602/9,602** with **four Julia threads**, matching
the CI thread setting: Julia 1.10.12 (5m46.9s) and Julia 1.13.1 (7m17.2s).
This exercises the existing concurrent-caller test with actual parallel workers.

Literate checks pass **12/12**. The complete Documenter build executes all 26
examples and passes doctests, references and exported API checks; built navigation
and figure tests pass **237/237**. The build passed on Julia 1.13.1 and, after the
ownership follow-up, Julia 1.10.12 with a freshly resolved temporary docs environment.
The existing ignored docs manifest targets Julia 1.13, so it could not be reused
on 1.10 (PrecompileTools used a newer Julia internal API). No dependency constraint
was changed. Runic and `git diff --check` pass. Commands:

```powershell
julia +1.10 --startup-file=no --project -e 'using Pkg; Pkg.test(; julia_args=["--startup-file=no", "--check-bounds=yes"])'
julia --startup-file=no --project -e 'using Pkg; Pkg.test(; julia_args=["--startup-file=no", "--check-bounds=yes"])'
julia +1.10 --threads=4 --startup-file=no --project -e 'using Pkg; @assert Threads.nthreads()==4; Pkg.test(; julia_args=["--threads=4","--startup-file=no","--check-bounds=yes"])'
julia --threads=4 --startup-file=no --project -e 'using Pkg; @assert Threads.nthreads()==4; Pkg.test(; julia_args=["--threads=4","--startup-file=no","--check-bounds=yes"])'
julia --startup-file=no --project="$env:TEMP/da-unified-ads" --check-bounds=yes -e 'using Test,DifferentialAlgebra,LinearAlgebra; include("test/domain_splitting.jl"); include("test/polygon_ads.jl"); include("test/ads_contracts.jl")'
julia --startup-file=no --project="$env:TEMP/da-unified-ads" --check-bounds=yes -e 'using Test,DifferentialAlgebra; include("test/ads_contracts.jl")'
julia --startup-file=no --project=docs docs/test_literate.jl
julia --startup-file=no --project=docs docs/make.jl
$adsDocsEnvironment = Join-Path $env:TEMP 'da-docs-110'
New-Item -ItemType Directory -Path $adsDocsEnvironment -Force | Out-Null
Copy-Item -LiteralPath docs/Project.toml -Destination (Join-Path $adsDocsEnvironment 'Project.toml')
julia +1.10 --startup-file=no --project="$adsDocsEnvironment" -e 'ENV["JULIA_PKG_PRECOMPILE_AUTO"]="0"; using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate(); include("docs/test_literate.jl"); include("docs/make.jl")'
julia --startup-file=no --project=@runic -e 'using Runic; exit(Runic.main(ARGS))' -- --check --diff --docstrings src test ext examples benchmark docs/make.jl docs/literate.jl docs/test_literate.jl docs/test_navigation.jl
git diff --check
julia --startup-file=no --project="$env:TEMP/da-ads-baseline-env" --check-bounds=yes "$env:TEMP/da-ads-baseline/test/runtests.jl"
julia --startup-file=no --project="$env:TEMP/da-ads-baseline-env" benchmark/unified_ads.jl
julia --startup-file=no --project="$env:TEMP/da-unified-ads" benchmark/unified_ads.jl
julia --startup-file=no --project="$env:TEMP/da-ads-baseline-env" benchmark/polygon_ads.jl
julia --startup-file=no --project="$env:TEMP/da-unified-ads" benchmark/polygon_ads.jl
julia --startup-file=no --project="$env:TEMP/da-ads-baseline-env" benchmark/interval_models.jl
julia --startup-file=no --project="$env:TEMP/da-unified-ads" benchmark/interval_models.jl
julia --startup-file=no --project="$env:TEMP/da-ads-baseline-env" benchmark/continuous_ads.jl
julia --startup-file=no --project="$env:TEMP/da-unified-ads" benchmark/continuous_ads.jl
```

The two temporary benchmark environments contain the baseline worktree and
current checkout respectively, with IntervalArithmetic. For another machine,
use an environment with the intended local revision and declared dependencies.
Package dependencies and the Ubuntu-only workflows are unchanged. Ubuntu CI
itself and other compatible IntervalArithmetic releases were not run locally.
TaylorModels/TaylorSeries comparisons were not added as dependencies or executed.
Certified differentiation, inversion, ODE integration and relative tolerance
remain unsupported. The simple interval tail bounder and polygon covers can be
conservative; the direction heuristic does not attribute a scalar remainder.

The records below describe earlier changes and measurements; their patch counts,
verification totals and commands are historical.

`interval_models.jl` is a separate benchmark, not part of the default test run.
It measures warmed construction/enclosure time, allocations, and interval widths
independently, for a six-variable uncertainty map at total degree three. The
minimum of ten timed calls is reported; compilation is excluded. Timing is
machine-dependent, and enclosure widths are not accuracy comparisons when the
objects have different contracts.

Run with an environment containing the local package and IntervalArithmetic:

```powershell
julia --project=examples -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia --project=examples benchmark/interval_models.jl
julia --project=examples examples/interval_models.jl
julia --project=examples benchmark/polygon_ads.jl
julia --project=examples examples/polygon_ads.jl
```

The example saves `results/interval_models.png` and `results/interval_fit_errors.png`
and is also executed by Literate
in the documentation. Its bands enclose entire displayed physical cells; sampled
function curves are illustrations only. At order three and tolerance `1e-5`, the
two-variable example used 16 patches. Its uniform fit-error bound decreased from
`0.001203788916394627` on the parent to at most `8.217570041330492e-6` on a child.

## Oriented polygons

The interval error criterion is selectable through the existing static ADS API:
`adaptive_map(f, box; estimator=IntervalBound(), ...)`. The other three methods
remain heuristic. `splitter=:oriented` selects 2D convex polygon geometry for
any of the four estimators. Automatic gradient/Hessian sensitivity or supplied
projection rows determine a fixed linear frame; exact rational clipping preserves
the physical partition. Every expansion covers its polygon's parallelogram.
The interval method validates on that complete cover and retains its remainder.

`polygon_ads.jl` compares warmed construction, allocations, enclosure costs,
patch counts and widths separately. The following run used Julia 1.13.1,
IntervalArithmetic 1.0.12 and a 32 MiB multiplication-table budget:

| Case / geometry | Patches | Construction | Construction bytes | Full width | Point width at (1/4,1/4) |
|:--|--:|--:|--:|--:|--:|
| `(x+y)^2`, axis boxes | 64 | 2.074 ms | 933,183 | 4.09375 | 0.15625 |
| `(x+y)^2`, oriented polygons | 8 | 4.281 ms | 816,280 | 4.0625 | 0.0625 |
| `exp(x+y)`, axis boxes | 14 | 0.778 ms | 389,736 | 2.405255 | 0.000672592 |
| `exp(x+y)`, oriented polygons | 4 | 1.248 ms | 486,232 | 2.365317 | 0.000442429 |

The quadratic uses order one, `[-1,1]^2` and `atol=1/16`; the exponential uses
order three, `[-1/2,1/2]^2` and `atol=1/1000`. The uniform fitting-error upper
bounds were respectively `1/16` for both geometries, `0.000823975` for the
exponential boxes and `0.000442429` for its polygons. Full-domain enclosure cost
was 0.202/0.200 ms for the quadratic and 0.083/0.096 ms for the exponential
(boxes/polygons). Both multiplication settings were measured by the script.
Timings are machine/load dependent: these measurements ran alongside validation.
Fewer polygons did **not** imply faster construction in this run.

`examples/polygon_ads.jl` writes `results/polygon_ads.png` (partitions, interval
bands, certified fitting errors) and `results/polygon_triangle.png` (a convex
triangle input domain). Both are embedded in the Literate documentation.

## Measurement on this checkout

Julia 1.13.1, Windows, IntervalArithmetic 1.0.12, Float64 endpoints:

| Operation | No multiplication table | 32 MiB table budget |
|:--|--:|--:|
| Six-variable model construction | 139.9 μs / 97,641 B | 138.9 μs / 97,609 B |
| Model enclosure | 8.4 μs / 5,280 B | 8.6 μs / 5,280 B |
| Certified box ADS construction | 256.8 μs / 129,367 B | 194.6 μs / 129,335 B |
| ADS enclosure | 13.1 μs / 7,792 B | 13.4 μs / 7,792 B |
| Ordinary polynomial construction | 2.7 μs / 8,272 B | 2.3 μs / 5,536 B |
| Reusable ordinary multiplication | 0.6 μs / **0 B** | 0.3 μs / **0 B** |
| Reusable ordinary evaluation | 0.1 μs / **0 B** | 0.1 μs / **0 B** |

Both model/ADS enclosures had width `0.062307791912239896`, with absolute remainder
width `7.919122357670971e-7`. The six-variable ADS case needed one patch at `atol=1e-6`;
its uniform fit-error bound was `4.083095687002589e-7`. The ordinary stored-polynomial
enclosure had width `0.062307000000003665`; it omits the original function's
truncation error and earlier floating coefficient rounding.

## Optional C0/C1/C2 continuous ADS

`continuous_map(f, fit; continuity=:c2, overlap=1//4)` constructs an owned
partition-of-unity surrogate from fresh fits on overlapping covers of the
source patches. It accepts ordinary boxes, interval boxes and oriented polygons.
Linear/cubic/quintic taper endpoint jets give C0/C1/C2 respectively. ForwardDiff
can differentiate the numeric surrogate; an absolute function remainder does
not certify original-function derivatives or a minimizer. Default ADS is
unchanged. Construction restores the caller's algebra, and snapshots survive
reinitialization. The optional ForwardDiff weak extension does not load it when
DifferentialAlgebra is loaded alone.

For an interval source, midpoint coefficients define the explicitly requested
numeric surrogate. Fresh certificates preserve coefficient widths, normalization,
domains and all remainders. Convex blending of the per-support function errors
proves the reported uniform `error_bounds`. `enclose` retains original-function
certificates, rather than discarding metadata or treating the numeric surrogate
as an exact function. Errors on enlarged supports are checked anew if `atol`
is supplied. The original source tolerance is not inherited. Numeric evaluation
roundoff is separate; enclosure endpoints need not be continuous.

Run:

```powershell
julia --startup-file=no --project=docs examples/ads_continuity.jl
julia --startup-file=no --project=docs examples/ads_optimization.jl
julia --startup-file=no --project=docs benchmark/continuous_ads.jl
julia +1.10 --startup-file=no --project=C:/Users/jyar540/AppData/Local/Temp/da-interval-julia110 --check-bounds=yes test/continuous_ads.jl
```

The new examples save four PNGs under `results/`: `ads_continuity.png`,
`ads_continuity_intervals.png`, `ads_optimization.png`, and
`ads_optimization_intervals.png`. They are also executed and embedded by Literate.
The first compares raw/C0/C1/C2 values, gradients and Hessians. Its deliberately
coarse cubic fit shows that continuity is separate from derivative accuracy.
The quadratic interval example has original patch error `[0,1/64]`; enlarged
supports and convex blending give `[0,9/256]`, checked against `atol=1/16`.
The optimization example uses a C2 objective, ForwardDiff gradient/Hessian and
damped Newton with an in-domain line search. Its analytical minimizer is an
independent comparison; the displayed value intervals do not certify an optimum.
Its oriented frame produces more patches than the box frame, demonstrating that
direction selection does not guarantee a performance improvement.

`benchmark/continuous_ads.jl` measures source/overlap construction, numeric
evaluation, ForwardDiff gradient/Hessian, original-function enclosure, bytes
and widths separately, with both multiplication settings. The six-variable
uncertainty case uses `[-1/100,1/100]^6`; an initial `[-1/4,1/4]^6` run exhausted
the default patch budget at `atol=1e-5`. This simple interval bounder can be
conservative in many dimensions; increasing limits is not a performance fix.
Overlapping queries currently scan supports, allocate and use exact geometry
for branch decisions. Existing ordinary polynomial buffer guarantees are
preserved; the new wrapper does not promise allocation-free evaluation.

The completed warmed benchmark on Julia 1.13.1 / IntervalArithmetic 1.0.12,
with a 32 MiB multiplication table budget, reported:

| Case | Patches | Overlap construction / bytes | Numeric query / bytes | Gradient / bytes | Hessian / bytes | Original-function point width |
|:--|--:|--:|--:|--:|--:|--:|
| 2D boxes | 29 | 4.685 ms / 1,535,175 | 139.8 μs / 396,664 | 221.8 μs / 398,904 | 243.3 μs / 403,551 | `1.23337e-6` |
| 2D oriented polygons | 20 | 3.063 ms / 1,189,671 | 73.8 μs / 207,760 | 82.1 μs / 209,904 | 92.8 μs / 214,592 | `1.68542e-7` |
| Six-variable box | 1 | 1.528 ms / 787,392 | 44.9 μs / 37,360 | 46.9 μs / 39,840 | 67.1 μs / 55,884 | `1.02190e-10` |

The six-variable full-domain enclosure width was `0.02069200590077569`, with
uniform blend-error upper bound `5.1092549499346795e-11`. Source construction
took 0.362 ms / 216,647 bytes; full original-function enclosure took
88.8 μs / 18,432 bytes. Both multiplication settings completed. Timings ran
alongside validation and are machine/load dependent. The 2D box and polygon
full-domain widths were respectively `1.062009` and `1.117271`; the tighter
point width of polygons does not imply a tighter full-domain bound.

The current continuity changes started from `19f9900`, which already contained
the interval foundation, common ADS estimator, polygon geometry and test-fixture
migration. The focused continuity file passed **444/444 on Julia 1.10.12**:
66 exact taper/one-sided jet oracles, 226 geometry/ownership/AD checks and
152 certified function-bound checks. Endpoint regularity checks use exact
rational polynomial/derivative identities; the quadratic inclusion checks use
the analytical Taylor remainder. Samples and plot smoothness are supplementary.

The final full `Pkg.test` command listed below passed **8,853/8,853** on
Julia 1.13.1 (11m52.7s), including all 444 continuity checks and the existing
Float32/Float64 reusable-buffer allocation tests. The full suite on Julia 1.10
was not repeated; its focused continuity file passed. Runic and `git diff --check`
passed. `docs/test_literate.jl` passed 6/6. An initial documentation run exposed
Literate treating indented single-hash comments inside the new loops as Markdown;
they now use escaped code comments. The corrected two pages executed with both
figures embedded, and all four saved PNGs were visually checked.
All **26 registered examples** then executed successfully and the complete
Documenter build passed, including doctests, cross-references and exported API
coverage. Already successful examples were retained while the corrected pages
and remaining interval/storage examples were rerun. The exact recovery commands
were:

```powershell
julia --startup-file=no --project=docs -e 'using Literate; include("docs/literate.jl"); for name in ("ads_continuity", "ads_optimization"); Literate.markdown("examples/" * name * ".jl", "docs/src/generated"; flavor=Literate.DocumenterFlavor(), execute=true, postprocess=m->example_images(m,"docs/src/generated",name)); end'
julia --startup-file=no --project=docs C:/Users/jyar540/AppData/Local/Temp/da-continuity-finish-docs.jl
```

The temporary recovery script reran `interval_models` and `serialization`,
asserted all pages existed and both new pages had two PNGs each, then executed
the unchanged `DocMeta`/`makedocs` body of `docs/make.jl`. A fresh build uses the
ordinary `docs/make.jl` command below; no recovery mode was added to the package.

## Verification commands

The initial checkout was `c01d589019a53ef7d2c8b87ac0ba448a49a5a56a` on `cleanup`,
newer than the reviewed `1c99682`. The working-order regression is separate:
scalar addition, unary plus and power one now truncate consistently with other
arithmetic; copying/conversion still preserve stored coefficients.

The commands used on this machine were:

```powershell
julia --startup-file=no --project=. -e 'using Pkg; Pkg.test(; julia_args=["--startup-file=no", "--check-bounds=yes"], allow_reresolve=true)'
julia --startup-file=no --project=C:/Users/jyar540/AppData/Local/Temp/da-interval-tests --check-bounds=yes test/runtests.jl
julia --startup-file=no --project=C:/Users/jyar540/AppData/Local/Temp/da-interval-tests --check-bounds=yes -e 'include("test/interval_models.jl"); include("test/validated_ads.jl")'
julia --startup-file=no --project=C:/Users/jyar540/AppData/Local/Temp/da-interval-tests --check-bounds=yes -e 'using DifferentialAlgebra,IntervalArithmetic; include("test/utils.jl"); include("test/coefficient_types.jl"); include("test/regressions.jl"); include("test/interval_models.jl"); include("test/validated_ads.jl")'
julia +1.10 --startup-file=no --project=C:/Users/jyar540/AppData/Local/Temp/da-interval-julia110 --check-bounds=yes -e 'include("test/interval_models.jl"); include("test/validated_ads.jl")'
julia +1.10 --startup-file=no --project=C:/Users/jyar540/AppData/Local/Temp/da-interval-julia110 --check-bounds=yes -e 'include("test/statistics.jl"); include("test/regressions.jl")'
julia +1.10 --startup-file=no --project=C:/Users/jyar540/AppData/Local/Temp/da-interval-julia110 --check-bounds=yes -e 'include("test/polygon_ads.jl"); include("test/interval_models.jl"); include("test/validated_ads.jl"); include("test/interval_polygon_ads.jl")'
julia --startup-file=no --project=C:/Users/jyar540/AppData/Local/Temp/da-interval-tests benchmark/interval_models.jl
julia --startup-file=no --project=C:/Users/jyar540/AppData/Local/Temp/da-interval-tests benchmark/polygon_ads.jl
julia --startup-file=no --project=docs examples/interval_models.jl
julia --startup-file=no --project=docs examples/polygon_ads.jl
julia --startup-file=no --project=docs docs/test_literate.jl
julia --startup-file=no --project=docs docs/make.jl
julia --startup-file=no --project=@runic -e 'using Runic; exit(Runic.main(ARGS))' -- --check --diff --docstrings src test ext examples benchmark docs/make.jl docs/literate.jl docs/test_literate.jl
git diff --check
```

The complete `Pkg.test` run passed **8,409/8,409** on Julia 1.13.1 (8m34.0s),
including the existing ordinary allocation tests, migrated Julia coefficient/
moment oracles and all interval/box/polygon ADS checks. It resolved the test environment from
the reduced package extras. The two focused interval files passed
**1,054/1,054** on Julia 1.13.1 and
1.10.12: interval coefficients/ranges/models (762) and certified ADS/snapshots
(292). Existing coefficient-type allocation tests also passed with the interval
extension loaded first. Their decisive
inclusion checks use independent rational polynomial coefficients, exact uniform
monomial tail bounds and analytical Taylor-theorem remainders. They also check
guarantee/decorations, both table settings, Float64/BigFloat, extension loading,
copying/lifetimes, working order, coefficient tolerance and failure cases.

Documentation/Literate checks passed (6/6), every registered example executed,
and the final complete documentation build succeeded with all 24 examples and
both polygon figures embedded. The fitting scripts also save their figures as
PNGs for review without building the documentation.
Formatting and `git diff --check` passed. Expected NaI warnings come from tests
of invalid interval data; method-overwrite warnings come from exercising and
restoring IntervalArithmetic's public rounding configuration.

The full suite on Julia 1.10 and differential comparisons with TaylorModels.jl
have not been run. TaylorModels/TaylorSeries are not runtime or test dependencies.
The interval implementation was exercised with IntervalArithmetic 1.0.12; other
allowed 1.x releases were not individually tested. Polygon ADS supports 2D
positive-area convex domains and a fixed frame. Higher-dimensional polyhedra,
nonconvex domains, per-child changing frames, polygon-constrained range
optimization, verified time integration/inversion/differentiation and
relative-tolerance certified ADS remain outside this implementation.

The subsequent test-fixture migration removed `test/reference`, `test/data`,
the external C-core generator, and the DelimitedFiles/TOML test dependencies.
Scalar ForwardDiff derivatives plus an independent multinomial expansion now
check all 84 coefficients of the twelve function cases with both multiplication
backends. Exact rational Gaussian contractions replace the moment tables, with
the small means/covariances kept inline. The migrated statistics/regression files
passed **3,991/3,991** on Julia 1.10.12 (1,584 moments and 2,407 regressions).
The contribution guide documents the Julia-oracle and CSV-fixture conventions.

After adding the common interval estimator and polygon geometry, the focused
Julia 1.10.12 run passed **2,281/2,281**: ordinary polygon ADS (897), interval
coefficients/ranges/models (762), certified box ADS (292) and certified polygon
ADS (330). The polygon checks use exact area/intersection oracles and the
analytical quadratic remainder, including nonorthogonal frames, shared faces,
subpolygons, degenerate point/line queries, BigFloat precision changes, ownership
and algebra restoration. The numerical plots are supplementary illustrations.
