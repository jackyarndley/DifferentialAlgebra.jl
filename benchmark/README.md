# Interval models and ADS

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
two-variable example used 16 leaves. Its uniform fit-error bound decreased from
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
leaf counts and widths separately. The following run used Julia 1.13.1,
IntervalArithmetic 1.0.12 and a 32 MiB multiplication-table budget:

| Case / geometry | Leaves | Construction | Construction bytes | Full width | Point width at (1/4,1/4) |
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
width `7.919122357670971e-7`. The six-variable ADS case needed one leaf at `atol=1e-6`;
its uniform fit-error bound was `4.083095687002589e-7`. The ordinary stored-polynomial
enclosure had width `0.062307000000003665`; it omits the original function's
truncation error and earlier floating coefficient rounding.

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
