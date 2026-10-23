# Automatic domain splitting

A single Taylor expansion may lose accuracy across a large input region. ADS
replaces it with local expansions on smaller regions [Wittig2015](@cite).
Choose the error method with `estimator` and the geometry with `splitter`;
construction and refinement use the same driver.

## Fit, inspect and refine

```@example ads
using DifferentialAlgebra
f(x) = [exp(x[1]+x[2]), sin(x[1]-x[2])]
fit = adaptive_map(f, [-0.5, -0.5], [0.5, 0.5]; order=4, atol=1e-6)
(value = fit([0.2, -0.1]), converged = fit.converged, patches = length(fit.patches))
```

The callback accepts a vector and returns a real scalar or nonempty real vector.
For ordinary ADS it must work with both numeric and polynomial coordinates.
Bounds select the promoted floating input coefficient type; equal bounds fix a
coordinate. Public queries use physical coordinates; each local polynomial uses
normalized coordinates on `[-1,1]`.

```@example ads
fine = adaptive_map(f, fit; atol=1e-8)
(coarse = length(fit.patches), refined = length(fine.patches))
```

Refinement preserves existing boundaries and leaves the source unchanged. It
inherits the estimator, requested order, tolerances and other construction options;
explicit keywords override compatible options. An ordinary box fit may upgrade to
`IntervalBound`, while a validated box fit must retain it. Polygon refinement may
switch estimator families; switching away from `IntervalBound` drops certified
bounds. Changing the estimator selects compatible guard and point-check defaults.
Always supply the **original callback**: reevaluating it on children recovers terms
that the parent expansion discarded.

## Select an error estimator

| Estimator | Degree computed | Acceptance evidence |
|:--|:--|:--|
| `GuardedTail()` (default) | `order + guard_order` | Discarded coefficient L1 norm plus one extrapolated degree |
| `ExtrapolatedTail()` | `order` | Exponential fit to coefficient magnitude by degree |
| `LastTerms()` | `order` | Largest coefficient in the last two nonlinear degrees |
| `IntervalBound()` | `order` | Uniform absolute remainder and retained coefficient uncertainty |

The first three are heuristic. Independent numeric checks at selected corners
and axis endpoints supplement them by default. `check_points=false` removes
those checks. Coefficient indicators can miss uncomputed terms; validate a
heuristic fit independently. Retained-tail estimators can also split an exactly
represented polynomial unnecessarily.

`GuardedTail` computes two extra degrees by default. A zero discarded tail adds
no extrapolation. `ExtrapolatedTail` fits nonzero nonlinear degree norms, falling
back to the last two degrees when too little information exists.
`LastTerms(degrees=3, safety=2)` changes its window and safety factor. Both retained-
tail methods use `guard_order=0`; their order-one indicators vanish, leaving the
point checks to assess error.

An ordinary output passes when

```math
\mathrm{indicator}_j \le \mathrm{atol}_j + \mathrm{rtol}\,s_j,
```

where `s_j` is its largest magnitude at the center and checked points. Supply
per-output `atol` values when units or scales differ. Defaults are `order=5`,
`atol=1e-8`, `rtol=0` and `splitter=:tail`.

For a rigorous static-function fit, load IntervalArithmetic and use model-compatible
arithmetic in the callback:

```@example certified_ads
using DifferentialAlgebra, IntervalArithmetic
certified = adaptive_map(v -> (v[1]+v[2])^2, [-1, -1], [1, 1];
    estimator=IntervalBound(), order=1, atol=1//16)
(point = certified([1//2, 1//4]), whole_domain = enclose(certified))
```

`IntervalBound` requires `rtol=0`, `guard_order=0` and `check_points=false`, selected
automatically. It accepts interval boxes or finite endpoint vectors and computes
new Taylor models on each child. See [Certified domain splitting](interval-models.md#Certified-domain-splitting)
for the acceptance bound, supported arithmetic and domain-failure handling.

## Split directions

- `splitter=:tail` chooses a coordinate using coefficient sensitivity and point
  discrepancies from unsatisfied outputs. Ties use relative box widths.
- `splitter=:width` bisects the longest side relative to the original box.
- `splitter=:oriented` builds a two-dimensional convex polygon partition.

Fixed or unrepresentably small coordinates cannot split. With `IntervalBound`,
direction scores use the last nonlinear degrees and coefficient uncertainty,
ignoring exact affine variation. The scalar remainder cannot be attributed to a
coordinate, so missing directional information falls back to relative widths.
These heuristics affect efficiency; the interval bound alone decides acceptance.

### Oriented and polygonal splitting

```@example certified_polygon
using DifferentialAlgebra, IntervalArithmetic
triangle = ConvexPolygon([(0,0), (1,0), (0,1)])
fit = adaptive_map(v -> exp(v[1]+v[2]), triangle;
    estimator=IntervalBound(), order=3, atol=1//1000, directions=[1 1; -1 1])
(frame = split_directions(fit), vertices = polygon_vertices(domain(first(fit.patches))))
```

Pass a positive-area 2D convex polygon, or select `:oriented` on a box. A
nonsingular 2×2 `directions` matrix specifies projection rows and may be oblique.
The default sensitivity probe chooses a fixed frame from nonlinear coefficients;
`directions=:axes` supplies an unrotated baseline. No optimal frame is promised.

Each split clips the physical polygon at a projected midpoint. Geometry uses
exact rationals; expansions cover the bounding parallelogram in that frame.
The callback must be valid on the **entire cover**, including points outside
the polygon. Interval arithmetic bounds that cover outward. Nonconvex polygons
and higher-dimensional polyhedra are unsupported; box ADS handles higher
dimensions and fixed coordinates.

The result is a `PiecewisePolygonMap`. Refinement can retain old frames while
introducing new ones. Public point/subdomain queries enforce the physical
polygon, and certified queries include every intersecting closed patch.
Ordinary evaluations retain promoted numeric types, including BigFloat and
ForwardDiff partials; routing uses primal coordinates. Use explicit BigFloat
inputs when precision matters near a boundary.

Compare [partitions and fitting bands](../generated/polygon_ads.md), an
[inclined orbital ribbon](../generated/orbit_polygon_ads.md), and
[B-plane flyby classification](../generated/flyby_uncertainty.md).

## Status, storage and evaluation

`max_depth` limits splits along a path and `max_patches` limits the final
partition. By default, `strict=true` throws when a tolerance remains unmet.
With `strict=false`, valid unresolved patches remain and report `:max_depth`,
`:max_patches` or `:roundoff`; inspect `fit.converged` and patch `status`.
An invalid interval function domain still throws at a splitting limit.

`max_order(fit)` reports the requested degree; `degree(fit)` reports the largest
stored degree. Ordinary patches expose heuristic `error_estimate`; certified
patches expose `error_bounds`. `domain`, `nvariables` and `noutputs` apply to both.
`domain` returns copied endpoints, interval boxes or polygons as appropriate.

`fit.patches` copies its vector: cache it before repeated indexing. Entries are
read-only; `copy(fit)` gives an independent map with copied patch storage. Numeric maps
survive algebra reinitialization. Ordinary queries choose the lower child at a
shared face; certified queries hull all touching closed patches. Independent
expansions need not agree on a face; [continuous ADS](continuous-ads.md) provides
optional C0/C1/C2 blending with newly assessed overlap errors.

Reuse buffers for ordinary box queries:

```@example ads
result = zeros(noutputs(fit))
normalized = zeros(nvariables(fit))
work = zeros(degree(fit) + 1)
evaluate!(result, fit, [0.2, -0.1], normalized, work)
```

Builders use private scoped algebras and restore the caller's binding on failure.
Independent builders may overlap; callbacks must not change their algebra settings.
Global reinitialization must not overlap readers of its existing polynomials.

## Propagation with checkpoints

Define `initial(parameters)` and `advance(state, (tstart, tend))`. Here the exact
flow of `u′=u²` isolates uncertainty truncation from time integration error:

```@example ads
initial(x) = x[1]
advance(u, span) = u / (1-(span[2]-span[1])*u)
flow = adaptive_flow(advance, initial, [0.0], [0.6], range(0, 1; length=5);
    order=4, atol=1e-6)
(approximation = flow([0.4]), exact = 0.4/0.6, patches = length(flow.patches))
```

The initial state and each checkpoint must pass. Unsplit patches continue from
their current state with guard coefficients. Failed patches split and restart
from the initial time to recover discarded terms. Late splits can therefore
be expensive. Decreasing times request backward propagation; no checks occur
between supplied times.

Passing `(t0, tf)` instead selects `advance(state, span, monitor)`. Call
`monitor(state, time)` after accepted steps and stop immediately if it returns
`true`; return `(; state, time)` at the stopping time or endpoint. The
[orbit-flow tutorial](../generated/ads_flow.md) implements both interfaces with
Vern9. Online point checks propagate independent numeric states and can be costly.

Every completed patch stores a final-time state, including unresolved patches
with `strict=false`. Its estimate is the maximum over monitored states. Flow ADS
uses heuristic box methods and does not bound errors between checks or the ODE
solver's integration error; choose solver tolerances separately.

## ADS terminology

A **subdomain** is an input region; a **patch** combines that region with its
local map/model and normalization. Subdomains form a partition with disjoint
interiors and shared boundaries. In DA set propagation their images are also
called a DA manifold representation; this name alone establishes neither a
regular manifold nor smooth joins. Continuous ADS enlarges supports into an
overlapping **cover**, while retaining the original partition as its cores.
