# Automatic domain splitting

Choose the error method with `estimator` and the geometry with `splitter`.
[`IntervalBound`](@ref) is the optional validated method for static maps; the
other three methods use coefficient indicators and optional point checks.
Neither those indicators nor their point checks implicitly perform interval
arithmetic. `splitter=:oriented` adds a 2D convex polygon partition to all four
static-map methods. Existing flow methods keep their heuristic box contract.

| Example calculation | Estimator comparison | Geometry and evidence |
|:--|:--|:--|
| [Gaussian and diagonal functions](../generated/ads_methods.md) | All four methods | Same-domain partitions; sampled discrepancies and separate interval bounds |
| [Explicit orbital phase states](../generated/ads_kepler.md#All-four-estimators-on-an-explicit-orbital-phase-map) | All four methods | Boxes; complete-cell interval bands, separate from the unverified implicit Kepler solve |
| [Initial-epoch energy and apses](../generated/ads_flow.md#Interval-ADS-on-initial-epoch-orbit-diagnostics) | All four methods | Boxes; static diagnostic certificates, separate from numerical time propagation |
| [Inclined orbital ribbon](../generated/orbit_polygon_ads.md) | All four methods × two geometries | Boxes and diagonal polygons; correlated hexagonal prior and original-state enclosures |
| [Mars B-plane flyby](../generated/flyby_uncertainty.md) | All four box methods plus interval polygons | Whole-cell altitude classification; guaranteed clear, intersecting, or unresolved cells |

A single Taylor expansion may lose accuracy over a large input box.
Automatic domain splitting (ADS) replaces it with a collection of local
expansions [Wittig2015](@cite). Each patch uses normalized coordinates on
`[-1, 1]`, while the public map accepts physical coordinates.

The original box methods support static functions, refinement and propagated
states. The interval and polygon extensions currently support static maps:

| Operation | Interface |
|:--|:--|
| Approximate a function | `adaptive_map(f, lower, upper; ...)` |
| Refine an existing partition | `adaptive_map(f, previous; ...)` |
| Check a flow at specified times | `adaptive_flow(advance, initial, lower, upper, times; ...)` |
| Check a flow after accepted solver steps | `adaptive_flow(advance, initial, lower, upper, (t0, tf); ...)` |

## ADS terminology

Use **subdomain** for an input region produced by splitting and **patch** for
the local Taylor map or Taylor model together with its subdomain and coordinate
normalization. Splitting a subdomain replaces its patch with two child patches.
The input subdomains form a **partition**: they cover the original domain and
have disjoint relative interiors, with shared boundaries. `map.patches` is the patch
collection; `max_patches` limits its size.

In astrodynamics and DA set propagation, a collection of local maps representing
a parameterized set is also called a **DA manifold representation**. The patches
describe portions of the mapped set, such as an orbital uncertainty ribbon.
[Wittig's introduction to DA](https://indico.esa.int/event/98/attachments/2928/3393/IntroductionToDA.pdf)
describes this connection between local DA maps and manifold charts (slides
21–22). This terminology describes the representation; a patch collection alone
does not establish that its image is a smooth manifold. Independently computed
ADS maps can disagree on shared boundaries, and arbitrary callbacks need not
parameterize a regular manifold.

For [`continuous_map`](@ref), enlarged supports overlap and form a **cover**
of the physical domain. The original ADS subdomains still form a partition.
The optional C0/C1/C2 blend controls continuity of the surrogate; it does not
establish manifold geometry or original-function derivative accuracy. Existing
`PiecewiseTaylorMap`, `PiecewiseTaylorModel` and `PiecewisePolygonMap` names
refer to the corresponding patch collections.

## Static maps

```@example ads
using DifferentialAlgebra

f(x) = [exp(x[1] + x[2]), sin(x[1] - x[2])]
map = adaptive_map(f, [-0.5, -0.5], [0.5, 0.5]; order = 4, atol = 1e-6)
map([0.2, -0.1])
```

The callback accepts a vector and returns a real scalar or a nonempty vector
of reals. It must work for numeric and polynomial coordinates. Equal lower
and upper bounds fix a coordinate. The promoted floating-point type of the
bounds determines the input coefficient type.

## Error estimators

Select an estimator with `estimator = ...`:

| Estimator | Expansion degree | Indicator |
|:--|:--|:--|
| `GuardedTail()` (default) | `order + guard_order` | Sum of absolute discarded coefficients plus an extrapolated next degree |
| `ExtrapolatedTail()` | `order` | Exponential fit to maximum coefficient magnitude by degree |
| `LastTerms()` | `order` | Maximum absolute coefficient in the last two retained degrees |
| `IntervalBound()` | `order` | Certified absolute remainder plus retained coefficient uncertainty |

`GuardedTail` defaults to two extra degrees. Their L1 norm bounds their
contribution on the normalized box, but not the uncomputed remainder.
No extrapolation is added when the discarded tail is zero.

`ExtrapolatedTail` fits the logarithms of nonzero degree norms and predicts
the next degree. If fewer than two nonzero degrees are available, it falls
back to the maximum in the last two degrees.
`LastTerms(degrees = 3, safety = 2)` changes the number of retained degrees
and multiplies the indicator by a safety factor. These two methods require
`guard_order = 0`, which is selected automatically.

The first three methods are **heuristic**, not rigorous remainder bounds. Retained-tail
methods may split an exactly represented polynomial unnecessarily. Any
coefficient-only method can miss terms absent from the computed expansion.

For those three methods, independent numeric evaluations at axis endpoints and selected
corners supplement the coefficient indicator. Repeated points are evaluated
only once. Setting `check_points = false` avoids these evaluations, but also
removes their ability to detect missing high-degree terms. Always validate
the map at independent points.

An output passes when

``\mathrm{indicator}_j \le \mathrm{atol}_j + \mathrm{rtol}\,s_j,``

where `s_j` is the largest absolute output at the center and checked points.
`atol` can be a scalar or a vector with one entry per output. Use separate
absolute tolerances for outputs with different units or scales.

`IntervalBound()` instead reevaluates the original expression with native
Taylor models and accepts a patch only when a rigorous uniform absolute error
bound meets `atol`. It includes coefficient widths as well as the remainder.
It requires IntervalArithmetic, `rtol=0`, `guard_order=0` and
`check_points=false` (selected automatically). Its callbacks accept Taylor-model
inputs, and returned values retain their remainder and domain:

```julia
using DifferentialAlgebra, IntervalArithmetic
box = fill(interval(Float64, -1, 1), 2)
fit = adaptive_map(v -> (v[1]+v[2])^2, box;
                   estimator=IntervalBound(), order=1, atol=1//16)
enclose(fit)                     # Original function on the complete domain
fit([1//2, 1//4])               # Interval at a physical point
```

Finite ordinary endpoint vectors can also be passed to `adaptive_map` with
this estimator. `validated_adaptive_map` remains a convenience entry point;
its default splitter is `:width`, whereas `adaptive_map` defaults to `:tail`.
Both apply the same acceptance test. See the
[interval manual](interval-models.md#Certified-domain-splitting) for its proof.

## Split directions

`splitter = :tail` scores each coordinate by its predicted reduction in
tail coefficients when halved, scaled by each output's tolerance. Axis
point discrepancies supplement the scores; ties use relative box widths.

`splitter = :width` bisects the longest side relative to the original
domain. This simple geometric alternative is useful as a baseline, though
it may spend patches on weakly varying coordinates. Fixed or unrepresentably
small coordinates are never selected.

Both children reevaluate the original callback. Translating an already
truncated polynomial alone would not recover omitted terms.

### Oriented and polygonal splitting

```julia
fit = adaptive_map(v -> (v[1]+v[2])^2, [-1, -1], [1, 1];
                   estimator=IntervalBound(), splitter=:oriented,
                   order=1, atol=1//16)
split_directions(fit)            # [1 1; -1 1] for this diagonal dependence
polygon_vertices(domain(first(fit.patches)))

triangle = ConvexPolygon([(0,0), (1,0), (0,1)])
fit = adaptive_map(v -> exp(v[1]+v[2]), triangle;
                   estimator=IntervalBound(), order=3, atol=1//1000,
                   directions=[1 1; -1 1])
```

Automatic directions come from gradient/Hessian sensitivity of an interior
degree-two probe. They are a heuristic frame choice, not a search for a global
optimum. A nonsingular 2x2 `directions` matrix supplies projection rows; it can
be oblique as well as rotated. `directions=:axes` is an unrotated polygon
baseline. The selected frame remains fixed during construction; each split
clips a physical polygon by a line at the midpoint of a projected extent.
`splitter=:width` on a polygon uses relative projected widths; `:oriented` or
`:tail` uses coefficient sensitivity to choose between the two frame directions.

Clipping, inverse coordinates, containment and area use exact rational
geometry. All function expansions cover each polygon's bounding parallelogram.
Certified arithmetic rounds outward on that cover and requires the function
to be valid on the entire cover, which can extend beyond the physical polygon.
No polygon-specific optimization bounder is implemented. Point and subdomain
queries enforce the original physical polygon. Interval queries hull every
overlapping patch, including shared edges. Polygon snapshots survive algebra
changes, and public geometry accessors own their data.

These methods return `PiecewisePolygonMap`. They support 2D positive-area convex
domains, scalar/vector outputs and refinement from an existing polygon map.
Higher-dimensional and fixed-coordinate domains continue to use box ADS.
The [polygon comparison](../generated/polygon_ads.md) plots the partitions and
interval bands and reports patch counts, time, allocations and widths separately.
The [orbital ribbon](../generated/orbit_polygon_ads.md) and
[flyby](../generated/flyby_uncertainty.md) give astrodynamics applications.
Their explicit callbacks avoid an unverified root solve or time integrator.
The Mars example bounds the stated point-mass scattering formula; it is not a
full atmospheric or mission-specific safety analysis.

## Optional continuity for optimization

`continuous_map(f, fit; continuity=:c2)` (also `:c0` or `:c1`) reevaluates the original
function on overlapping box or projected polygon covers and blends the local
fits. C0 gives continuous values, C1 adds continuous gradients, and C2 adds
continuous Hessians of the surrogate. This is separate from ADS error-method
selection and preserves the original map. With IntervalBound sources it retains
function certificates and recomputes uniform errors on the overlaps. The source
tolerance is not inherited; optional `atol` checks the new bounds.
See [continuous ADS](continuous-ads.md) for the partition-of-unity proof,
ForwardDiff usage, accuracy contract and optimization examples. Absolute function
remainders do not certify original-function derivatives or an optimum.

## Refinement

An existing map can supply the starting partition for a tighter tolerance,
different order, or another estimator:

```@example ads
fine = adaptive_map(f, map; atol = 1e-8, estimator = GuardedTail())
(coarse_patches = length(map.patches), fine_patches = length(fine.patches))
```

Existing boundaries are preserved and the input map is unchanged. Supply the
original function as `f`. Using the old map itself would only approximate its
existing polynomials. Depth and patch limits apply to the complete tree,
including the starting partition.

## Checkpointed propagation

Define `initial(parameters)` and `advance(state, (tstart, tend))`. For an
ODE, `advance` can call an adaptive solver and return its final state.
This example uses the exact flow of `du/dt = u²` to isolate spatial
truncation from numerical integration error:

```@example ads
initial(x) = x[1]
advance(u, span) = u / (1 - (span[2] - span[1]) * u)
times = range(0, 1; length = 5)
flow = adaptive_flow(advance, initial, [0.0], [0.6], times; order = 4, atol = 1e-6)
(approximation = flow([0.4]), exact = 0.4 / 0.6, patches = length(flow.patches))
```

The initial state and each supplied checkpoint must pass. On an unsplit
patch, the next interval starts from the current state, retaining guard
coefficients. On failure, the children restart from their initial conditions
at the first time. This recomputation recovers terms discarded in the parent
trajectory. It also means splitting late in propagation can be expensive.

A decreasing time vector requests backward propagation. Checks do not occur
between the supplied times. Choose enough checkpoints to resolve changes in
the flow; choose the ODE solver's tolerances separately.

## Accepted-step propagation

Passing a tuple `(t0, tf)` selects a three-argument propagator:

```julia
advance(state, (tstart, tend), monitor)
```

After each accepted step, call `monitor(state, time)`. A return value of
`true` requests an immediate stop; `false` permits continuation. Return
`(; state, time)` at the stopping time or final endpoint. With OrdinaryDiffEq,
a `DiscreteCallback` can invoke the monitor and call `terminate!`.

The [orbit splitting example](../generated/ads_flow.md) implements both
propagation interfaces with Vern9. Point checks in online mode propagate
independent numeric trajectories to each monitored time. They can be costly;
coefficient-only monitoring with a separate final validation grid is an
alternative.

All completed patches store final-time states. With `strict = false`, even
a patch that fails a check is propagated to the endpoint and marked
unresolved. A flow patch's `error_estimate` is the maximum indicator over
monitored states; each check uses its own relative scale. These indicators
do not bound error between checks or the ODE solver's integration error.

## Inspecting and evaluating results

```@example ads
(converged = map.converged, patch_count = length(map.patches),
 largest_estimate = maximum(maximum(p.error_estimate) for p in map.patches))
```

Each [`TaylorPatch`](@ref) stores bounds, center, radius, depth, status,
a compiled polynomial and per-output indicators. Treat these fields as
read-only. The compiled polynomial accepts normalized coordinates; the
piecewise map normalizes physical coordinates and selects the patch through
a binary tree. Shared faces belong to the lower child. Separate expansions
need not agree exactly on a face.

Reuse buffers for repeated evaluation:

```@example ads
result = zeros(noutputs(map))
normalized = zeros(nvariables(map))
work = zeros(degree(map) + 1)
evaluate!(result, map, [0.2, -0.1], normalized, work)
```

`max_depth` limits splits along a path; `max_patches` limits all patches.
The default `strict = true` raises an error at an unresolved limit.
With `strict = false`, the map covers the entire box and reports
`:max_depth`, `:max_patches` or `:roundoff` on unresolved patches.
Check `map.converged` before using the result.

Construction restores the caller's algebra and existing polynomials even
on failure. Callbacks must not change algebra settings. Construction must
not overlap other polynomial calculations using the global algebra;
numeric evaluation of completed maps is independent of it.

See the [method comparison](../generated/ads_methods.md) for partitions and
measured errors, and the [analytic Kepler map](../generated/ads_kepler.md)
for a nonlinear uncertainty transformation.
