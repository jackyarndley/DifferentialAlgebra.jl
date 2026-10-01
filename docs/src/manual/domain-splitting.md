# Automatic domain splitting

[`adaptive_map`](@ref) approximates a function over a box using a collection of
Taylor maps. It splits domains whose estimated truncation error is too large,
recomputing the original function on each child. This is a map-based form of
automatic domain splitting (ADS), following the coefficient-norm approach of
[Wittig2015](@citet). It can wrap an analytic map or an entire numerical propagation;
it does not insert splitting events into an ODE solver.

```@example ads
using DifferentialAlgebra
f(x) = [inv(1 - x[1]), exp(x[2])]
map = adaptive_map(f, [-0.5, -0.5], [0.5, 0.5]; order = 4, atol = 1e-6)
map([0.2, -0.1])
```

The callback receives a vector and returns a scalar or a nonempty vector. It must
work with ordinary numbers and Taylor polynomials. Bounds are physical coordinates;
equal lower and upper bounds fix a parameter. Call the resulting map with the same
coordinates. Points outside the box are rejected.

## Error estimates and splitting

For each box, the callback is expanded to `order + guard_order` (default: two guard
degrees). The returned patch retains only `order` degrees. Its error estimate adds
the L1 norm of discarded coefficients to a log-linear extrapolation of the
homogeneous coefficient norms to the next degree. If all discarded coefficients
are zero, the coefficient estimate is zero. Numeric checks at axis endpoints and
four selected corners supplement this estimate; they can detect some terms beyond
the computed expansion. This uses a number of checks proportional to dimension.

The split direction favors the largest predicted reduction of discarded terms
after halving a coordinate. Endpoint discrepancies also contribute. Ties favor
the widest remaining fraction of the original domain. Each child reruns the
callback; translating the truncated parent polynomial would not recover missing
terms.

Every output is tested against `atol + rtol * scale`. The scale is the largest
absolute value of that output at the center and checked points. `atol` accepts a
vector when outputs need different tolerances or have different units.

These estimates are **not certified uniform bounds**. A finite expansion and a
finite set of checks can both miss behavior between checked points. Validate the
result independently, particularly near singularities or branches in a callback.
`check_points = false` skips numeric checks and requires particular care with
functions whose first nonzero terms exceed the working order.

## Inspecting and evaluating patches

```@example ads
(patches = length(map.patches), converged = map.converged,
 largest_estimate = maximum(maximum(p.error_estimate) for p in map.patches))
```

Each [`TaylorPatch`](@ref) records its bounds, center, radius, depth, status,
compiled polynomial and per-output error estimate. Treat these fields as read-only.
The compiled polynomial uses normalized coordinates; the piecewise map handles
normalization and finds the patch through a binary tree. A shared face belongs to
the lower child. Separately expanded patches need not be exactly continuous there.

For repeated evaluations, reuse buffers:

```@example ads
result = zeros(noutputs(map))
normalized = zeros(nvariables(map))
work = zeros(degree(map) + 1)
evaluate!(result, map, [0.2, -0.1], normalized, work)
```

`max_depth` limits splits along a path; `max_patches` limits the total number of
leaves. By default, an unresolved patch at a limit raises an error. With
`strict = false`, the map still covers the entire input box, and each unresolved
patch has status `:max_depth`, `:max_patches` or `:roundoff`. Always check
`map.converged` before using such a result.

Construction restores the caller's algebra and existing polynomials, even after
an exception. The callback must not change the algebra configuration. Construction
uses the global algebra and must not overlap other polynomial calculations;
numeric evaluation of completed maps is independent of it.

The [Kepler example](@ref "Automatic domain splitting for a Kepler map") compares
a single expansion with ADS, validates both against scalar propagation, and plots
the domain partition and propagated uncertainty.
