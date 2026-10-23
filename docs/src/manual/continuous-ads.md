# Continuous ADS surrogates

Independent ADS expansions can disagree at shared faces. For a smooth objective,
reevaluate the original function on overlapping supports and blend the local fits
with [`continuous_map`](@ref). Ordinary boxes, certified boxes and 2D polygons
are supported; the source partition is unchanged.

```@example continuous_ordinary
using DifferentialAlgebra, ForwardDiff
f(v) = exp(v[1]+v[2]) + (v[1]-v[2])^2
fit = adaptive_map(f, [-1.0, -1.0], [1.0, 1.0]; order=3, atol=1e-4)
smooth = continuous_map(f, fit; continuity=:c2, order=4)
(value = smooth([0.2, -0.1]),
 gradient = ForwardDiff.gradient(smooth, [0.2, -0.1]),
 hessian = ForwardDiff.hessian(smooth, [0.2, -0.1]))
```

| Setting | Surrogate regularity | Taper on `0 ≤ t ≤ 1` |
|:--|:--|:--|
| `:c0` | Continuous values | `1-t` |
| `:c1` | Continuous values and gradients | `(1-t)^2*(1+2t)` |
| `:c2` (default) | Continuous values, gradients and Hessians | `(1-t)^3*(1+3t+6t^2)` |

Smoothness concerns the ideal real surrogate defined by stored coefficients.
Floating evaluation has rounding error. C1/C2 alone says nothing about the
accuracy of derivatives of the original function; compare them independently.

## Construction and continuity

The original subdomains supply core rectangles, projected into the fixed frame
for polygon ADS. `overlap=1//4` adds one quarter of each full core width on both
sides, clipped to the projected root rectangle. The callback must be valid on
every complete support, which may extend outside a physical polygon.
Fixed coordinates have weight one.

Each weight is one in its core, tapers to zero across the overlap and is zero
outside its support. Coordinate weights multiply to give `w_i`. Normalize them:

```math
\lambda_i(x)=\frac{w_i(x)}{\sum_jw_j(x)},\qquad
S(x)=\sum_i\lambda_i(x)P_i(x).
```

The Ck taper's first k derivatives vanish at its endpoints. Every physical point
lies in a core with weight one, so the denominator is at least one. Products,
normalization and polynomial blending preserve Ck at faces, T-junctions and
corners, relative to the declared domain. Independently computed jets need not match.

The blend is generally piecewise rational. `max_order(smooth)` reports local
polynomial order; `degree(smooth)` is undefined. `blend_weights(smooth, point)`
returns all normalized weights. Queries require exactly one finite coordinate
per variable inside the physical domain; ForwardDiff supports nested duals.

## Errors on the enlarged supports

Fresh overlap fits have fresh errors. An accepted source tolerance is **not
inherited**. Use `atol` to check the new errors; without it, construction only
reports them. On failure, refine the source, increase local order or reduce
overlap. Automatic overlap refinement and `rtol` are unsupported.

Ordinary sources report heuristic `smooth.error_estimate`. Certified sources
retain original-function certificates:

```@example continuous_certified
using DifferentialAlgebra, IntervalArithmetic
f(v) = (v[1]+v[2])^2
fit = adaptive_map(f, fill(interval(-1, 1), 2);
    estimator=IntervalBound(), splitter=:oriented, order=1, atol=1//64)
smooth = continuous_map(f, fit; continuity=:c2, atol=1//16)
(surrogate = smooth([1//4, 1//4]), error = smooth.error_bounds,
 original_function = enclose(smooth, [1//4, 1//4]))
```

Coefficient midpoints define the explicitly requested numeric surrogate; the
certificates retain coefficient widths and full remainders. If
`f(x)-P_i(x) ∈ E_i` on support i, nonnegative weights summing to one give

```math
f(x)-S(x)=\sum_i\lambda_i(x)(f(x)-P_i(x))\in\operatorname{hull}_i E_i.
```

`smooth.error_bounds` encloses this ideal-surrogate error per output, excluding
floating evaluation roundoff. `enclose(smooth, query)` instead hulls the fresh
certificates to enclose the original function. It accepts points, interval
subboxes and, in 2D, convex subpolygons. Interval endpoints need not be continuous.
Function bounds certify neither derivatives nor a minimizer.

## Storage and examples

The wrapper owns compiled fits, geometry and certificates, survives algebra
reinitialization and copies data returned by domain/error accessors. Internal
arrays are read-only. Construction uses a private scoped algebra and restores
the caller's binding on failure. Query trees exclude distant supports; evaluation
cost depends on supports visited, while `blend_weights` returns a full vector.

The [continuity tutorial](../generated/ads_continuity.md) compares values,
gradients and Hessians across joins. The [optimization tutorial](../generated/ads_optimization.md)
compares box/polygon partitions and damped Newton convergence against an
analytical minimizer, with function certificates shown separately.
