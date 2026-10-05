# Continuous ADS surrogates

Independent ADS expansions need not agree on a face. For optimization, use
[`continuous_map`](@ref) to construct an optional overlapping approximation of
the original static function. It supports ordinary box ADS, interval box ADS
and 2D oriented/polygon ADS. The original partition and its evaluation behavior
are preserved.

```julia
using DifferentialAlgebra, ForwardDiff
f(v) = exp(v[1] + v[2]) + (v[1] - v[2])^2
fit = adaptive_map(f, [-1., -1.], [1., 1.]; order=3, atol=1e-4)
smooth = continuous_map(f, fit; continuity=:c2, order=4)
smooth([0.2, -0.1])
ForwardDiff.gradient(smooth, [0.2, -0.1])
ForwardDiff.hessian(smooth, [0.2, -0.1])
```

| Option | Regularity of the blended surrogate | Compact taper on `0 ≤ t ≤ 1` |
|:--|:--|:--|
| `continuity=:c0` | Continuous values | `1-t` |
| `continuity=:c1` | Continuous values and first derivatives | `(1-t)^2*(1+2t)` |
| `continuity=:c2` (default) | Continuous values, first and second derivatives | `(1-t)^3*(1+3t+6t^2)` |

These guarantees concern the ideal real surrogate defined by exact rational
geometry and stored local coefficients. Numeric evaluation, including
ForwardDiff, has ordinary floating rounding. C1/C2 regularity does **not** imply
accurate gradients/Hessians of the original function. Enlarging supports can
make derivative accuracy worse; higher local order or a finer source partition
can help. No original-function derivative error is certified.

## Construction and continuity proof

In box coordinates, each old leaf is a core rectangle. For polygons, the core
is its bounding rectangle in the source's fixed projected frame `z=B*x`.
`overlap=1//4` adds a quarter of the **full core width** on either side,
clipped to the projected root rectangle. The callback is reevaluated on each
new support, using the native polynomial engine. The same frame is retained;
the inverse `B⁻¹` maps projected inputs to physical inputs. Supports may extend
beyond a physical polygon, so the callback must be valid on their entire covers.
Fixed box coordinates have constant weight one and zero normalized coordinate.

For each coordinate, set the weight to one inside its core, taper it to zero
across either overlap margin, and extend it by zero outside the support. Take
the product across coordinates to obtain `w_i`, and define

```math
\lambda_i(x)=\frac{w_i(x)}{\sum_j w_j(x)},\qquad
S(x)=\sum_i\lambda_i(x)P_i(x).
```

The tapers are nonnegative. For Ck, their first k derivatives vanish at both
endpoints; their values join the constant-one plateau and constant-zero exterior.
Their products therefore are Ck. Every physical point lies in an original
leaf, and hence in a core with weight one, so the denominator is at least one.
Normalization by this positive Ck function and multiplication by polynomials
preserve Ck regularity, including shared faces, T-junctions and polygon corners.
Continuity is relative to the declared physical domain, including one-sided
behavior at its boundary. The construction does not require independently
computed patch jets to match.

The result is a partition-of-unity function, generally a piecewise rational
function rather than a polynomial. `max_order(smooth)` gives its retained local
order; `degree(smooth)` fails explicitly. `blend_weights(smooth, point)` returns
the normalized weights for inspection. Public evaluation requires exactly one
finite scalar coordinate per variable and rejects points outside the physical
domain. ForwardDiff's optional extension supports nested dual queries for
gradients and Hessians. Loading DifferentialAlgebra alone loads neither
ForwardDiff nor IntervalArithmetic and does not initialize an algebra.

## Function accuracy and interval certificates

For ordinary sources, `smooth.error_estimate` reports the recomputed heuristic
overlap indicators. For interval sources, use:

```julia
using DifferentialAlgebra, IntervalArithmetic, ForwardDiff
f(v) = (v[1] + v[2])^2
fit = adaptive_map(f, fill(interval(-1,1),2);
                   estimator=IntervalBound(), splitter=:oriented,
                   order=1, atol=1//64)
smooth = continuous_map(f, fit; continuity=:c2, atol=1//16)
smooth([1//4,1//4])                 # Numeric surrogate S
smooth.error_bounds                 # Uniform original-function error f-S
enclose(smooth, [1//4,1//4])         # Original function f at this point
enclose(smooth)                      # Original function f on the full domain
```

Each fresh interval-coefficient Taylor model retains its domain and entire
absolute remainder. Stored coefficient midpoints define the **explicitly
requested numeric surrogate**; the model certificates retain coefficient widths
and all remainders. This does not introduce a midpoint-coefficient model backend
or use midpoint extraction in validated arithmetic.

On support i, the previously documented ADS error formula proves
`f(x)-P_i(x) ∈ E_i`, including retained coefficient uncertainty. Since
the normalized weights are nonnegative and sum to one,

```math
f(x)-S(x)=\sum_i\lambda_i(x)(f(x)-P_i(x))
\in\operatorname{hull}_i E_i.
```

This gives the uniform `smooth.error_bounds` separately for each output.
It bounds error relative to the ideal real surrogate, not floating evaluation
roundoff. `enclose(smooth, query)` instead encloses the **original function** by
hulling fresh certificates on all intersecting supports, keeping their full
remainders. It accepts physical points, interval subboxes and, for two variables,
convex subpolygons; it rejects wrong dimensions and outside-domain queries.
No continuity of the resulting interval endpoints is claimed.

Enlarged supports have new errors: an accepted source tolerance does not imply
the same tolerance for the blend. Supply `atol` to check the newly computed
bounds (or heuristic estimates for ordinary sources). Failure throws with
suggestions to refine the original map, increase local order or reduce overlap.
Omitting `atol` still reports errors but makes no new tolerance assertion.
`rtol` and automatic refinement of overlap fits are not implemented. Domain
failures on a complete overlap cover also throw.

The wrapper owns its compiled polynomials, certificates and geometry. `copy`
and domain/error accessors own their storage; changing a source map cannot
change the wrapper. Snapshots survive algebra reinitialization. Internal arrays
are read-only implementation data. Construction restores the caller's algebra;
do not overlap construction with other global-algebra calculations.
Numeric evaluation scans supports and allocates; it does not change the existing
ordinary DA buffer/allocation guarantees.

The [continuity example](../generated/ads_continuity.md) plots values, gradients,
Hessians and original-function interval bands. The
[optimization example](../generated/ads_optimization.md) compares box and
oriented partitions and uses a C2 objective with damped Newton steps. Its
function-value bounds do not prove a minimizer or certify its gradients.
`benchmark/continuous_ads.jl` measures construction, numeric evaluation,
AD derivatives and enclosure width separately, including six-variable box ADS.
