# Intervals and Taylor models

Install and load the optional IntervalArithmetic dependency (compatible series
`1.0.12` and later `1.x`). DifferentialAlgebra alone neither loads it nor
initializes an algebra. TaylorSeries and TaylorModels are not dependencies.

```julia
using DifferentialAlgebra, IntervalArithmetic
```

Validated enclosure and model operations require the default
`IntervalArithmetic.configure(rounding=:correct)` policy. Other rounding modes
are rejected, including non-rigorous `:none`; the initial implementation also
leaves `:ulp` unsupported. The extension checks the setting without changing it.
Do not change interval rounding or algebra configuration during a calculation.

## Three distinct calculations

`enclose(p, box)` rigorously encloses the **stored polynomial** on an explicit
box of decorated intervals. It does not recover terms previously truncated by
DA or bound an unknown generating function. For example, `enclose(sin(x), box)`
encloses the computed finite Taylor polynomial of sine, not sine itself.

`variables(Interval{Float64}, n; order)` gives polynomials with interval-valued
coefficients. Construction, promotion, sums, products, division at regular
constant intervals, integer powers, calculus, substitution, composition,
evaluation, `exp`, `log`, `sin`, `cos`, and `sqrt` support coefficient arithmetic.
These polynomials still have ordinary truncation semantics and no function
remainder. Other coefficient operations depend on scalar methods and may fail;
there is no midpoint fallback. Only coefficients proved to be valid, guaranteed
exact zero are skipped. An interval containing zero is retained. Approximate
coefficient tolerance never filters interval coefficients.

Taylor models retain an absolute remainder and a validity domain:

```julia
box = [interval(Float64, -1//4, 1//4), interval(Float64, -1//8, 1//8)]
x, y = taylor_models(box; order = 3, names = (:x, :y))
m = exp(x + y) * cos(x * y)
enclose(m)                       # Uniform enclosure of the original expression
evaluate(m, [0, 0])              # Interval enclosure at a physical point
enclose(m, [interval(Float64, 0, 1//4), interval(0)])
polynomial(m)                    # Independent copy in normalized coordinates
remainder(m)                     # Absolute remainder, including truncation
domain(m)                        # Independent copy of the physical box
```

## Input semantics and guarantees

Float64 and BigFloat interval endpoints are supported. BareInterval is not a
real coefficient type and is deliberately unsupported. Boxes must be finite,
nonempty decorated intervals of the correct dimension. NaI, empty, unbounded
and loose-evaluation (`trv`) intervals fail explicitly.

An ordinary floating value denotes its stored binary value. In particular,
`interval(0.1)` is a thin enclosure of the binary Float64 `0.1`, which is greater
than the decimal real `1/10`. Use `interval(1//10)` or `I"0.1"` to enclose that
decimal real. Already-rounded ordinary coefficients cannot recover earlier
rounding or truncation errors.

The explicit operations `enclose` and `TaylorModel(value, reference)` use those
declared stored-value semantics. Exact integers, rationals and internal degree
factors use explicit interval construction. Ordinary mixed float arithmetic on
interval-coefficient polynomials follows IntervalArithmetic's implicit conversion
semantics and may produce NG; for example `0.1 * variable(1, Interval{Float64})`.
Use an explicitly constructed interval scalar to keep its guarantee.
`TaylorPolynomial{Interval{Float64}}(0.1)` also preserves the NG conversion flag.
Converting or reconstructing an existing interval never repairs NG or upgrades
its decoration. Polynomial enclosure preserves flags; model construction and
arithmetic reject NG or invalid/degraded data.

See IntervalArithmetic's [guarantee rules](https://juliaintervals.github.io/IntervalArithmetic.jl/stable/manual/guarantee/)
and [interval construction](https://juliaintervals.github.io/IntervalArithmetic.jl/stable/manual/construction/).
These rules are distinct from a proof of uniform truncation error.

## Range bounding

The initial range bounder sums interval monomials with outward rounding. It uses
the existing exponent matrix and caches integer powers of each coordinate.
The public `pown` operation preserves even powers: `[-1,1]^2` is bounded by
`[0,1]`, not `[-1,1]`. Sums can overestimate because repeated occurrences of
the same coordinate are evaluated independently. The bounder does not change
the algebra, its working order or the polynomial.

The older `DifferentialAlgebra.bounds(p)` uses ordinary accumulation and assumes
`[-1,1]^n` implicitly. Its floating endpoints are not certified. Wrapping those
endpoints in an interval does not make their prior calculation rigorous.

## Representation, normalization and ownership

For physical coordinates `c + r .* ξ`, a model stores `(P, R)` with

```math
f(c+r\odot\xi)\in P(\xi)+R.
```

Midpoints are chosen as explicit stored binary values. Each radius is the upper
endpoint of an outward-rounded enclosure of both distances to the physical
endpoints. This conservatively covers the physical domain by `[-1,1]^n`.
The symmetric normalized box may cover extra physical points when a midpoint
rounds; domain validation on that larger box is conservative. Fixed coordinates
have zero radius and normalized value zero. Physical queries normalize with
interval subtraction and division and intersect with the declared normalized
box. Subboxes keep the entire original remainder.

Each coordinate family has an identity, domain, normalization, algebra owner
and retained order. Binary arithmetic requires the same family and order.
`taylor_models` follows `variables`' algebra-lifetime rules: reinitialization
invalidates earlier polynomials and models. Model arithmetic rejects a changed
working order rather than silently reducing the polynomial. Restore the recorded
order with `with_order(max_order(model))`. Enclosure and evaluation still use the
stored model when only the working order has changed.

Model constructors and `copy` own coefficients. The public accessors return
independent copies, including BigFloat endpoints. Assignment shares a model;
internal fields are read-only implementation data. No supported mutation API
can modify a model's polynomial, remainder or normalization in place.

`TaylorModel(p, R, reference)` explicitly asserts a supplied absolute error in
the reference coordinates, copies its inputs, and rejects terms above the
reference order. It represents the stored polynomial plus that asserted error.
Attaching zero or an estimated remainder does **not** certify the generating
function of an existing truncated DA polynomial. To obtain a function enclosure,
reevaluate the original expression using model inputs.

## Arithmetic and its justification

Constants, independent coordinates, unary signs, addition, subtraction,
multiplication, integer powers, reciprocal/division, `exp`, `log`, `sin`, `cos`,
and `sqrt` are supported. Retained coefficients use outward-rounded intervals.
No tolerance pruning or midpoint extraction is used.

For `A=P+R_A` and `B=Q+R_B`, multiplication retains the polynomial product up to
the recorded total degree and encloses

```math
R_{AB}=\operatorname{bound}(\text{discarded products})
 +\operatorname{bound}(P)R_B+\operatorname{bound}(Q)R_A+R_AR_B.
```

Every ordered coefficient product above that degree is bounded on the normalized
box, including mixed-variable terms. This enumerates the tail using the original
basis metadata and powers up to twice the order, without constructing another
algebra. Retained multiplication uses the existing table-backed or table-free
kernel. This reference tail algorithm is quadratic in the active term counts.

Elementary functions first bound the **whole** input model by `X`, choose a thin
binary expansion point `s` inside `X`, then evaluate the finite scalar Taylor
polynomial in `A-s` using remainder-aware Horner arithmetic. Taylor's theorem
adds the absolute interval

```math
\frac{f^{(n+1)}(X)}{(n+1)!}(X-s)^{n+1}.
```

The entire segment from `s` to each possible input lies in `X`. Normalized
derivatives use interval formulas: exponential divided by factorial, the sine/
cosine derivative cycle divided by factorial, `(-1)^k/X^{k+1}` for reciprocal,
`(-1)^{k-1}/(kX^k)` for logarithm, and
`binomial(1/2,k)*sqrt(X)/X^k` for square root. Integer and factorial factors are
constructed or accumulated as exact enclosed constants. This also propagates
incoming remainders and coefficient rounding. No fitted decay, sample discrepancy
or epsilon padding is part of the proof.

Reciprocal must exclude zero in `X`; log and the implemented square-root expansion
require `inf(X)>0`. The identically zero enclosure has an exact square root.
Other square-root inputs touching zero fail. A valid constant term is insufficient
when the full enclosure violates these conditions; loose interval evaluation
does not establish validity. Wide enclosures may reject a function that is
actually valid; use smaller input domains and reevaluate the original expression.
Finite stored coefficients and remainders are required; overflow fails clearly.

## Certified domain splitting

Select `estimator=IntervalBound()` in `adaptive_map` for the validated ADS
method. It reevaluates the original function with Taylor-model inputs on every
child and uses rigorous interval bounds to decide whether a fit needs splitting.
`validated_adaptive_map` is also available as a convenience entry point:

```julia
f(v) = exp(v[1] + v[2]) * cos(v[1] * v[2])
a = adaptive_map(f, box; estimator = IntervalBound(), order = 3, atol = 1e-5)
a.converged
length(a.patches)
evaluate(a, [0, 0])              # Enclosure, including the selected remainders
enclose(a)                       # Enclosure on the entire partition
enclose(a, [interval(Float64, 0, 1//4), interval(0)])
```

For each output, acceptance bounds the difference from the polynomial with
stored midpoint coefficients. If `m_k` is a stored midpoint and `C_k` its retained
interval coefficient, the uniform error bound is

```math
E=R+\sum_k (C_k-[m_k,m_k])\,\operatorname{bound}(\xi^{\alpha_k}).
```

All operations in this bound round outward. A patch is accepted only when
`sup(abs(E)) ≤ atol` for every output. Checking `R` alone would omit coefficient
uncertainty. Midpoints here define the comparison polynomial for acceptance;
there is no midpoint-coefficient model arithmetic backend. Scalar or per-output
absolute tolerances are supported; relative tolerances are not yet implemented.
The deterministic callback must return models or real constants and must not
change algebra settings.

On boxes, `splitter=:width` bisects the longest side relative to the original
box, skipping fixed or unsplittable coordinates. `:tail` selects a side using
retained coefficient sensitivity, with relative widths breaking ties. These
heuristic direction scores affect efficiency, not validity: only the uniform
interval error decides acceptance. `max_depth` and `max_patches` bound the work. Limits throw unless
`strict=false`, which retains valid unresolved models and records their status.
Failed full-enclosure function-domain checks still throw; the algorithm does
not assume validity at the center or automatically retry a failed expansion.

`compile(model)` produces an owned `CompiledTaylorModel` containing coefficients,
monomial exponents, normalization, domain and remainder. It supports numeric
enclosure/evaluation, not arithmetic or composition. Certified ADS stores these
snapshots in its patches, uses a temporary algebra during construction, and restores
the caller's configuration on success or failure. Returned maps and snapshots
remain usable after global algebra changes. Public domain/remainder accessors
return independent copies; internal arrays are read-only. Queries intersect
every overlapping patch and hull their enclosures with preserved interval flags,
so faces and subboxes spanning several patches are covered. Initial lookup is
linear in the number of patches.

GuardedTail, ExtrapolatedTail, LastTerms and `adaptive_flow`
retain their heuristic contracts. They do not acquire certified error bounds by
passing intervals to ordinary DA. Restricting an existing Taylor model can tighten
its polynomial range but keeps its full remainder. Certified splitting needs
fresh evaluation of the **original function** on both children.

## Non-axis-aligned domains

Use `splitter=:oriented` on a 2D box, or pass a `ConvexPolygon` directly to
`adaptive_map`. All four error methods share this geometry. Projection rows
come from an automatic sensitivity probe or an explicit nonsingular matrix
`directions=B`. The exact rational inverse `A=B⁻¹` defines `f(A*z)` for model
construction; a numerical transpose is never substituted for the inverse.
Each polygon is enclosed by its exact projected extents in `z=B*x`. Rational
endpoints and inverse entries are enclosed by outward-rounded intervals before
any polynomial arithmetic. The model proof holds over this whole parallelogram,
so it holds on the physical polygon contained within it.

Each split exactly clips the parent's closed polygon into two closed children.
Their union equals the parent; interiors are disjoint and the shared cut belongs
to both. Physical queries use exact rational containment. Query intersections
are projected and outwardly bounded, then evaluated using stored snapshots and
their entire remainders. Shared faces and degenerate point/line subboxes are
included. `enclose(map, subpolygon)` supports a convex subpolygon; interval
subboxes must lie entirely inside the original physical polygon.

Automatic directions use midpoint coefficients only for a heuristic choice of
frame. This does not enter certified coefficients, errors or acceptance. There
is no guarantee of an optimal direction or improved speed: fewer patches may
still cost more because exact geometry is more expensive. Frames stay fixed
during a construction. General nonconvex polygons, changing local frames,
higher-dimensional polyhedra, polygon-constrained range optimization and verified
flow propagation are not implemented. A parallelogram cover can be too wide
to establish a valid function domain; this throws rather than accepting a
valid center. See the [polygon example](../generated/polygon_ads.md) for plots
and the [ADS manual](domain-splitting.md#Oriented-and-polygonal-splitting) for usage.

## Remaining limits

Whole-model comparisons, scalar conversion, certified differentiation/integration,
verified map inversion and model-vector compilation are unsupported. An absolute
value remainder gives no derivative remainder. Independent time and uncertainty
orders and the existing time solver are unchanged. TaylorMethod, adaptive_flow
and taylor_expand are not validated integrators.

Uncertainty-domain validation does not bound time integration error. Certified
ADS currently supports deterministic static maps on boxes and 2D convex polygons
with absolute tolerances.

See the [Literate example](../generated/interval_models.md) and
`benchmark/interval_models.jl` for separate timing, allocation and width reports,
including a six-coordinate uncertainty map. Benchmark widths measure different
objects: stored-polynomial range versus original-function enclosure.

The remainder invariant follows the usual [Taylor-model definition](https://juliaintervals.github.io/TaylorModels.jl/dev/);
the implementation uses DifferentialAlgebra's native polynomial engine.
