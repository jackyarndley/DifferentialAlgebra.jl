# Intervals and Taylor models

Load optional IntervalArithmetic (`1.0.12` or later `1.0.x`) for interval
coefficients and function enclosures. Validated operations require its default
`rounding=:correct`; other modes are rejected without changing the setting.
Do not change rounding or algebra settings during a calculation.

## Choose what to enclose

| Calculation | Meaning |
|:--|:--|
| `enclose(p, box)` | Range of the stored polynomial, including its stored coefficients |
| `TaylorPolynomial{Interval{T}}` | Truncated polynomial with interval coefficient arithmetic |
| `TaylorModel` | Original function represented by a polynomial plus an absolute remainder |

The order-one square makes the distinction explicit:

```@example interval_distinction
using DifferentialAlgebra, IntervalArithmetic
box = [interval(-1, 1)]
p, = variables(Interval{Float64}, 1; order=1)
enclose(p*p, box)              # [0,0]: degree-two terms were discarded
```

```@example interval_distinction
x, = taylor_models(box; order=1)
m = x*x
(polynomial(m), remainder(m), enclose(m))  # zero polynomial, [0,1], [0,1]
```

Interval coefficients enclose supported coefficient arithmetic, but do not bound
truncation. For example, enclosing an ordinary `sin(p)` polynomial does not enclose
the original sine function. Reevaluating an expression with model inputs propagates
its remainder at every operation.

## Build and query a model

```@example interval_models
using DifferentialAlgebra, IntervalArithmetic
box = [interval(Float64, -1//4, 1//4), interval(Float64, -1//8, 1//8)]
x, y = taylor_models(box; order=3, names=(:x, :y))
m = exp(x+y)*cos(x*y)
(whole_domain = enclose(m), point = m([0, 0]),
 subbox = enclose(m, [interval(Float64, 0, 1//4), interval(0)]))
```

`polynomial(m)` returns an independent polynomial in normalized coordinates;
`remainder(m)` returns the absolute interval error and `domain(m)` the physical
box. Restricting a query keeps the entire remainder. To improve it, reevaluate
the original expression on a smaller domain.

`compile(m)` creates an owned `CompiledTaylorModel` for numeric evaluation and
enclosure after algebra reinitialization. It does not support arithmetic or
composition. Live models follow `variables`' lifetime rules; `taylor_models`
reinitializes the algebra. Binary model arithmetic requires the same coordinate
family and retained order. If the working order changed, wrap arithmetic in
`with_order(max_order(m)) do ... end`; enclosure queries still inspect the stored model.

## Constants, domains and guarantee flags

Supported domains use finite, nonempty, guaranteed decorated intervals with
Float64 or BigFloat endpoints. NaI, empty, unbounded, loose (`trv`) and bare
intervals are unsupported. Interval metadata is preserved, including fixed and
unused coordinates; copying an interval never upgrades its guarantee or decoration.

An explicit floating interval denotes its stored binary value:

| Constant | Meaning |
|:--|:--|
| `interval(0.1)` or `IntervalArithmetic.exact(0.1)` | The binary Float64 value `0.1` |
| `interval(1//10)` or `I"0.1"` | An enclosure of the exact decimal real 1/10 |
| `x*(1//10)` in a model callback | Exact rational converted outward |
| `x*0.1` in a model callback | Implicit floating conversion; NG and rejected |

Use explicitly declared interval constants in validated callbacks. Ordinary
interval-coefficient polynomials preserve NG flags rather than reject every such
operation; model construction/arithmetic reject NG and invalid metadata.
See IntervalArithmetic's [construction](https://juliaintervals.github.io/IntervalArithmetic.jl/stable/manual/construction/)
and [guarantee rules](https://juliaintervals.github.io/IntervalArithmetic.jl/stable/manual/guarantee/).

`enclose(p, box)` uses the stored values of an existing polynomial. Likewise,
`TaylorModel(p, R, reference)` **asserts** the supplied absolute remainder in the
reference coordinates, copies its inputs and rejects excessive degree. Neither
operation reconstructs earlier rounding or discarded terms. Attaching a zero or
estimated remainder cannot certify an unknown generating function.

## Representation and range bounds

A model stores a polynomial and remainder satisfying

```math
f(c+r\odot\xi)\in P(\xi)+R,\qquad \xi\in[-1,1]^n.
```

Centers are stored binary midpoints; radii bound both endpoint distances outward.
The normalized cover can include extra physical points, so function-domain checks
must hold throughout that cover. Fixed coordinates have zero radius and normalized
value zero. Physical queries normalize with outward interval subtraction/division.

Range bounding sums interval monomials, caching coordinate powers. Even powers
retain their sign information (`[-1,1]^2=[0,1]`), but repeated-coordinate dependency
can still widen bounds. This does not alter the polynomial or working order.
`DifferentialAlgebra.bounds(p)` instead uses ordinary arithmetic on an implicit
unit box; wrapping those floating endpoints in an interval adds no certificate.

Elementary operations also keep an independent interval image of the complete
input enclosure. Queries intersect it with `P+R`, and arithmetic propagates it
outward. This can retain range facts such as `sin(X) ⊆ [-1,1]`; it does not reduce
`R`, certify derivatives or change ADS errors. Explicit scalar/model constructors
do not inherit the reference's generating-function certificate.

Models, copies and snapshots own coefficients, remainders and certificates,
including BigFloat endpoints. Accessors return independent data. Assignment shares
a model; internal arrays are read-only and no public in-place mutation is provided.

## Supported arithmetic and remainder bounds

Models support constants, coordinates, unary signs, addition, subtraction,
multiplication, integer powers, division/reciprocal, `exp`, `log`, `sin`, `cos`,
`sincos` and `sqrt`. Retained coefficients round outward; tolerance pruning and
midpoint arithmetic are not used. Interval-coefficient polynomials additionally
support their ordinary calculus, substitution and composition operations, with
coefficient arithmetic rather than function-remainder guarantees.

For `A=P+R_A` and `B=Q+R_B`, multiplication retains the product through the
recorded degree and bounds every discarded ordered product:

```math
R_{AB}=\operatorname{bound}(\text{discarded products})
+\operatorname{bound}(P)R_B+\operatorname{bound}(Q)R_A+R_AR_B.
```

Unit-box monomial ranges use exponent parity and fixed-axis zeros; general boxes
use cached powers. Both multiplication backends share this invariant. Tail
enumeration is quadratic in the numbers of active terms.

Elementary functions first enclose the whole input by `X`, choose a thin center
`s` within it, then compose with remainder-aware Horner arithmetic. Exponential,
sine and cosine add the Taylor-theorem tail

```math
\frac{f^{(n+1)}(X)}{(n+1)!}(X-s)^{n+1}.
```

Reciprocal and logarithm use `U=(A-s)/s` to avoid overflowing powers of `s`.
With the interval counterpart `U_X=(X-s)/s`, their additional tails satisfy

```math
T_{\mathrm{inv}}=\frac{(-U_X)^{n+1}}{X},\qquad
T_{\log}\subseteq\frac{(-1)^n U_X^{n+1}}{(n+1)(X/s)}.
```

If the logarithm tail overflows or its rounded denominator contains zero, an
independent interval residual `log(X)-log(s)-P_n(U_X)` supplies the bound.
For square root, choose `s ≥ sup(X)/2`, so `|U_X| ≤ r ≤ 1`, and use

```math
|T_{\sqrt{\cdot}}|\le\sqrt{s}\,\frac{\binom{2n}{n}}{4^n}\,r^{n+1}.
```

For `r<1`, the first omitted binomial coefficient divided by `1-r` supplies an
additional bound; use the smaller. The tail is nonpositive for odd `n`.
All integer/factorial factors are enclosed. Incoming remainders and rounding are
already included by model arithmetic. Samples or fitted decay play no role.

Reciprocal requires `0 ∉ X`; log and nonconstant square root require `inf(X)>0`.
An identically zero enclosure has an exact square root. A regular center alone
is insufficient. A wide enclosure may fail despite a valid original function;
refine and reevaluate. When elementary polynomial intermediates overflow or a
constant image yields a tighter remainder, the result can fall back to that
constant model while preserving the invariant. Otherwise nonfinite results fail.

## Certified domain splitting

```@example interval_models
f(v) = exp(v[1]+v[2])*cos(v[1]*v[2])
fit = adaptive_map(f, box; estimator=IntervalBound(), order=3, atol=1e-5)
fine = adaptive_map(f, fit; atol=1e-6)
(converged = fine.converged, patches = length(fine.patches), range = enclose(fine))
```

Each child reevaluates the deterministic callback with fresh models. For a stored
midpoint `m_k` of coefficient interval `C_k`, its uniform fit error is

```math
E=R+\sum_k(C_k-[m_k,m_k])\,\operatorname{bound}(\xi^{\alpha_k}).
```

A patch passes only when `sup(abs(E)) ≤ atol` for every output. This includes
coefficient uncertainty; checking `R` alone would not. Midpoints define the
comparison polynomial only. Relative tolerances are unsupported; use scalar or
per-output absolute tolerances.

Analytic domain failures, including callback `DomainError`s, trigger subdivision.
Use a different exception for callback errors that should propagate immediately.
If no valid model exists at a splitting/roundoff limit, construction throws even
with `strict=false`. Otherwise `strict=false` retains valid unresolved models
with their status. The [ADS manual](domain-splitting.md) defines refinement,
direction scores and resource limits.

Queries intersect relevant closed patches and hull their enclosures. A split-face
point visits both children; a straddling subbox visits both branches. Snapshots
survive algebra changes, and public accessors own their returned storage.

## Non-axis-aligned domains

`splitter=:oriented` or a `ConvexPolygon` selects 2D polygon ADS. A nonsingular
projection matrix `B` defines `z=B*x`; the exact inverse maps back to physical
coordinates. Models enclose the full projected parallelogram using outward
interval endpoints. Exact clipping preserves the partition and shared boundaries.

`enclose(fit, subpolygon)` accepts convex subpolygons; interval subboxes must lie
inside the physical root polygon. All intersecting patches contribute, including
closed faces. Automatic frame selection is heuristic and does not affect the
certificate. Exact geometry can cost more even when it reduces patch counts.
See [polygon usage](domain-splitting.md#Oriented-and-polygonal-splitting).

## Limits

No derivative remainder, verified inversion or validated time integration is
provided. Whole-model comparisons, scalar conversion, model calculus and
model-vector compilation are unsupported. `TaylorMethod`, `taylor_expand` and
`adaptive_flow` retain numerical/heuristic error contracts.

The [interval tutorial](../generated/interval_models.md) plots the three distinct
range targets and the effect of fresh child fits. The remainder invariant follows
the [Taylor-model definition](https://juliaintervals.github.io/TaylorModels.jl/stable/)
using this package's native polynomial engine.
