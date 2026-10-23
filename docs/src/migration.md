# API migration

## 0.2.0 (unreleased)

Requires Julia 1.12 or later. Optional integrations remain optional; the
ForwardDiff integration now requires version 1. Deprecated wrappers and redundant
SpecialFunctions aliases have been removed.

## Replace removed calls

| Previous use | Replacement |
|:--|:--|
| `push_order!(n)` / `pop_order!()` | `with_order(n) do ... end` restores order on success and failure |
| `validated_adaptive_map(f, box)` | `adaptive_map(f, box; estimator=IntervalBound())` |
| `DifferentialAlgebra.erf(p)` and other special-function aliases | Import the operation from `SpecialFunctions`; polynomial methods remain available |
| `round(p)`, `trunc(p)` | Apply to `constant_term(p)` for a scalar result |
| `Float64(p)` to discard higher terms | `Float64(constant_term(p))` |
| `TaylorPolynomial{Real}(c)` | A concrete coefficient type, such as `TaylorPolynomial{Float64}(c)` |

The former validated wrapper used `order=3, atol=1e-6, splitter=:width`.
Specify these to preserve its behavior; the unified defaults are `order=5`,
`atol=1e-8`, `splitter=:tail`.

## Make inputs explicit

- Exponent/repeated-calculus vectors require one nonnegative count per variable:
  use `[1,0]` for a two-variable coefficient. Numeric evaluation retains its
  documented padding behavior.
- Time orders must be integers: `TaylorMethod(20)`.
- With uncertain polynomial ODE parameters, use polynomial initial states too:
  `ODEProblem(f, [one(x)], span, 1+x)`. Interval ODE states remain unsupported.
- Supply variable names once: `variables((:x,:y); order=5)`.
- Scalar conversion, including `Bool`, rejects nonconstant polynomials.
  Arithmetic comparisons still inspect constants; `isequal`/`hash` inspect all terms.
- Polynomial-map operations require one owning algebra, including zero-weight
  linear transforms and diagonal eigenproblems.
- `central_moments(mgf,n)` requires constant coefficient one and a sufficient
  working cutoff; use `with_order(n) do ... end` if needed. Raw moments inspect
  all stored orders.

Working order and coefficient tolerance belong to the owning algebra and reach
child tasks. Configure before parallel work and await children inside `with_order`.
Distributed workers require explicit initialization. ADS uses private scoped
algebras; global reinitialization still invalidates earlier default polynomials.
See [Configuration](manual/configuration.md).

## Adaptive maps and snapshots

Refinement now inherits the estimator, requested order, tolerances and resource
limits. Override only the settings to change:
`adaptive_map(f, previous; atol=..., max_depth=...)`.
Patch errors consistently use one entry per output; use `only(...)` for scalars.

`map.patches` copies its vector, so cache it before repeated indexing. Entries
are read-only; use `copy(map)` for independent patch storage. Irrational query coordinates
use BigFloat for routing; supply explicit BigFloat values when boundary precision
matters. See [Domain splitting](manual/domain-splitting.md).

## Validated arithmetic

Declare floating constants explicitly: `x*interval(0.1)` means the stored binary
value; `x*(1//10)` means the exact decimal rational. Implicit `x*0.1` is NG and
rejected by model arithmetic. Similarly, use an interval remainder in
`TaylorModel(p, interval(0), reference)`; that remainder is an assertion, not a
reconstruction of prior errors.

Certified ADS subdivides on analytic domain failures, including callback
`DomainError`s. A failure remaining at a limit throws even with `strict=false`;
use another exception type for errors that should propagate immediately.

Range queries can intersect `P+R` with an independent interval certificate.
Wide elementary compositions may return a constant image model when tighter or
needed to avoid overflowing intermediates. The absolute remainder invariant and
input metadata remain intact. See [Intervals and Taylor models](manual/interval-models.md).
