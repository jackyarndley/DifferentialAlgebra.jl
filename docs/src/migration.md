# API migrations

## 0.2.0 (unreleased)

Julia 1.12 or later is required. Upgrade older Julia installations before upgrading
the package. Optional integrations remain optional.

### Removed APIs

- `push_order!(n)` / `pop_order!()` → `with_order(n) do … end`, which restores
  the working order after an exception.
- `validated_adaptive_map(f, box)` →
  `adaptive_map(f, box; estimator=IntervalBound(), order=3, atol=1e-6, splitter=:width)`.
  These explicit options preserve the wrapper's defaults. The canonical defaults
  remain order five, tolerance `1e-8`, and `:tail`.
- `DifferentialAlgebra.erf(p)` → `SpecialFunctions.erf(p)`. Likewise import `erfc`,
  `gamma`, `loggamma`, `besselj`, `bessely`, `besseli`, `besselk`, `besselix`, and
  `besselkx` from SpecialFunctions; their polynomial methods remain supported.
- `round(p)` → `round(constant_term(p))` when a scalar result is intended.
- `trunc(p)` → `trunc(constant_term(p))` when a scalar result is intended.
  Julia's fallbacks do not provide polynomial rounding. `mod` keeps its semantics.

### Changed semantics

- `Float64(p)` → `Float64(constant_term(p))` to discard derivatives explicitly.
  Scalar conversion, including `Bool`, rejects nonconstant polynomials.
  Structural `isequal` and `hash` inspect all stored coefficients; ordinary `==`
  still compares constant parts.
- `TaylorPolynomial{Real}(c)` → `TaylorPolynomial{Float64}(c)` or another concrete
  coefficient type. `float(p)` now matches `float(typeof(p))`.
- `central_moments(mgf, n)` → `with_order(n) do; central_moments(mgf, n); end` if
  the requested order exceeds the working cutoff. The constant coefficient must
  be one. Raw moments still inspect all stored orders.

Working orders and coefficient tolerances belong to the owning algebra and reach
child tasks. Configure them before parallel calculations and await children inside
`with_order`. Distributed workers need explicit configuration. ADS uses private
scoped algebras; `initialize!` inside a construction replaces only that binding.
Global initialization still invalidates the previous default's polynomials.

Repeated derivatives and Hessians preserve intermediate degrees needed for the
requested result; both `abs` branches use the working cutoff. Norms use stable
floating accumulation for ordinary coefficients and outward bounds for intervals.
Diagonal eigenpairs need no optional scalar solver, including BigFloat and
Rational{BigInt}. Root and hypotenuse coefficient types follow their scalar
operations. NaN constants propagate through polynomial analytic functions.
`TimeSeries(vector)` owns its storage, including input views and ranges.

### New errors

- `coefficient(p, [1])` → `coefficient(p, [1, 0])` for a two-variable polynomial.
  Exponent and repeated-calculus vectors require one nonnegative count per variable;
  remove extra entries explicitly. Numeric evaluation keeps its documented padding.
- `TaylorMethod(20.0)` → `TaylorMethod(20)`; time orders must be integers.
- `ODEProblem(f, [1.0], span, 1+x)` → `ODEProblem(f, [one(x)], span, 1+x)` for
  uncertain polynomial parameters. Interval states and parameters remain unsupported.
- `variables((:x, :y); names=...)` → `variables((:x, :y))`: supply names once.

Unsupported axes fail before indexing. Eigenpairs require symmetry through the
working order. Converting `Monomial{T}(c, integer_exponents)` is supported again;
negative exponents and values beyond `UInt32` are rejected.

### ADS

- `adaptive_map(f, previous)` → `adaptive_map(f, previous; atol=…, max_depth=…)`
  when changing saved settings. Refinement now inherits tolerances and resource limits.
- `map.error_estimate` → `only(map.error_estimate)` for a scalar interval-polygon
  estimate; output estimates consistently use vectors.
- Repeated `map.patches[i]` → `patches = map.patches; patches[i]`: the public
  accessor copies the vector. Snapshot entries are read-only; copy them to own storage.
- Queries such as `[π]` → `[BigFloat(π)]` when choosing boundary precision
  explicitly. Irrational containment and routing now consistently use BigFloat.

Certified ADS subdivides on analytic domain failures, including callback
`DomainError`s. A failure remaining at a splitting limit throws even with
`strict=false`. Use another exception type for errors that should propagate at once.

### Interval models

- `x * 0.1` → `x * interval(0.1)` or `x * IntervalArithmetic.exact(0.1)` to declare
  the stored binary constant. Use `x * (1//10)` for an exact decimal rational.
- `TaylorModel(p, 0.0, reference)` → `TaylorModel(p, interval(0), reference)`.
  Polynomial coefficients have stored-binary semantics; the supplied remainder is
  an assertion, not a reconstruction of earlier rounding or truncation errors.

Queries intersect `P+R` with an independent range certificate when available.
Certificates retain the full normalized cover and propagate through arithmetic;
copies and compiled snapshots own them. They do not reduce the absolute remainder.
Wide elementary compositions may return a constant image model when its bound is
tighter or polynomial intermediates overflow. The uniform remainder invariant still
holds. Validated arithmetic requires correct interval rounding and preserves metadata.
