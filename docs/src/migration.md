# API migrations

## 0.2.0 (unreleased)

- Julia 1.12 or later is required; update older Julia installations before upgrading.

- Working orders and filtering tolerances belong to the owning algebra and reach
  child tasks and threaded loops. Configure settings before parallel calculations;
  await child tasks inside `with_order`. Configure Distributed workers explicitly.
  ADS construction uses private algebras whose bindings are inherited by child
  tasks. Global `initialize!` still invalidates the old default's polynomials;
  inside a temporary construction it replaces only that scope's algebra.

- Certified ADS retries analytic domain failures on child cells. A failure that
  remains at a splitting limit raises an `ErrorException` with its cause, including
  with `strict=false`; only valid models can form unresolved patches.
- Taylor-model scalar operands and explicit remainders reject implicit floating
  values as NG. Use `interval(value)`, `IntervalArithmetic.exact(value)`, or exact
  integers/rationals in model arithmetic and `IntervalBound` callbacks.
- Taylor-model queries intersect the polynomial-plus-remainder bound with an
  independent function-range certificate when available. The stored absolute
  remainder and ADS uniform fit-error contract are unchanged; compiled snapshots
  retain the range certificate after algebra reinitialization.

- Compound analytic operations, integer powers (including negative literals),
  and polynomial exponentiation filter coefficients only after completing the
  result. The algebra's tolerance and truncation order are unchanged.
- Interval coefficient norms require correct interval rounding and return
  outward-rounded enclosures, including for
  intervals containing zero, and preserve coefficient metadata. Ordinary-number
  norm types and stable accumulation are unchanged.
- Diagonal polynomial eigenpairs derive scalar working types without a dense
  eigensolve; BigFloat and Rational{BigInt} diagonals need no optional solver.
- Polynomial `sqrt`, `cbrt`, and `hypot` use their scalar operation's working
  coefficient type, including zero roots and axis-aligned hypotenuses with
  integer or rational coefficients. `nthroot(p, 1)` still copies all stored terms
  with their original type.

- Import `erf`, `erfc`, `gamma`, `loggamma`, `besselj`, `bessely`, `besseli`,
  `besselk`, `besselix`, and `besselkx` from SpecialFunctions. Their polynomial
  methods remain supported; the redundant `DifferentialAlgebra` aliases are removed.
- Coefficient norms return floating norm types, including for integer coefficients.
  Converting `Monomial{T}(coefficient, integer_exponents)` constructors are restored;
  negative exponents and values beyond `UInt32` are rejected.
- `eigenpairs` requires symmetry for scalar and polynomial matrices, inspecting
  coefficients through the working order. Polynomial outputs use the scalar
  eigensolver's types.
- Horner composition, normalization, quotient construction, and repeated
  derivatives preserve intermediate coefficients until filtering the completed result.
  Repeated derivatives and Hessians retain intermediate degrees needed
  to produce the requested result at the working order. Both `abs` branches
  discard terms above that order.

- Replace `push_order!`/`pop_order!` with exception-safe `with_order(order) do … end`.
- Exponent and repeated calculus count vectors now require exactly one entry per
  variable, with nonnegative machine-sized counts. Pad omitted exponents with zeros
  and remove unused entries explicitly. Coefficient exponents must fit the basis.
  Unsupported array axes fail before indexing. Numeric evaluation retains its
  documented zero padding and ignored extra coordinates.
- Replace `validated_adaptive_map(f, box; …)` with
  `adaptive_map(f, box; estimator=IntervalBound(), order=3, atol=1e-6,
  splitter=:width, …)` to preserve the old defaults. The canonical defaults remain
  order five, tolerance `1e-8`, and `:tail`.
- Derivative-preserving `round(p)` and `trunc(p)` are removed. Round
  `constant_term(p)` explicitly when a scalar result is intended. Julia's generic
  fallbacks do not provide polynomial rounding. `mod` retains its existing semantics.
- Scalar conversion rejects nonconstant polynomials with `InexactError`; use
  `constant_term(p)` to discard derivatives explicitly. Structural `isequal` and
  `hash` include stored coefficients; ordinary `==` uses the constant part.
- `float(p)` converts coefficients to match `float(typeof(p))`. Public storage
  constructors require concrete coefficient types. `TimeSeries(vector)` owns a
  vector of copied coefficients, including when given a view or range.
- `central_moments(mgf, order)` requires constant coefficient one and rejects orders
  above the working cutoff. Use `with_order` explicitly when needed; inspection
  never changes the algebra's cutoff. Raw moments still inspect all stored orders.

The internal, test-only `filled` helper is removed. Documented coefficient
inspection, construction, and norm helpers remain available in `DifferentialAlgebra`.

Variable names supplied both positionally and by keyword now raise
`ArgumentError`. Polynomial displays separate variable factors with `*` when
any variable name has more than one character, and parenthesize the components
of `Complex{TaylorPolynomial}`. Mixed ForwardDiff dual/polynomial arithmetic
preserves coefficient partials without requiring explicit polynomial constants.

ADS refinement now inherits the original tolerances and construction options,
including resource limits. Supply keywords explicitly to change them.

ADS patch collections, Taylor-model patch outputs and error bounds use vectors.
Public `patches` and `models` accessors return independent containers; their
snapshot entries remain read-only. Use `copy(snapshot)` for independent snapshot
storage. Interval polygon `error_estimate` now returns a vector consistently
with ordinary polygon estimates.

The API catalog distinguishes supported qualified operations from internal
coefficient and buffer helpers. Internal methods remain available for package
implementation and tests; their signatures are documented in source.
