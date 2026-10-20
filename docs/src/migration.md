# API migrations

## 0.2.0 (unreleased)

- Certified ADS retries analytic domain failures on child cells. A failure that
  remains at a splitting limit raises an `ErrorException` with its cause, including
  with `strict=false`; only valid models can form unresolved patches.
- Taylor-model scalar operands and explicit remainders reject implicit floating
  values as NG. Use `interval(value)`, `IntervalArithmetic.exact(value)`, or exact
  integers/rationals in model arithmetic and `IntervalBound` callbacks.

From `24dc1a1`:

- Compound analytic operations, integer powers (including negative literals),
  and polynomial exponentiation filter coefficients only after completing the
  result. Global tolerance and truncation order are unchanged.
- Interval coefficient norms return outward-rounded enclosures, including for
  intervals containing zero, and preserve coefficient metadata. Ordinary-number
  norm types and stable accumulation are unchanged.
- Diagonal polynomial eigenpairs derive scalar working types without a dense
  eigensolve; BigFloat and Rational{BigInt} diagonals need no optional solver.

From `720b46c`:

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

From `17f74cc`:

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
- `float(p)` converts coefficients to match `float(typeof(p))`. Public storage
  constructors require concrete coefficient types. `TimeSeries(vector)` owns a
  vector of copied coefficients, including when given a view or range.
- `central_moments(mgf, order)` requires constant coefficient one and rejects orders
  above the working cutoff. Use `with_order` explicitly when needed; inspection
  never changes the global cutoff. Raw moments still inspect all stored orders.

The internal, test-only `filled` helper is removed. Documented coefficient
inspection, construction, and norm helpers remain available in `DifferentialAlgebra`.

Variable names supplied both positionally and by keyword now raise
`ArgumentError`. Polynomial displays separate variable factors with `*` when
any variable name has more than one character, and parenthesize the components
of `Complex{TaylorPolynomial}`. Mixed ForwardDiff dual/polynomial arithmetic
preserves coefficient partials without requiring explicit polynomial constants.

ADS refinement now inherits the original tolerances and construction options,
including resource limits. Supply keywords explicitly to change them.
