# API migration from `17f74cc`

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
