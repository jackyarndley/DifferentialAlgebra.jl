# Release notes

## 0.2.0 (unreleased)

This release removes legacy APIs and tightens numerical and ownership contracts.
See the [migration guide](docs/src/migration.md) for replacement calls and defaults.

Julia 1.12 or later is required. The optional ForwardDiff integration requires
ForwardDiff 1; support for 0.10 has been removed.

- Temporary construction algebras are scoped and inherited by child tasks.
  Working orders and coefficient tolerances belong to their owning algebra.
  Configure them before parallel calculations; independent ADS builds retain
  private settings. Global initialization still invalidates the old default.
- Scalar conversion rejects nonconstant polynomials. Structural `isequal` and
  `hash` include stored coefficients; ordinary comparisons retain their documented
  constant-part semantics. Unsigned subtraction and reciprocal expansions preserve
  signs through suitable working types.
- Analytic composition, map inversion, eigenpair lifting and moment centering
  filter completed results. Repeated derivatives and Hessians preserve intermediate
  degrees needed at the requested working order. Scalar division avoids overflowing
  reciprocals; shared recurrences preserve finite extreme-scale coefficients.
  Mixed BigFloat calculations retain the promoted precision.
- Bessel J/I expansions remain stable near zero; scaled K expansions remain stable
  at large positive centers. Boundary inverse functions, large/small hypotenuses
  and real-domain failures have independent regression coverage. Narrow and unsigned
  integer Bessel orders are normalized before recurrence arithmetic.
- Coefficient norms use stable ordinary-number accumulation and outward interval
  bounds with metadata preserved. Diagonal eigenpairs avoid an unnecessary dense
  solve, including BigFloat and Rational{BigInt} without GenericLinearAlgebra.
- Raw moments scale coefficients by small factorial factors, avoiding premature
  overflow and `0*Inf`. Exact arithmetic stays exact; fixed-width integer overflow
  is explicit. Unsigned polynomial bounds preserve negative endpoints.
- Taylor models preserve decorations and guarantees. Implicit floating scalar
  operands are rejected; use interval construction or `IntervalArithmetic.exact`.
  Zero powers and zero square roots retain input metadata. Certified point queries
  enclose irrational coordinates outward and preserve stored BigFloat precision.
  Validated operations require correct interval rounding.
- ADS excludes affine variation from tail estimates, retries analytic domain
  failures on children, preserves refinement options, and reports roundoff limits.
  Shared box/polygon/continuous construction retains geometry and estimator adapters.
- Polygon frames normalize direction rows without changing physical cuts.
  Exact affine-cover construction avoids overflow before cancellation and rejects
  covers that the source scalar type cannot represent.
- Compiled snapshots own their coefficients and use concrete patch vectors without
  retaining obsolete multiplication tables. Conservative geometry filters and
  support trees reduce repeated exact routing and blending work.
- ODE mutation, dense interpolation and discontinuity handling use the public
  solver interfaces while preserving buffer ownership and coefficient precision.
  Adaptive-flow callbacks may reuse scratch vectors; retained trajectories own
  their state. Scalar `save_idxs` gives scalar dense values and derivatives.
- Eager and recorded time expansions share numerical kernels and preserve exponent
  sensitivities carried by ForwardDiff. Zero-denominator and nonanalytic-center
  checks are consistent across their arithmetic paths. Polynomial seeds and time
  intermediates retain coefficients until the public calculation is complete.
- Orbital tutorials share Kepler dynamics, element conversion and complete
  polynomial moment calculations, with independent analytical checks.

Optional integrations remain supported. SciML integration environments
require OrdinaryDiffEqCore 4.18, DiffEqBase 7.21.2 and SciMLBase 3.53 or later compatible
releases. IntervalArithmetic is limited to compatible 1.0 patches from 1.0.12
while its rounding configuration is checked.
