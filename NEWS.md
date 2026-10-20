# Release notes

## 0.2.0 (unreleased)

This release removes legacy APIs and tightens numerical and ownership contracts.
See the [migration guide](docs/src/migration.md) for replacement calls and defaults.

- Temporary construction algebras, working orders and coefficient tolerances are
  task-local. Independent ADS builds may overlap read-only polynomial calculations.
  Global initialization still invalidates polynomials of the old default algebra.
- Scalar conversion rejects nonconstant polynomials. Structural `isequal` and
  `hash` include stored coefficients; ordinary comparisons retain their documented
  constant-part semantics. Unsigned subtraction and reciprocal expansions preserve
  signs through suitable working types.
- Analytic composition filters completed results. Repeated derivatives and Hessians
  preserve intermediate degrees needed at the requested working order. Mixed
  BigFloat reciprocal and power calculations normalize at the promoted precision.
- Bessel J/I expansions remain stable near zero; scaled K expansions remain stable
  at large positive centers. Boundary inverse functions, large/small hypotenuses
  and real-domain failures have independent regression coverage.
- Coefficient norms use stable ordinary-number accumulation and outward interval
  bounds with metadata preserved. Diagonal eigenpairs avoid an unnecessary dense
  solve, including BigFloat and Rational{BigInt} without GenericLinearAlgebra.
- Taylor models preserve decorations and guarantees. Implicit floating scalar
  operands are rejected; use interval construction or `IntervalArithmetic.exact`.
  Validated operations require correct interval rounding.
- ADS excludes affine variation from tail estimates, retries analytic domain
  failures on children, preserves refinement options, and reports roundoff limits.
  Shared box/polygon/continuous construction retains geometry and estimator adapters.
- Compiled snapshots own their coefficients and use concrete patch vectors without
  retaining obsolete multiplication tables. Conservative geometry filters and
  support trees reduce repeated exact routing and blending work.
- ODE mutation, dense interpolation and discontinuity handling use the public
  solver interfaces while preserving buffer ownership and coefficient precision.
- Tests isolate subjects with SafeTestsets, capture fresh-process failures, and
  exercise independent scalar/combinatorial oracles and both multiplication paths.
  CI covers Julia 1.10, current stable on Linux/Windows/macOS, and prereleases.

Julia 1.10 and optional integrations remain supported. SciML integration environments
require OrdinaryDiffEqCore 4.18, DiffEqBase 7.21.2 and SciMLBase 3.53 or later compatible
releases. Standalone ForwardDiff supports 0.10.36 and 1. IntervalArithmetic is limited
to compatible 1.0 patches from 1.0.12 while its rounding configuration is checked.
