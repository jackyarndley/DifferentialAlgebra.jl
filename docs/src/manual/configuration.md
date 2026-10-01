# Configuration and storage

## Algebra lifetime

`variables(n; order)` creates a global algebra. The lower-level
`DifferentialAlgebra.initialize!(order, n)` initializes the same configuration without
constructing variables. Reinitializing invalidates existing polynomials, even
when the dimensions are unchanged.

Set up the algebra before concurrent calculations. Operations may share read-only
inputs, but each task must own its outputs. Configuration changes and mutation
of shared polynomials must not run concurrently with calculations.

A polynomial owns mutable coefficient storage. Assignment shares that object;
`copy(p)` creates independent storage. A compiled map owns its coefficients, so
numeric evaluation continues to work after reinitialization. Polynomial
composition still requires the original algebra.

Julia manages polynomial memory automatically. A standard array copy, `copy(A)`,
shares its polynomial entries; use `copy.(A)` for independent coefficients.
When creating mutable work buffers, a comprehension such as
`[zero(p) for _ in 1:n]` creates independent entries. `fill(p, n)` shares `p`.

## Truncation and filtering

`DifferentialAlgebra.set_truncation_order!(order)` changes the working order without changing the
maximum order chosen at initialization. It returns the previous working order.
`with_order(order) do ... end` restores the previous order even when the calculation throws.

`DifferentialAlgebra.set_coefficient_tolerance!(tolerance)` discards sufficiently small floating-point
coefficients. The default threshold is zero. Filtering and truncation are
approximations, not rigorous error bounds.

## Memory and repeated calculations

The `table_bytes` keyword of `variables` and `initialize!` bounds the multiplication
lookup tables. Larger bases calculate indices as needed. This budget excludes
polynomial coefficients and basis metadata.

Reuse compiled maps when evaluating many points. The in-place operations
`evaluate!`, `LinearAlgebra.mul!`, `DifferentialAlgebra.add!` and
`DifferentialAlgebra.scale!` can reuse storage. Scalar types with allocating
arithmetic may still allocate inside these operations.

## Package compilation

PrecompileTools records a small, two-variable workload for common arithmetic,
elementary functions, differentiation and compiled-map evaluation. It uses
ordinary vectors and a single coefficient type. Polynomial order and variable
count are runtime data, so they do not create a separate family of types.
Other scalar types and application callbacks compile when first used.

Plotting packages and ODE solvers are optional application
dependencies; none is loaded or exercised by the package's precompile workload.
Loading DifferentialAlgebra leaves the algebra uninitialized.
