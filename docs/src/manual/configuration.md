# Configuration and storage

## Algebra lifetime

`variables(n; order)` creates a global algebra. The lower-level
`DifferentialAlgebra.init(order, n)` initializes the same configuration without
constructing variables. Reinitializing invalidates existing polynomials, even
when the dimensions are unchanged.

Set up the algebra before concurrent calculations. Operations may share read-only
inputs, but each task must own its outputs. Configuration changes and mutation
of shared polynomials must not run concurrently with calculations.

A polynomial owns mutable coefficient storage. Assignment shares that object;
`copy(p)` creates independent storage. A compiled map owns its coefficients, so
numeric evaluation continues to work after reinitialization. Polynomial
composition still requires the original algebra.

## Truncation and filtering

`DifferentialAlgebra.setTO(order)` changes the working order without changing the
maximum order chosen at initialization. It returns the previous working order.
`pushTO` and `popTO` provide a stack for temporary changes.

`DifferentialAlgebra.setEps(tolerance)` discards sufficiently small floating-point
coefficients. The default threshold is zero. Filtering and truncation are
approximations, not rigorous error bounds.

## Memory and repeated calculations

The `table_bytes` keyword of `variables` and `init` bounds the multiplication
lookup tables. Larger bases calculate indices as needed. This budget excludes
polynomial coefficients and basis metadata.

Reuse compiled maps when evaluating many points. The in-place operations
`evaluate!`, `LinearAlgebra.mul!`, `DifferentialAlgebra.add!` and
`DifferentialAlgebra.scale!` can reuse storage. Scalar types with allocating
arithmetic may still allocate inside these operations.

## Persistence

`DifferentialAlgebra.toString(p)` and `parse(DA{T}, text)` provide text round
trips at the coefficient type's precision. `write(io, p)` and `read(io, DA)`
interoperate with the DACE binary format. Binary records contain double-precision
coefficients and 32-bit exponent indices.

The text parser accepts DACE records and whitespace-separated COSY-style
records, including Fortran `D` exponents. Terms outside the initialized order
or variable count are discarded when reading.
