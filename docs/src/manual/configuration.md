# Configuration and storage

## Algebra lifetime

`variables(n; order)` creates the global default algebra;
`initialize!(order, n)` creates the same configuration without variables.
Reinitializing invalidates all live polynomials belonging to the previous
algebra, even when its dimension and order are unchanged.

Initialize before parallel work. Read-only operations can share polynomials;
mutation and reinitialization must not overlap readers. Each task owns its outputs.
Working orders and coefficient tolerances belong to the algebra and are shared
by child tasks. Configure them before parallel calculations and await children
inside `with_order`. Distributed workers need their own initialization.

ADS builders use private scoped algebras inherited by their child tasks. Other
tasks retain their default binding. Object-based operations always use the
object's owning algebra; logging scopes do not change that ownership.

## Working order and coefficient tolerance

```@example configuration
using DifferentialAlgebra
x, = variables(1; order=6)
p = exp(x)
low_order = with_order(2) do
    exp(x)
end
(stored_degree = degree(p), low_degree = degree(low_order), restored = truncation_order())
```

`set_truncation_order!(n)` changes the working cutoff within the initialized
maximum and returns the previous value. `with_order(n) do ... end` restores it
on success or failure. Changing the cutoff does not mutate existing polynomials.
Arithmetic results obey it, including scalar addition, unary signs and `p^1`.
Copies, coefficient conversion and inspection preserve stored terms.

`set_coefficient_tolerance!(tol)` discards small floating coefficients from
arithmetic results. The default is zero; a BigFloat tolerance keeps its precision.
Other tolerance types must convert to a finite Float64; pass a BigFloat explicitly
when a finite threshold exceeds that range.
Filtering and truncation are approximations, not error bounds. Taylor models
record their retained order and reject arithmetic under a different working
order; their queries preserve the stored polynomial and remainder.

## Ownership and repeated evaluation

| Operation | Storage behavior |
|:--|:--|
| `q = p` | Shares the mutable polynomial |
| `copy(p)` | Independent polynomial coefficients |
| `copy(A)` | New array, shared polynomial entries |
| `copy.(A)` | Independent polynomial entries |
| `[zero(p) for _ in 1:n]` | Independent work buffers |
| `fill(p, n)` | Every entry shares `p` |
| `CompiledMap(A)` | Owns coefficients; numeric evaluation survives reinitialization |

Reuse compiled maps for many evaluations. `evaluate!` and `LinearAlgebra.mul!`
reuse buffers, although scalar arithmetic may allocate. Polynomial composition
with a compiled map still requires the original algebra.

`table_bytes` bounds multiplication lookup tables in `variables`/`initialize!`.
Larger bases compute indices when needed. This budget excludes coefficients and
basis metadata; table-backed and table-free products have the same semantics.

## Package compilation

A small PrecompileTools workload covers Float64 arithmetic, elementary functions,
calculus and compiled evaluation. Variable count/order are runtime data; other
scalar types and application callbacks compile when first used. Loading the
package leaves its algebra uninitialized and loads no plotting or ODE package.
