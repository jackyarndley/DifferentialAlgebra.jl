# Time expansions

```@meta
CurrentModule = DifferentialAlgebra
```

[`TimeSeries`](@ref) represents a local expansion in the time offset `τ = t-t₀`.
Its coefficients are ordinary powers, so coefficient `k` is the kth time
derivative divided by `k!`. Numeric coefficients need no global algebra.

```@example time
using DifferentialAlgebra
τ = TimeSeries([0.0, 1.0, zeros(6)...])
s = exp(τ)
coefficient(s, 3), s(0.2), max_order(s)
```

Arithmetic between series requires equal orders. Scalar arithmetic preserves
the order. Differentiation reduces the order by one; integration increases it
by one and supplies a zero constant. `TimeSeries` supports arithmetic, integer
and fractional powers, roots, `exp`, `log`, `sin`, `cos`, `tan`, `sinh`, `cosh`
and `tanh`, where the expansion is analytic. It is a `Number`, not a `Real`:
comparison-based branches in an ODE's right-hand side are not supported.

## Expand an ODE solution

[`taylor_expand`](@ref) constructs a local expansion for an out-of-place right
hand side `f(u, parameters, t)`. Evaluate each component at an offset to advance
the state. For example, `u′ = u`, `u(2) = 1` has solution `exp(t-2)`:

```@example time
expansion = taylor_expand((u, p, t) -> u, [1.0], 2.0; order=16)
[component(0.1) for component in expansion] # state at t = 2.1
```

The recurrence is `u[k+1] = coefficient(f(u,p,t),k)/(k+1)`. The right-hand side
is evaluated once at the expansion point to determine coefficient types, then
with recording scalars to build its arithmetic graph. Each intermediate time
coefficient is computed once and reused at later degrees. The graph is rebuilt
at each new expansion point, so changed parameters are respected. No macro or
source rewriting is required. `taylor_expand` constructs a local series; it
does not select step sizes or guarantee accuracy away from the expansion point.

## SciML solvers

Load an OrdinaryDiffEq solver package to enable [`TaylorMethod`](@ref). It uses
SciML's adaptive controller, solution storage and callback machinery:

```@example time
using OrdinaryDiffEqVerner, SciMLBase
problem = ODEProblem((u, p, t) -> [u[2], -u[1]], [0.0, 1.0], (0.0, 2π))
solution = solve(problem, TaylorMethod(18); abstol=1e-12, reltol=1e-12)
solution(0.5), solution(0.5, Val{1})
```

Both in-place `f!(du,u,p,t)` and out-of-place `f(u,p,t)` vector right-hand
sides are accepted. Options include `saveat`, `dtmax`, `tstops`, and callbacks;
`init` and `step!` provide the usual integrator interface. Dense output evaluates
the local time polynomial, including its derivatives. For fixed steps use
`adaptive=false, dt=...`.

The state dimension stays fixed during a solve. When parameters are uncertain
polynomials, construct polynomial initial states as well (for example,
`TaylorPolynomial.(u0)`) so SciML's solution type can retain the uncertainty.

The last two retained time terms estimate local error; scalar or component-wise
`abstol` and `reltol` scale it. For uncertain states, each error component uses
the sum of absolute polynomial coefficients before applying `internalnorm`.
This includes uncertain coefficients in step selection. The estimate is not
a rigorous error bound, and the method is intended for analytic, nonstiff ODEs.
The [Taylor integration tutorial](../generated/taylor_integration.md) shows
forward/backward Kepler propagation, independent validation with Vern9, and a
polynomial uncertainty map.

Further tutorials cover [dense output, derivatives and continuous events](../generated/taylor_dense_output.md),
[200-revolution CR3BP propagation and Jacobi drift](../generated/taylor_cr3bp.md),
and [CR3BP flow maps and uncertainty moments](../generated/taylor_cr3bp_maps.md).
The latter validates polynomial predictions against independently integrated
trajectories and distinguishes uncertainty truncation from integration error.

## Independent time and uncertainty orders

Use multivariate polynomials as time coefficients to propagate uncertainties:

```@example time
δ, = variables((:δ,); order=2)
expansion = taylor_expand((u, p, t) -> u, [(1 + δ)^2], 0.0; order=20)
coefficient(expansion[1], 20), expansion[1](0.1)
```

Here `order=20` controls time, while `variables(...; order=2)` controls the
total degree in the uncertainty variables. The mixed term `τ²⁰ δ²` is retained.
For `m` uncertainty variables, storage per state is at most
`(time_order+1) * binomial(m+uncertainty_order, m)` coefficients. A single total
degree across time and uncertainties would either discard these mixed terms
or require a much larger basis.

This separation does not add arbitrary per-variable limits to
`TaylorPolynomial`. All its uncertainty variables still share a total-degree
limit. Time series neither change nor reinitialize that algebra, and their
polynomial coefficients follow its usual lifetime rules.

The state and right-hand-side coefficients determine the storage type through
Julia's promotion rules. Integers acquire a division-compatible type; rational
coefficients remain exact for rational arithmetic. Choose the scalar precision
before expansion, and use compatible parameter and time types.
Constructing or copying a `TimeSeries` copies its coefficients, including
mutable polynomial coefficients. Reuse the resulting local series for dense
evaluation at several offsets inside its convergence region.
