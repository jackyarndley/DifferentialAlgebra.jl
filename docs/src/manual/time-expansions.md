# Time expansions

```@meta
CurrentModule = DifferentialAlgebra
```

Use a univariate time series to advance an ODE, and multivariate coefficients
to retain its dependence on uncertain initial conditions. These two expansion
orders are independent.

## A local time series

[`TimeSeries`](@ref) stores ordinary powers of `τ = t-t₀`; coefficient `k` is
the kth time derivative divided by `k!`. Numeric coefficients need no algebra.

```@example time
using DifferentialAlgebra
τ = TimeSeries([0.0, 1.0, zeros(6)...])
s = exp(τ)
(coefficient = coefficient(s, 3), value = s(0.2), order = max_order(s))
```

Series arithmetic requires equal orders; scalar arithmetic preserves order.
Differentiation reduces it by one and integration increases it by one, with
zero integration constant. Arithmetic, roots, powers, `exp`, `log`, trigonometric
and hyperbolic functions are supported at regular analytic centers.
`TimeSeries` is a `Number`, not a `Real`; comparison-based RHS branches are unsupported.

## Expand an ODE solution

[`taylor_expand`](@ref) accepts `f(u, parameters, t)` and constructs a local
solution. For `u′=u`, `u(2)=1`, evaluate at offset `0.1` to reach time `2.1`:

```@example time
expansion = taylor_expand((u, p, t) -> u, [1.0], 2.0; order=16)
(approximation = expansion[1](0.1), exact = exp(0.1))
```

The recurrence is `u[k+1] = coefficient(f(u,p,t), k)/(k+1)`. An arithmetic graph
records the RHS and reuses each intermediate coefficient. It is rebuilt at each
expansion point, so parameter changes are respected. The local series itself
neither chooses a step size nor guarantees accuracy away from the center.

## Integrate with SciML

Loading an OrdinaryDiffEq solver package enables [`TaylorMethod`](@ref):

```@example time
using OrdinaryDiffEqVerner, SciMLBase
problem = ODEProblem((u, p, t) -> [u[2], -u[1]], [0.0, 1.0], (0.0, 2π))
solution = solve(problem, TaylorMethod(18); abstol=1e-12, reltol=1e-12)
(state = solution(0.5), derivative = solution(0.5, Val{1}))
```

Use a dense `Vector` for the state, with numeric or polynomial elements and a
fixed dimension. Both `f!(du,u,p,t)` and `f(u,p,t)` RHS forms work. SciML supplies
adaptive control, `saveat`, `dtmax`, `tstops`, callbacks, `init` and `step!`.
Dense output evaluates the local polynomial and its derivatives. Use
`adaptive=false, dt=...` for fixed steps.

The last two time terms estimate local error. Scalar or component tolerances
scale each error; for polynomial states the coefficient L1 norm includes uncertain
terms before `internalnorm` is applied. This is a numerical estimate for analytic,
nonstiff ODEs, not a verified integration bound.

## Independent uncertainty order

```@example time
δ, = variables((:δ,); order=2)
expansion = taylor_expand((u, p, t) -> u, [(1+δ)^2], 0.0; order=20)
coefficient(expansion[1], 20), expansion[1](0.1)
```

`order=20` controls time and `variables(...; order=2)` controls uncertainty.
The term `τ²⁰ δ²` is retained. With `m` uncertainty variables and uncertainty
order `q`, each state needs at most `(time_order+1)*binomial(m+q,m)` coefficients.
All uncertainty variables still share a total-degree limit.

Choose precision before constructing coefficients. Promotion determines storage;
rational arithmetic can remain exact. If parameters are polynomial, make initial
states polynomial too (`TaylorPolynomial.(u0)`) so solution storage retains them.
Time series own their coefficients, including copied polynomial storage, and do
not reinitialize the uncertainty algebra.

## Worked checks

- [Taylor integration](../generated/taylor_integration.md): forward/backward Kepler
  propagation against Vern9 and independent uncertainty-map trajectories.
- [Dense output and events](../generated/taylor_dense_output.md): state, derivative
  and event-time errors against an analytical solution.
- [CR3BP flow maps](../generated/taylor_cr3bp_maps.md): independent sample
  propagation, Jacobi residuals and uncertainty moments.

Time integration error and uncertainty truncation error require separate checks.
