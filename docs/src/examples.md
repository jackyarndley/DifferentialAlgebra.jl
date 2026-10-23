# Tutorials

Every tutorial is a Julia script runnable from the checkout. Shared physical
models live in `examples/support/`. The documentation runs each script's assertions
and displays numerical reports and Makie figures. Choose a
problem below; the result column states what the example lets you check.

For a first pass, use `taylor_expansions.jl`, `implicit_maps.jl`,
`ode_integration.jl`, `interval_models.jl` and `ads_optimization.jl`.
The larger orbit examples then combine these ideas.

## Taylor arithmetic and calculus

| Tutorial | Result |
|:--|:--|
| [Taylor expansions and truncation](generated/taylor_expansions.md) | Coefficient identities and approximation error versus order |
| [Rational functions and convergence](generated/rational_functions.md) | Physical displacement coordinates and useful expansion range |
| [Differentiation and integration](generated/calculus.md) | Gaussian probability checked against `erf` |
| [Multivariate derivatives and removable singularities](generated/multivariate_calculus.md) | Regularized origin series and an analytical gradient/tangent check |
| [Inverse maps](generated/inverse_maps.md) | Scalar and multivariate inverses with composition residuals |
| [Polynomial linear algebra](generated/linear_algebra.md) | Local solves compared with numerical solves and conditioning |
| [Orthogonal polynomials](generated/orthogonal_polynomials.md) | A total-degree Legendre basis checked by its Gram matrix |

## Implicit equations

| Tutorial | Result |
|:--|:--|
| [Implicit Kepler maps](generated/implicit_maps.md) | One inverse map, analytical sensitivities and finite root comparisons |
| [Newton and fixed-point iteration](generated/implicit_solvers.md) | Polynomial root lifting at fixed time or fixed period fraction |
| [Pseudo-arclength continuation](generated/continuation.md) | Local solution curves, unit tangent and residual-selected range |

## Time integration

| Tutorial | Result |
|:--|:--|
| [Adaptive Runge–Kutta](generated/adaptive_integration.md) | A transparent solver propagating a Taylor map |
| [Taylor integration in time](generated/taylor_integration.md) | Separate time/uncertainty orders, reverse flow and Vern9 comparison |
| [Taylor dense output and events](generated/taylor_dense_output.md) | States, derivatives and event times checked analytically |

## Orbit propagation and uncertainty

| Tutorial | Result |
|:--|:--|
| [Orbit integration and state transition matrices](generated/ode_integration.md) | Six-state STM, symplectic residual and quadratic prediction improvement |
| [Long-term CR3BP propagation](generated/taylor_cr3bp.md) | 200 revolutions with independent trajectory checks and Jacobi drift |
| [CR3BP flow maps and uncertainty propagation](generated/taylor_cr3bp_maps.md) | Nonlinear maps, moments and independently integrated samples |
| [Monte Carlo with a Taylor map](generated/damc_kepler.md) | 10,000 samples after 30 revolutions at orders 2, 4 and 8 |
| [Six-variable orbital-element uncertainty](generated/six_variable_ephemeris.md) | Polynomial moments checked by quadrature and three enclosure targets |
| [Six-variable injection uncertainty under J2 gravity](generated/j2_uncertainty.md) | Ground tracks, local-frame dispersion and force/energy checks |

The J2 force model follows [Murri2015](@citet). Time propagation is numerical;
static-function interval certificates do not cover integration error.

## ADS methods and geometry

| Tutorial | Result |
|:--|:--|
| [Comparing estimators and split directions](generated/ads_methods.md) | Four error methods, patch counts, samples and separate certified bounds |
| [Oriented polygon ADS and interval fitting bands](generated/polygon_ads.md) | Axis/oblique frames, exact area preservation and whole-cell bands |
| [An analytic Kepler map](generated/ads_kepler.md) | Large-domain fits, smooth derivatives and a separate explicit certified map |
| [Splitting during orbit propagation](generated/ads_flow.md) | Checkpoint/online splitting and separate static energy certificates |
| [Polygon ADS for an inclined orbital uncertainty ribbon](generated/orbit_polygon_ads.md) | Correlated hexagonal prior mapped into a 3D orbital ribbon |
| [B-plane flyby uncertainty and interval safety classification](generated/flyby_uncertainty.md) | Whole-cell altitude classification and an outgoing-velocity cone |

These examples distinguish sampled discrepancies from uniform bounds. Implicit
root solves and time integrators retain their own numerical error contracts.
See [ADS terminology](manual/domain-splitting.md#ADS-terminology) for partitions,
patches and mapped uncertainty sets.

## Continuous ADS and optimization

| Tutorial | Result |
|:--|:--|
| [C0, C1 and C2 continuity across ADS fits](generated/ads_continuity.md) | Values, gradients and Hessians across joins |
| [Optimization with a C2 oriented ADS surrogate](generated/ads_optimization.md) | Damped Newton against an analytical minimizer and separate function bounds |

## Validated intervals and Taylor models

| Tutorial | Result |
|:--|:--|
| [Interval evaluation and Taylor-model enclosures](generated/interval_models.md) | Stored-polynomial range, original-function range and child-fit error reduction |

## Orbit determination

| Tutorial | Result |
|:--|:--|
| [DAIOD: angles-only orbit maps and uncertainty](generated/daiod.md) | Gauss/Lambert initialization, nonlinear moments and a follow-up measurement update |

## Storage

| Tutorial | Result |
|:--|:--|
| [Serializing compiled maps](generated/serialization.md) | Owned scalar/vector snapshots that survive algebra reinitialization |

## Run locally

From the repository root:

```sh
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
julia --project=examples examples/ads_methods.jl
```

For interactive figures, start `julia --project=examples` and
`include("examples/ads_methods.jl")`. Terminal scripts print their reports.
Four plotting tutorials also save PNGs under `results/`; other figures can be
exported with `CairoMakie.save("figure.png", fig)`. The documentation build
exports every returned tutorial figure automatically.
