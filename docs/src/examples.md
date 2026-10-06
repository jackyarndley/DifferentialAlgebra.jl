# Tutorials

Each tutorial is a standalone Julia script. The documentation executes the
scripts with Literate.jl, showing their numerical results and Makie figures.
The same topic groups appear in the documentation sidebar. Start with polynomial
arithmetic, then choose static maps, time integration or validated intervals
according to the calculation you need.

## Taylor arithmetic and calculus

- [Taylor expansions and truncation](generated/taylor_expansions.md)
- [Rational functions and convergence](generated/rational_functions.md)
- [Differentiation and integration](generated/calculus.md)
- [Multivariate derivatives and removable singularities](generated/multivariate_calculus.md)
- [Inverse maps](generated/inverse_maps.md)
- [Polynomial linear algebra](generated/linear_algebra.md)
- [Orthogonal polynomials](generated/orthogonal_polynomials.md)

Values, derivatives, coefficient decay, matrix conditioning and approximation
errors connect the algebra to numerical behavior.

## Implicit equations

- [Implicit Kepler maps](generated/implicit_maps.md)
- [Newton and fixed-point iteration](generated/implicit_solvers.md)
- [Pseudo-arclength continuation](generated/continuation.md)

## Time integration

- [Adaptive Runge–Kutta](generated/adaptive_integration.md)
- [Taylor integration in time](generated/taylor_integration.md)
- [Taylor dense output and events](generated/taylor_dense_output.md)

Time-step control and dense output are numerical integration features. Their
error estimates do not become interval certificates for the integrated flow.

## Orbit propagation and uncertainty

- [Orbit integration and state transition matrices](generated/ode_integration.md)
- [Long-term CR3BP propagation](generated/taylor_cr3bp.md)
- [CR3BP flow maps and uncertainty propagation](generated/taylor_cr3bp_maps.md)
- [Monte Carlo with a Taylor map](generated/damc_kepler.md)
- [Six-variable orbital-element uncertainty](generated/six_variable_ephemeris.md)
- [Six-variable injection uncertainty under J2 gravity](generated/j2_uncertainty.md)

The element-to-state example compares linear and cubic maps, exact formulas for
retained-polynomial moments, and three distinct interval enclosure targets.
The J2 example adds inclined low-Earth orbits, rotating-Earth ground tracks and
radial/along-track/cross-track dispersion from six uncertain injection coordinates.
The Cartesian J2 force model is documented by [Murri2015](@citet).
Its time integration is numerical; the static-map interval certificates do not
cover integrated trajectories.

## ADS methods and geometry

- [Comparing estimators and split directions](generated/ads_methods.md)
- [Oriented polygon ADS and interval fitting bands](generated/polygon_ads.md)
- [An analytic Kepler map](generated/ads_kepler.md)
- [Splitting during orbit propagation](generated/ads_flow.md)
- [Polygon ADS for an inclined orbital uncertainty ribbon](generated/orbit_polygon_ads.md)
- [B-plane flyby uncertainty and interval safety classification](generated/flyby_uncertainty.md)

Compare heuristic coefficient methods with `IntervalBound`, and box partitions
with oriented/polygon geometry. The implicit Kepler and ODE tutorials explain
why their numerical solvers remain outside the validated static-function layer.
Both now compare all four estimators on related explicit orbital calculations.
The inclined-orbit example compares all four estimators with both geometries
and maps a correlated hexagonal prior into an orbital ribbon. The Mars flyby
example classifies whole polygon cells with guaranteed periapsis-altitude bounds
and plots an outgoing-velocity cone. A boundary cell can remain unresolved;
samples do not resolve its interval certificate.

## Continuous ADS and optimization

- [C0, C1 and C2 continuity across ADS fits](generated/ads_continuity.md)
- [Optimization with a C2 oriented ADS surrogate](generated/ads_optimization.md)

Optional blends give continuous values, gradients or Hessians of the surrogate.
Their function-value certificates do not certify derivatives or a minimizer.

## Validated intervals and Taylor models

- [Interval evaluation and Taylor-model enclosures](generated/interval_models.md)

Distinguish stored-polynomial ranges, interval-valued coefficients and original-
function enclosures with absolute remainders. Bands cover entire plotted cells.

## Orbit determination

- [DAIOD: angles-only orbit maps and uncertainty](generated/daiod.md)

## Storage

- [Serializing compiled maps](generated/serialization.md)

## Run locally

From the repository root, install the example environment:

```sh
julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=examples examples/ads_methods.jl
```

For interactive plotting, start Julia with `julia --project=examples` and run:

```julia
include("examples/ads_methods.jl")
```

Each plotting block ends with its figure, so the interactive frontend or
documentation renderer displays it. Terminal scripts print their results
without opening a viewer. `interval_models.jl`, `polygon_ads.jl`,
`ads_continuity.jl` and `ads_optimization.jl` also save PNGs under `results/`.
For other figures, call CairoMakie's `save("orbit.png", ribbon_fig)` with the
figure you want to export; the documentation build exports every tutorial figure.
