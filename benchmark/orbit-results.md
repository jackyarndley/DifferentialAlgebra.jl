# Orbit propagation benchmarks

Run `benchmark/setup.jl`, then `julia --project=benchmark benchmark/orbits.jl output.toml`
from the repository root. The full suite measures 80 cases with seven samples;
`--quick` measures 52 cases with three samples. Every timed backend is checked
for accuracy before timing. The scripts require no plotting dependencies.

## Problems

All models use nondimensional Cartesian coordinates, gravitational parameter
μ = 1, semi-major axis a = 1, and six independent initial-state variables.
Polynomial variables represent perturbations scaled by 0.001. The outer state
is an `SVector` for every backend; polynomial storage uses each package's defaults.

| Orbit | Duration | Solver | Absolute and relative tolerance |
|:--|:--|:--|--:|
| Circular | One period | Tsit5 | 1e-9 |
| Circular | One period | Vern9 | 1e-11 |
| Eccentric, e = 0.4, inclination = 0.5 rad | Three periods | Vern9 | 1e-11 |
| Same eccentric initial state with J₂ R² = 1.08263e-3 × 0.5² | One nominal period | Vern9 | 1e-11 |

The initial step is 0.01 and subsequent steps are adaptive. All backends use
the RMS norm of the nominal state for error control. This makes their controllers
comparable, but **solver tolerances do not directly control STM or higher-order
coefficient errors**. Accepted/rejected steps and RHS evaluations are recorded.
These statistics come from separate untimed solves, using an equivalent fully
seeded dual state for the AD case.
See SciML's [solver options](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/).

The Kepler cases use an independent elliptic Lagrange f/g propagator over the
same circular and eccentric time spans. Newton iteration first solves the scalar
eccentric-anomaly equation, then lifts the solution to polynomial order. Both
polynomial packages execute the same generic implementation. Adaptive solves
use orders 1, 2, 4, 6 and 8; Kepler maps also include orders 10 and 12.
With six variables, order eight has 3,003 coefficients per component and order
twelve has 18,564.

## Comparisons and validation

- Float64 computes only the final state. Its timing gives the cost of propagation
  without sensitivities.
- At first order, DifferentialAlgebra and TaylorSeries return the final state
  and the 6 × 6 STM, including coefficient extraction in their timings.
- DifferentiationInterface's prepared `value_and_jacobian` with
  `AutoForwardDiff(chunksize=6)` differentiates through the same ODE solver and
  returns the final state and STM in one call. The nominal-only controller
  excludes derivatives from step-size selection. Preparation is outside timing;
  no Jacobian values are cached. See the
  [DifferentiationInterface API](https://juliadiff.org/DifferentiationInterface.jl/DifferentiationInterface/stable/api/).
- Higher-order cases return the complete polynomial flow map, providing information
  beyond a first-order STM.
- An independently integrated 42-component state/variational system checks STMs.
  Its Jacobian is obtained by differentiating the RHS, and its error controller
  includes the STM. Unperturbed cases also use the analytic Kepler solution.
- Every Taylor coefficient is compared between packages, grouped by total degree,
  divided by `0.001^degree`, and normalized by the largest coefficient magnitude
  of that degree (or one). This prevents small perturbation scaling from hiding
  high-order errors. Adaptive maps are also checked against Kepler maps or, for
  J₂, a solve with twenty times tighter tolerances.
- A nonzero perturbed initial state is propagated separately and checked against
  the evaluated polynomial map. Solver return codes must indicate success.

## Timing method

Compilation, polynomial-basis construction and AD preparation are excluded.
ODE problem creation, solver setup and output extraction are included. After
warmup, BenchmarkTools measures batches sized for the fast cases, rotating
backend execution order between samples. Garbage collection is enabled, with
an explicit collection before each sample. Reported bytes are cumulative
allocations for one propagation, not peak resident memory.

The output records every timing sample, allocation counts, per-degree errors,
solver statistics, package versions, CPU, operating system and a hash of library
source files. TaylorSeries is the unmodified registered version 0.22.8.

## Results

Measurements use Julia 1.13.1 on Linux/WSL, an Intel Alder Lake CPU, one Julia
compute thread and one BLAS thread, pinned to logical CPU 2. These measurements
describe warmed execution on this machine; compilation and other workloads can
have different tradeoffs.

Measured on 2026-09-30 with DifferentiationInterface 0.7.21, ForwardDiff 1.4.6,
OrdinaryDiffEqTsit5 2.1.5, OrdinaryDiffEqVerner 2.4.2 and SciMLBase 3.57.0.
[Raw results](results-orbits.toml) contain all 80 cases and seven timing samples
per case. DifferentialAlgebra was faster than TaylorSeries in every measured case.

### Final state and first-order STM

Times are milliseconds. Float64 returns the state only; the other columns return
both the state and STM. ForwardDiff is the fastest of the three STM implementations
in these cases; higher-order polynomial maps provide additional information.

| Orbit / propagator | Float64 state | DI + ForwardDiff | DifferentialAlgebra | TaylorSeries |
|:--|--:|--:|--:|--:|
| Circular / Tsit5 | 0.03776 | 0.09001 | 0.9567 | 10.63 |
| Circular / Vern9 | 0.02787 | 0.07185 | 0.8896 | 18.70 |
| Eccentric, three periods / Vern9 | 0.05792 | 0.1660 | 2.836 | 113.6 |
| J₂ / Vern9 | 0.03427 | 0.08397 | 1.502 | 27.64 |
| Circular / analytic Kepler | 0.00007449 | 0.0002836 | 0.003187 | 0.02466 |
| Eccentric, three periods / analytic Kepler | 0.00007850 | 0.0002832 | 0.003297 | 0.02490 |

For circular Vern9, the Float64, ForwardDiff, DifferentialAlgebra and TaylorSeries
cases allocate approximately 6.2 kB, 19.8 kB, 5.02 MB and 56.2 MB respectively.
The polynomial implementations allocate storage for intermediate polynomial values;
ForwardDiff's fixed-size six-partial dual numbers avoid that representation cost.

All four backends take the same first-order step counts: 129/0 accepted/rejected
for circular Tsit5, 40/0 for circular Vern9, 124/4 for the three-period eccentric
case, and 45/1 for J₂. They perform 775, 640, 2,048 and 736 RHS evaluations respectively.

### Higher-order maps

Times are milliseconds for complete propagation. The final column is TaylorSeries
time divided by DifferentialAlgebra time. Order-one comparisons alone do not
predict the cost of these larger maps.

| Orbit / propagator | Order | DifferentialAlgebra | TaylorSeries | Ratio |
|:--|--:|--:|--:|--:|
| Circular / Tsit5 | 4 | 23.82 | 113.1 | 4.75× |
| Circular / Tsit5 | 8 | 589.6 | 2,120 | 3.60× |
| Circular / Vern9 | 4 | 21.72 | 152.6 | 7.02× |
| Circular / Vern9 | 8 | 534.2 | 2,081 | 3.90× |
| Eccentric, three periods / Vern9 | 4 | 86.02 | 355.1 | 4.13× |
| Eccentric, three periods / Vern9 | 8 | 1,558 | 6,593 | 4.23× |
| J₂ / Vern9 | 4 | 55.05 | 221.0 | 4.01× |
| J₂ / Vern9 | 8 | 1,104 | 3,506 | 3.18× |
| Circular / analytic Kepler | 8 | 4.954 | 6.864 | 1.39× |
| Circular / analytic Kepler | 12 | 90.63 | 131.7 | 1.45× |
| Eccentric, three periods / analytic Kepler | 8 | 3.525 | 7.236 | 2.05× |
| Eccentric, three periods / analytic Kepler | 12 | 63.68 | 137.3 | 2.16× |

The largest normalized state and STM errors are `3.65e-9` and `5.21e-9`.
The largest between-package coefficient difference is `8.88e-10`, for eighth-order
eccentric Vern9. Both maps also agree with the independent Kepler reference within
`1.80e-9` under the per-degree normalization. Across all cases, the largest
coefficient error against the analytic or tighter-tolerance reference is `1.94e-8`.
Nonzero-perturbation checks pass too; the largest error, `4.76e-5`, comes from a
first-order map where nonlinear terms are intentionally absent.

### Changes from the previous implementation

The baseline is commit `6126689`. The
[baseline orbit measurements](results-orbits-before.toml) use three samples;
the current measurements use seven, with the same scenarios and solver settings.
These are separate runs, so small time differences should not be treated as precise
speedups. Allocations provide an additional check on the work removed.

| Orbit / solver | Order | Before, ms | After, ms | Time ratio | Allocation reduction |
|:--|--:|--:|--:|--:|--:|
| Circular / Vern9 | 1 | 8.117 | 0.8896 | 9.12× | 64.9% |
| Eccentric, three periods / Vern9 | 1 | 25.03 | 2.836 | 8.83× | 64.9% |
| J₂ / Vern9 | 1 | 9.524 | 1.502 | 6.34× | 57.7% |
| Circular / Vern9 | 8 | 585.9 | 534.2 | 1.10× | 35.8% |
| Eccentric, three periods / Vern9 | 8 | 2,050 | 1,558 | 1.32× | 35.8% |
| J₂ / Vern9 | 8 | 1,287 | 1,104 | 1.17× | 28.9% |

The main change supplies SciML with the underlying coefficient type as well as
the nominal value of a polynomial. Runge–Kutta tableau constants can then use that
scalar type. Previously Vern9 built polynomial tableau entries through its generic
high-precision coefficient-construction path. This fix removes that work and lets
stage arithmetic use scalar coefficients. TaylorSeries is measured unchanged;
part of the solver-level advantage therefore comes from this integration trait,
not just from polynomial arithmetic.

Weighted sums now allow SIMD, and division and square-root recurrences determine
their coefficient-pair ranges before reducing them. Coefficient types remain generic.
The recurrences are tested against independent Horner composition, including sparse
inputs, reduced working order, machine floats, arbitrary precision and exact rational
division. Broader trigonometric and convolution changes were discarded after measurement.

A separate [interleaved arithmetic comparison](results-kernels.toml) loads both source
revisions in the same Julia process, checks their coefficients and takes eleven
samples per implementation. At order twelve it measured 1.31× for dense division,
1.06× for square roots, and 1.22× for weighted sums. Unchanged multiplication was
1.00×. Whole analytic Kepler maps improved only 1–3% in that run; order-eight circular
Kepler was 11% slower. Short-kernel results also varied substantially, including the
unchanged multiplication control. These measurements support the large adaptive-solver
gains and selected dense-kernel gains, but **not a uniform analytic-Kepler speedup**.

To reproduce the arithmetic comparison with an instantiated older checkout:

```sh
julia --project=benchmark benchmark/compare_kernels.jl /path/to/older/checkout kernels.toml
```

This comparison loads only the older arithmetic module, without its SciML extension.
Use `orbits.jl` in each checkout to compare complete solver integrations.
