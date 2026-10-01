# Orbit propagation benchmarks

Run `benchmark/setup.jl`, then `julia --project=benchmark benchmark/orbits.jl output.toml`
from the repository root. The full suite measures 80 cases with seven samples;
`--quick` measures 52 cases with three samples. Every timed backend is checked
for accuracy before timing. The scripts require no plotting dependencies.

## Problems

All models use nondimensional Cartesian coordinates, gravitational parameter
μ = 1, semi-major axis a = 1, and six independent initial-state variables.
Polynomial variables represent perturbations scaled by 0.001. The outer state
is an ordinary `Vector` for every backend; polynomial storage uses each package's defaults.
The ODE right-hand side fills its derivative buffer in place for every backend.

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

Measured on 2026-10-01 with Julia 1.13.1 on Linux/WSL,
an Intel Core i7-12700, one Julia compute thread and one BLAS thread,
pinned to logical CPU 2. These measurements describe warmed execution on this
machine; compilation and other workloads can have different tradeoffs.

The environment used DifferentiationInterface 0.7.21,
ForwardDiff 1.4.6, OrdinaryDiffEqTsit5 2.1.5,
OrdinaryDiffEqVerner 2.4.2 and SciMLBase 3.57.0.
[Raw results](results-orbits.toml) contain all 80 cases, seven samples per case,
accuracy checks, and hashes of the library and benchmark sources.

### Final state and first-order STM

Times are milliseconds. Float64 returns the state only; the other columns return
both the state and STM. ForwardDiff is the fastest STM implementation in these
cases; higher-order polynomial maps provide additional information.

| Orbit / propagator | Float64 state | DI + ForwardDiff | DifferentialAlgebra | TaylorSeries |
|:--|--:|--:|--:|--:|
| Circular / Tsit5 | 0.0421 | 0.1118 | 1.058 | 11.92 |
| Circular / Vern9 | 0.03361 | 0.0945 | 1.053 | 19.38 |
| Eccentric, three periods / Vern9 | 0.07081 | 0.2258 | 3.564 | 109.7 |
| J₂ / Vern9 | 0.04011 | 0.1139 | 1.706 | 28.94 |
| Circular / analytic Kepler | 0.0001888 | 0.0006638 | 0.003602 | 0.02576 |
| Eccentric, three periods / analytic Kepler | 0.0001998 | 0.0007063 | 0.003783 | 0.02603 |

The polynomial implementations allocate storage for intermediate polynomial
values. ForwardDiff's six-partial dual numbers avoid that representation cost.
Allocation counts and solver statistics for every case are available in the raw
results.

### Higher-order maps

Times are milliseconds for complete propagation. The ratio is TaylorSeries time
divided by DifferentialAlgebra time.

| Orbit / propagator | Order | DifferentialAlgebra | TaylorSeries | Ratio |
|:--|--:|--:|--:|--:|
| Circular / Tsit5 | 4 | 27.09 | 164.1 | 6.06× |
| Circular / Tsit5 | 8 | 713.1 | 2277 | 3.19× |
| Circular / Vern9 | 4 | 21.98 | 163.2 | 7.43× |
| Circular / Vern9 | 8 | 623 | 2163 | 3.47× |
| Eccentric, three periods / Vern9 | 4 | 97.69 | 458.6 | 4.69× |
| Eccentric, three periods / Vern9 | 8 | 1680 | 6965 | 4.15× |
| J₂ / Vern9 | 4 | 42.23 | 216.7 | 5.13× |
| J₂ / Vern9 | 8 | 1303 | 3714 | 2.85× |
| Circular / analytic Kepler | 8 | 5.496 | 7.37 | 1.34× |
| Circular / analytic Kepler | 12 | 105.9 | 141.2 | 1.33× |
| Eccentric, three periods / analytic Kepler | 8 | 3.688 | 7.577 | 2.05× |
| Eccentric, three periods / analytic Kepler | 12 | 71.52 | 147.7 | 2.07× |

Across the 34 paired polynomial cases, TaylorSeries time divided by
DifferentialAlgebra time ranges from 1.33× to 30.79×.
These results apply to the recorded problems and environment.

The largest normalized state and STM errors are 3.652e-9 and
5.214e-9. The largest between-package coefficient difference is
8.878e-10; the largest coefficient error against the
analytic or tighter-tolerance reference is 1.941e-8.
Nonzero-perturbation checks also pass, with a largest error of
0.00004756, including first-order maps that omit nonlinear terms.

DifferentialAlgebra supplies SciML with its underlying coefficient type and the
nominal value of a polynomial. This lets Runge–Kutta tableau entries use scalar
coefficients. TaylorSeries is measured unchanged, so the solver-level comparison
includes this integration difference as well as polynomial arithmetic.

## Comparing revisions

To compare arithmetic and analytic Kepler maps with an instantiated older checkout:

```sh
julia --project=benchmark benchmark/compare_kernels.jl /path/to/older/checkout kernels.toml
```

This loads both arithmetic modules in one process, checks their coefficients,
and interleaves eleven timing samples. It does not load the older SciML extension.
Use `orbits.jl` in each checkout to compare complete solver integrations.
