# Orbit integration comparison

From the repository root:

```sh
julia --startup-file=no benchmark/setup.jl
julia --startup-file=no --project=benchmark benchmark/integration.jl comparison.toml
```

`--quick` uses three samples instead of seven. Setup develops this checkout and
pins TaylorSeries 0.22.8, the registered release used for this comparison.
The output file records every timing sample, median time, allocated bytes and
coefficient error. Julia, TaylorSeries, OS and CPU versions are recorded too.
Setup is compatible with Julia 1.10 and current Julia.

## Workload and checks

Both packages run exactly the same six-state, six-variable two-body problem:
unit gravitational parameter, initial circular orbit `(1,0,0,0,1,0)`, and six
independent initial-state perturbations scaled by 0.001. Integration spans one
period, `2pi`.

- RK4 uses 256 fixed steps and the same Julia function for both libraries.
- Float64 polynomial orders are 2, 3 and 5; Float32 and 256-bit BigFloat use order 3.
- OrdinaryDiffEq's Vern9 comparison uses order 3 Float64, 64 fixed steps and the
  same RHS, solver settings and constant-part norm for both libraries. Adaptivity
  is disabled, so both perform the same requested integration steps.
- Every stored Taylor coefficient is compared. Errors are grouped by total degree
  and divided by `0.001^degree` before normalization, so small high-order
  perturbation coefficients cannot hide errors.
- The nominal orbit is checked against its exact period return. A nonzero
  perturbation is independently propagated as a scalar orbit and compared with
  the polynomial map. RK4's scalar error must decrease by approximately 16 when
  halving the step size.

Initialization and compilation are excluded from timing. Each implementation is
warmed twice. Execution order alternates between samples; garbage collection is
enabled during propagation, with an explicit collection before each sample.
Results include state/stage allocation and, for Vern9, solver setup. They do not
include package import or polynomial-basis initialization.

TaylorSeries uses `variables!` so OrdinaryDiffEq can construct typed zero values
in the correct default JetSpace. The standalone RK4 also works with an explicit
JetSpace. No library source or internal multiplication implementation is patched.

## Results

Measured on 2026-09-30, Linux/WSL x86-64, Intel Alder Lake, pinned to logical CPU 2,
with one Julia compute thread. TaylorSeries was 0.22.8 and OrdinaryDiffEqVerner was
2.4.2 on both Julia versions. Times below are medians of seven samples for one
complete orbit; the last column gives the corresponding Julia 1.10 speedup.

| Integrator / coefficients | Order | DifferentialAlgebra, Julia 1.13.1 | TaylorSeries, Julia 1.13.1 | Speedup | Speedup, Julia 1.10.12 |
|---|---:|---:|---:|---:|---:|
| RK4 / Float64 | 2 | 1.827 ms | 6.687 ms | 3.66× | 4.12× |
| RK4 / Float64 | 3 | 4.494 ms | 11.365 ms | 2.53× | 3.38× |
| RK4 / Float64 | 5 | 40.799 ms | 88.537 ms | 2.17× | 2.34× |
| RK4 / Float32 | 3 | 4.100 ms | 10.602 ms | 2.59× | 5.70× |
| RK4 / BigFloat, 256 bits | 3 | 700.433 ms | 993.927 ms | 1.42× | 1.22× |
| Vern9 / Float64 | 3 | 21.867 ms | 106.123 ms | 4.85× | 3.14× |

Raw timings and allocations: [Julia 1.10](results-linux-110.toml),
[Julia 1.13](results-linux-113.toml). DifferentialAlgebra allocated fewer bytes in every case.
On Julia 1.13 it used 50–65% fewer bytes for the Float32/Float64 RK4 cases, 42.5%
fewer for BigFloat, and 72.8% fewer for Vern9. BigFloat scalar arithmetic still
allocates substantially; reusable Float32/Float64 multiplication and numeric
evaluation buffers retain their tested zero-allocation behavior.

The largest normalized coefficient differences were `1.35e-14` for Float64,
`3.01e-6` for Float32 and `8.39e-76` for BigFloat. All independent orbit checks
passed. These are warm-runtime comparisons for this workload and machine, not
a guarantee for every polynomial order, sparsity pattern or application.
