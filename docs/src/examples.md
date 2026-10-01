# Examples

These 26 executable tutorials cover Taylor arithmetic, map inversion, orthogonal
polynomials, implicit equations, numerical integration and automatic domain splitting.
Each page is generated from a standalone Julia script with Literate.jl. Numerical
reports and CairoMakie figures appear beside the code that produces them.

## Run locally

From the repository root, install the example environment once:

```sh
julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```

Run a script:

```sh
julia --project=examples examples/ex16_kepler_eq_coupled.jl
```

For an interactive session, start Julia with `julia --project=examples`, then:

```julia
include("examples/ex16_kepler_eq_coupled.jl")
```

Plots end with `fig`, allowing the REPL, notebook or documentation renderer to
display the result. Running a script in a noninteractive terminal executes its
checks and prints its reports; it does not open a plot viewer or write image files.

## Tutorials

| Example | Script |
|:--|:--|
| Elementary functions | [`ex01_basic_da.jl`](generated/ex01_basic_da.md) |
| Trigonometric identity | [`ex02_trig_identity.jl`](generated/ex02_trig_identity.md) |
| Rational expansion | [`ex03_rational_expansion.jl`](generated/ex03_rational_expansion.md) |
| Truncation | [`ex04_trig_polynomial.jl`](generated/ex04_trig_polynomial.md) |
| Calculus | [`ex05_diff_integral.jl`](generated/ex05_diff_integral.md) |
| Gaussian integral | [`ex06_gaussian_integral.jl`](generated/ex06_gaussian_integral.md) |
| Sombrero at the origin | [`ex07_sombrero_origin.jl`](generated/ex07_sombrero_origin.md) |
| Sombrero gradient | [`ex08_sombrero_gradient.jl`](generated/ex08_sombrero_gradient.md) |
| Inverse sine | [`ex09_sinx_inversion.jl`](generated/ex09_sinx_inversion.md) |
| Inverse maps | [`ex10_direct_inverse_map.jl`](generated/ex10_direct_inverse_map.md) |
| Matrices and vectors | [`ex11_matrix_vector.jl`](generated/ex11_matrix_vector.md) |
| Linear algebra | [`ex12_linearalgebra.jl`](generated/ex12_linearalgebra.md) |
| Legendre basis | [`ex13_legendre_basis.jl`](generated/ex13_legendre_basis.md) |
| Gaussian ADS | [`ex14_ads_gaussian.jl`](generated/ex14_ads_gaussian.md) |
| Sombrero ADS | [`ex15_ads_sombrero.jl`](generated/ex15_ads_sombrero.md) |
| Implicit Kepler map | [`ex16_kepler_eq_coupled.jl`](generated/ex16_kepler_eq_coupled.md) |
| Pseudo-arclength continuation | [`ex17_kepler_eq_pseudo_arc_length.jl`](generated/ex17_kepler_eq_pseudo_arc_length.md) |
| Newton and fixed point | [`ex18_kepler_eq_fixedpoint.jl`](generated/ex18_kepler_eq_fixedpoint.md) |
| Adaptive Runge–Kutta | [`ex19_kepler_flow_da.jl`](generated/ex19_kepler_flow_da.md) |
| OrdinaryDiffEq flow | [`ex20_kepler_flow_ordinarydiffeq.jl`](generated/ex20_kepler_flow_ordinarydiffeq.md) |
| Integrated Kepler ADS | [`ex21_kepler_flow_ads.jl`](generated/ex21_kepler_flow_ads.md) |
| Picard time expansion | [`ex22_kepler_flow_picard.jl`](generated/ex22_kepler_flow_picard.md) |
| Serialization | [`ex23_serialization.jl`](generated/ex23_serialization.md) |
| State transition matrix | [`ode_integration.jl`](generated/ode_integration.md) |
| Kepler Monte Carlo | [`damc_kepler.jl`](generated/damc_kepler.md) |
| Analytic Kepler ADS | [`ads_kepler.jl`](generated/ads_kepler.md) |

## Sources and adaptations

The numbered tutorials adapt all 23 numbered examples from
[DACEjl](https://github.com/arma1978/DACEjl/tree/c5d062d277a28b02c89e587e4eed098fd3331036/examples).
They use DifferentialAlgebra's API, ordinary Julia arrays, independent numerical
checks and Makie figures. Their Apache-2.0 license and attribution notices are
retained in the examples directory. Each source file links to its original.

The origin sombrero expansion retains every term through the requested order.
The ADS sombrero uses the removable-singularity series without regularization.
Continuation uses a unit tangent and checks residuals across an interval.
Fixed-point iteration checks all coefficients. The adaptive Runge–Kutta tutorial
implements Dormand–Prince 5(4), which is the tableau used by DACEjl's
`RK78Integrator` at the linked revision. The OrdinaryDiffEq tutorial compares
Vern7 with Vern9. The persistence tutorial uses Julia Serialization for compiled
evaluation maps; DACE text files are not supported.

The initial-coordinate, patch-geometry and pointwise-accuracy checks from the
ancillary ADS validation scripts are incorporated into the ADS tutorials.
Profiling and benchmark scripts are not part of the tutorial collection.

The sine, gradient, inverse-map, state-transition and Monte Carlo tutorials also
cover the five examples in
[DACE.jl](https://github.com/UoA-AstroGroup/DACE.jl/tree/7271c372ce86350e80a95a9c2c0d7e8739d09e83/examples).
The analytic Kepler ADS example complements the integrated flow example.
