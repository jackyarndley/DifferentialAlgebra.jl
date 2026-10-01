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
| Elementary functions | [`ex01_basic_da.jl`](ex01_basic_da.jl) |
| Trigonometric identity | [`ex02_trig_identity.jl`](ex02_trig_identity.jl) |
| Rational expansion | [`ex03_rational_expansion.jl`](ex03_rational_expansion.jl) |
| Truncation | [`ex04_trig_polynomial.jl`](ex04_trig_polynomial.jl) |
| Calculus | [`ex05_diff_integral.jl`](ex05_diff_integral.jl) |
| Gaussian integral | [`ex06_gaussian_integral.jl`](ex06_gaussian_integral.jl) |
| Sombrero at the origin | [`ex07_sombrero_origin.jl`](ex07_sombrero_origin.jl) |
| Sombrero gradient | [`ex08_sombrero_gradient.jl`](ex08_sombrero_gradient.jl) |
| Inverse sine | [`ex09_sinx_inversion.jl`](ex09_sinx_inversion.jl) |
| Inverse maps | [`ex10_direct_inverse_map.jl`](ex10_direct_inverse_map.jl) |
| Matrices and vectors | [`ex11_matrix_vector.jl`](ex11_matrix_vector.jl) |
| Linear algebra | [`ex12_linearalgebra.jl`](ex12_linearalgebra.jl) |
| Legendre basis | [`ex13_legendre_basis.jl`](ex13_legendre_basis.jl) |
| Gaussian ADS | [`ex14_ads_gaussian.jl`](ex14_ads_gaussian.jl) |
| Sombrero ADS | [`ex15_ads_sombrero.jl`](ex15_ads_sombrero.jl) |
| Implicit Kepler map | [`ex16_kepler_eq_coupled.jl`](ex16_kepler_eq_coupled.jl) |
| Pseudo-arclength continuation | [`ex17_kepler_eq_pseudo_arc_length.jl`](ex17_kepler_eq_pseudo_arc_length.jl) |
| Newton and fixed point | [`ex18_kepler_eq_fixedpoint.jl`](ex18_kepler_eq_fixedpoint.jl) |
| Adaptive Runge–Kutta | [`ex19_kepler_flow_da.jl`](ex19_kepler_flow_da.jl) |
| OrdinaryDiffEq flow | [`ex20_kepler_flow_ordinarydiffeq.jl`](ex20_kepler_flow_ordinarydiffeq.jl) |
| Integrated Kepler ADS | [`ex21_kepler_flow_ads.jl`](ex21_kepler_flow_ads.jl) |
| Picard time expansion | [`ex22_kepler_flow_picard.jl`](ex22_kepler_flow_picard.jl) |
| Serialization | [`ex23_serialization.jl`](ex23_serialization.jl) |
| State transition matrix | [`ode_integration.jl`](ode_integration.jl) |
| Kepler Monte Carlo | [`damc_kepler.jl`](damc_kepler.jl) |
| Analytic Kepler ADS | [`ads_kepler.jl`](ads_kepler.jl) |

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
