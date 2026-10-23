# DifferentialAlgebra.jl

Compute a Taylor expansion by running ordinary Julia code on polynomial inputs.
The resulting coefficients give local sensitivities and a reusable approximation
for nearby inputs. Applications include implicit equations, orbit uncertainty,
automatic domain splitting (ADS) and numerical integration.

## Start here

Install with Julia 1.12 or later:

```julia
using Pkg
Pkg.add(url = "https://github.com/jackyarndley/DifferentialAlgebra.jl")
```

```@example introduction
using DifferentialAlgebra
δx, δy = variables((:δx, :δy); order=6)
p = sin(1 + δx) * exp(δy)
(value = p([0.1, 0.2]), gradient = constant_term(gradient(p)))
```

The variables represent perturbations about `(1, 0)`. All retained monomials
have total degree at most six. See [Getting started](manual/getting-started.md)
for coefficients, derivatives and algebra lifetime.

## Choose a workflow

| Goal | Read | Worked result |
|:--|:--|:--|
| Build and reuse local maps | [Polynomial maps](manual/maps.md) | [Implicit Kepler solution and sensitivities](generated/implicit_maps.md) |
| Propagate uncertain states | [Time expansions](manual/time-expansions.md) | [Orbit map and state transition matrix](generated/ode_integration.md) |
| Cover a large input region | [Domain splitting](manual/domain-splitting.md) | [Estimator and geometry comparison](generated/ads_methods.md) |
| Enclose the original function | [Intervals and Taylor models](manual/interval-models.md) | [Certified fitting and error bands](generated/interval_models.md) |
| Use smooth local fits in an optimizer | [Continuous ADS](manual/continuous-ads.md) | [C2 objective and analytical minimizer](generated/ads_optimization.md) |

Polynomial truncation, numerical integration error and interval-certified
remainders are distinct. Each manual states the contract of its calculation.

The [tutorial catalog](examples.md) includes six-variable uncertainty, CR3BP
flow maps, polygonal orbital ribbons and angles-only orbit determination.
Use the [API reference](api.md) for signatures and the
[migration guide](migration.md) when updating existing code.

![Polygon uncertainty mapped into an inclined orbital ribbon](generated/orbit_polygon_ads-figure-2.png)
