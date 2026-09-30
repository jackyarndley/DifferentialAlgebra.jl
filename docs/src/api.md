# API reference

```@meta
CurrentModule = DifferentialAlgebra
```

## Construction

```@docs
variables
variable
DA
init
```

## Coefficients and calculus

```@docs
coefficient
constant_term
differentiate
integrate
gradient
jacobian
hessian
```

## Maps

```@docs
evaluate
evaluate!
CompiledMap
compile
invert
```

## Containers and additional operations

```@docs
AlgebraicVector
AlgebraicMatrix
Monomial
eigh
norm
getRawMoments
getCentralMoments
getCoefficient
```

## Other mathematical operations

Standard Julia arithmetic and elementary functions act on polynomials, including
powers, roots, exponentials, logarithms, trigonometric and hyperbolic functions.
Error functions, gamma functions and integer-order Bessel functions have
qualified names such as `DifferentialAlgebra.erf(p)`,
`DifferentialAlgebra.gamma(p)` and `DifferentialAlgebra.besselj(n, p)`.
Their availability for a coefficient type depends on the required scalar methods.

The following functions are available with the `DifferentialAlgebra.` prefix:

| Function | Purpose |
|:--|:--|
| `setCoefficient!(p, exponents, value)` | Change one coefficient in place |
| `getMonomials(p)` | Extract nonzero coefficients and their exponents |
| `trim(p, low, high)` | Retain a range of total degrees |
| `plug(p, i, value)` | Partially evaluate one variable |
| `replaceVariable(p, i, j, scale)` | Substitute `scale*x[j]` for `x[i]` |
| `translateVariable(p, i, scale, shift)` | Substitute `scale*x[i] + shift` |
| `orderNorm(p)` | Group coefficient norms by degree |
| `bound(p)` | Bound the polynomial over coordinates in `[-1, 1]` |
| `estimNorm(p)` / `convRadius(p, tolerance)` | Estimate coefficient growth and convergence heuristically |
| `monomial(exponents, value)` | Construct a single monomial |
| `toString(p)` / `fromString(text, T)` | Convert between polynomials and text |

## Compatibility names

The compatibility names `cons`, `deriv`, `getCoefficient`, and `compiledDA`
remain available. Prefer `constant_term`, `differentiate`, `coefficient`, and
`CompiledMap` in new code. Unlike `getCoefficient`, `coefficient` validates that
exactly one exponent per variable is supplied and that the total degree is valid.
