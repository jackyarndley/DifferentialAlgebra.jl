# API reference

```@meta
CurrentModule = DifferentialAlgebra
```

## Construction and algebra

```@docs
TaylorPolynomial
variables
variable
initialize!
coefficient_type
nvariables
max_order
truncation_order
set_truncation_order!
with_order
coefficient_tolerance
set_coefficient_tolerance!
```

## Coefficients and calculus

```@docs
coefficient
set_coefficient!
constant_term
linear_part
Monomial
monomials
monomial
exponents
degree
coefficient_norm
differentiate
integrate
gradient
jacobian
hessian
```

## Maps and linear algebra

```@docs
evaluate
evaluate!
CompiledMap
compile
noutputs
nnodes
invert
eigenpairs
```

## Time expansions

```@docs
TimeSeries
taylor_expand
TaylorMethod
```

## Intervals and Taylor models

```@docs
enclose
TaylorModel
CompiledTaylorModel
taylor_models
polynomial
remainder
domain
validated_adaptive_map
PiecewiseTaylorModel
TaylorModelPatch
```

## Automatic domain splitting

```@docs
adaptive_map
adaptive_flow
GuardedTail
ExtrapolatedTail
LastTerms
IntervalBound
PiecewiseTaylorMap
TaylorPatch
ConvexPolygon
polygon_vertices
domain_area
split_directions
PiecewisePolygonMap
PolygonPatch
```

## Statistical moments

```@docs
raw_moments
central_moments
random_polynomial
```

## Other operations

Julia arithmetic and elementary functions act on polynomials, including powers,
roots, exponentials, logarithms, trigonometric and hyperbolic functions. Load
SpecialFunctions.jl for error functions, gamma functions and integer-order
Bessel functions. Availability depends on the corresponding scalar methods.

These additional operations use the `DifferentialAlgebra.` prefix:

| Function | Purpose |
|:--|:--|
| `trim(p, low, high)` | Retain a range of total degrees |
| `substitute(p, i, value)` | Evaluate one variable at a scalar |
| `replace_variable(p, i, j, scale)` | Substitute `scale*x[j]` for `x[i]` |
| `translate_variable(p, i, scale, shift)` | Substitute `scale*x[i] + shift` |
| `scale_variable(p, i, scale)` | Scale one coordinate |
| `divide_variable(p, i, power)` | Divide by a variable power when every term is divisible |
| `nthroot(p, n)` | Expand an integer-order root |
| `degree_norms(p)` | Group coefficient norms by degree |
| `bounds(p)` | Heuristic monomial bounds on `[-1, 1]`; floating rounding is not enclosed |
| `estimate_norms(p)` / `convergence_radius(p, tolerance)` | Estimate growth and convergence heuristically |
| `coefficient_product(p, q)` / `coefficient_dot(p, q)` | Multiply or sum corresponding coefficients |
| `filter_terms(p, mask)` | Keep the monomials present in a mask |
| `multiindices(order, n)` | Enumerate exponent vectors through a total degree |
| `nterms(p)` | Count nonzero monomials |
| `hessian_tensor(map)` | Stack component Hessians along a third dimension |
