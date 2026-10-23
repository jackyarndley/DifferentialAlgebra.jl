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

Use `DifferentialAlgebra.eigenpairs` for polynomial eigenpairs.

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

`TaylorPolynomial{Interval{T}}` uses ordinary truncated polynomial arithmetic
with interval coefficients. `TaylorModel` wraps that polynomial with an absolute
remainder and validity-domain/normalization data. Construct adaptive partitions
with `adaptive_map(...; estimator=IntervalBound())`.

```@docs
enclose
TaylorModel
CompiledTaylorModel
taylor_models
polynomial
remainder
domain
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
continuous_map
ContinuousTaylorMap
blend_weights
```

## Statistical moments

These operations use the `DifferentialAlgebra.` prefix.

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

The following operations also use the `DifferentialAlgebra.` prefix.

```@docs
trim
substitute
translate_variable
nthroot
degree_norms
bounds
estimate_norms
multiindices
hessian_tensor
```
