# Polynomial maps

```@meta
CurrentModule = DifferentialAlgebra
```

A vector of polynomials represents a map. Ordinary Julia vectors are sufficient;
`AlgebraicVector` is available when elementwise vector arithmetic is desired.

## Evaluation and composition

```@example maps
using DifferentialAlgebra
x, y = variables(2; order=5)
f = [x + y^2, y]
evaluate(f, [0.1, 0.2])
```

Numeric coordinates evaluate the map. Polynomial coordinates compose it:

```@example maps
evaluate(f, [x + y, y])
```

Missing coordinates are treated as zero; extra coordinates are ignored. Supply
one coordinate per independent variable when evaluating a complete state.

## Compile repeated evaluations

[`CompiledMap`](@ref) stores a shared evaluation tree for repeated use.

```@example maps
map = CompiledMap(f)
map([0.1, 0.2])
```

For repeated evaluations into existing buffers, use [`evaluate!`](@ref):

```@example maps
point = [0.1, 0.2]
result = zeros(2)
work = zeros(DifferentialAlgebra.getOrd(map) + 1)
evaluate!(result, map, point, work)
```

The inputs, output and workspace must not overlap. Allocation behavior also
depends on the scalar coefficient arithmetic.

## Invert a map

```@example maps
inverse = invert(f)
residual = evaluate(f, inverse) .- [x, y]
maximum(DifferentialAlgebra.norm, residual)
```

The linear part must be nonsingular. The inverse is constructed locally about
the input origin and includes the shift by the map's constant value, so its
arguments are output coordinates. For a zero-constant map, composition gives
the identity through the configured order. See the
[inversion example](../generated/polynomial_inversion.md).

## Linear algebra

`jacobian(f)` gives polynomial partial derivatives. Standard Julia matrix
operations apply to arrays of `DA` values. `DifferentialAlgebra.eigh(A)`
computes Taylor eigenpairs for symmetric matrices with distinct constant
eigenvalues; constant and diagonal matrices are also supported. Other repeated
eigenvalues do not generally define a unique Taylor eigenbasis.
