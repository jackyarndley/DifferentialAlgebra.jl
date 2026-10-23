# Polynomial maps

```@meta
CurrentModule = DifferentialAlgebra
```

A vector of polynomials is a map. Ordinary Julia arrays, broadcasting and linear
algebra apply; no map-specific container is needed for arithmetic.

## Evaluation and composition

```@example maps
using DifferentialAlgebra
x, y = variables(2; order=5)
f = [x + y^2, y]
(numeric = evaluate(f, [0.1, 0.2]), composed = evaluate(f, [x+y, y]))
```

Numeric coordinates evaluate; polynomial coordinates compose. Missing coordinates
are zero and extra coordinates are ignored, so normally supply exactly one per
variable. Numeric evaluation preserves every stored degree, even under a lower
working order. Polynomial composition follows the current arithmetic policy.

## Compile repeated evaluations

[`CompiledMap`](@ref) owns coefficients and a shared evaluation tree:

```@example maps
map = CompiledMap(f)
map([0.1, 0.2])
```

Reuse output and workspace buffers with [`evaluate!`](@ref):

```@example maps
point = [0.1, 0.2]
result = zeros(noutputs(map))
work = zeros(degree(map) + 1)
evaluate!(result, map, point, work)
```

For this four-argument method, inputs, output and workspace must not overlap.
Scalar types with allocating arithmetic may still allocate. Numeric compiled
evaluation survives algebra
reinitialization; polynomial composition requires the original algebra.

## Invert a map

```@example maps
inverse = invert(f)
residual = evaluate(f, inverse) .- [x, y]
maximum(coefficient_norm, residual)
```

The linear part must be nonsingular. The local inverse includes the original
map's constant shift, so its arguments are output coordinates. A zero-constant
map composes with its inverse to give the identity through the working order.
The residual above uses that order and coefficient tolerance; see
[Compare polynomials](getting-started.md#Compare-polynomials).

Use an augmented square map to solve implicit equations: retain the uncertain
parameters as outputs and invert the equation residual alongside them. The
[Kepler tutorial](../generated/implicit_maps.md) checks this against scalar roots
and analytical sensitivities; [DAIOD](../generated/daiod.md) applies it to noisy
angles and an implicitly determined orbit.

## Linear algebra and calculus

```@example maps
using LinearAlgebra
A = [2+x y; y 4-x]
(determinant = det(A), solution = A \ f, jacobian = jacobian(f))
```

`A*f` is matrix multiplication; `A.*A` and `sin.(f)` are elementwise.
Calculus methods accept arrays and views, preserving shape. `normalize(f)` and
`normalize!(f)` require a nonzero constant norm.

`DifferentialAlgebra.eigenpairs(A)` computes Taylor eigenpairs of symmetric
matrices with distinct constant eigenvalues. Constant and diagonal matrices are
also supported; other repeated eigenvalues generally lack a unique Taylor
basis. The [linear algebra tutorial](../generated/linear_algebra.md) compares
polynomial solves with numerical solves and their conditioning.
