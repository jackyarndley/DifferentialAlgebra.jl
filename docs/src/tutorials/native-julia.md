# Julia engine and precision

The polynomial engine is implemented in Julia. There are no DACE C bindings,
CxxWrap types, Eigen dependency or native build step. Elementary functions use
Taylor recurrences; special functions use scalar values from `SpecialFunctions.jl`
to start their recurrences.

## Coefficient types

`DA{T}` stores real coefficients of type `T`. `DA()` and `DifferentialAlgebra.identity()` default
to `Float64`. `DA(c)` preserves the floating-point type of `c`; integers become
`Float64`. Choose an explicit type with `DA{T}(c)`, `DA{T}(i,c)`,
`DifferentialAlgebra.variable(i,T)` or `DifferentialAlgebra.identity(T)`.

```julia
using DifferentialAlgebra
DifferentialAlgebra.init(6, 2)
x, y = DifferentialAlgebra.identity(Float32)
p = exp(x + y/3)              # DA{Float32}
DifferentialAlgebra.evaluate(p, Float32[0.1, 0.2])

setprecision(256) do
    x, y = DifferentialAlgebra.identity(BigFloat)
    p = log(1 + x + y*y)
    DifferentialAlgebra.evaluate(p, BigFloat[big"0.1", big"0.2"])
end
```

Operations promote coefficient types using Julia's promotion rules. Use `0.1f0`
for a Float32 literal, and construct high-precision constants from strings or
integers inside `setprecision`. Converting an already rounded Float64 constant
cannot recover its missing digits. `setEps(BigFloat(...))` preserves that
threshold's precision; the default threshold is zero.

Arithmetic, elementary functions, calculus, composition, inversion, moments and
all DACE special-function families support all three types. BigFloat gamma and
polygamma use arbitrary-precision scalar recurrences; modified Bessel I uses a
convergent series or a scaled asymptotic expansion, and Bessel K uses its positive
integral with `QuadGK`. Float32 special functions retain Float32 coefficients.
No BigFloat calculation silently converts its coefficients to Float64.
Load `GenericLinearAlgebra` to compute
eigenpairs of non-diagonal BigFloat matrices. Other real coefficient types may
work when their scalar operations are defined, but are not part of the tested
precision contract.

## Compatibility and ownership

The scalar arithmetic, elementary and special functions, coefficient access,
monomials, norms, bounds, factories and truncation settings remain available.
Compiled maps, map inversion, Jacobians/Hessians, moments, symmetric eigenpairs,
`AlgebraicVector` and `AlgebraicMatrix` are retained. Ordinary Julia arrays work.

Use `DifferentialAlgebra.evaluate` for polynomial evaluation. `evalScalar`, `compile` and
`compiledDA` remain available. Missing coordinates are zero; extra coordinates
are ignored. Comparisons retain constant-part semantics; use `DifferentialAlgebra.norm(p-q)`
to compare all coefficients. `iszero(p)` checks the entire polynomial.

A polynomial owns a Julia coefficient vector. Assignment shares the mutable
object, while `copy(p)` creates independent coefficient storage. Array-wrapper
copies also copy their scalars. `close(p)` optionally releases storage early.
Reinitialization invalidates existing polynomials and compiled-map composition;
numeric compiled maps own their coefficients and still work after `init`.

Configure order, dimensions and epsilon before launching concurrent work.
Arithmetic has no global scratch space or lock; concurrent operations may share
read-only inputs but must own their output buffers. Configuration changes must
run without concurrent polynomial calculations. Both `eigh` and map inversion
use local degree limits and leave the configured truncation order unchanged.

`eigh` supports distinct constant eigenvalues, constant matrices and diagonal
polynomial matrices. Other repeated or numerically indistinguishable eigenvalues
are rejected because a unique Taylor eigenbasis is not generally defined.

## Performance

Coefficients occupy a degree-ordered dense prefix. Constants use one coefficient;
arithmetic reserves only the degree range it can produce. Multiplication uses
packed monomial lookup tables and a symmetric kernel for squares. Elementary
functions use homogeneous coefficient recurrences. Composition tracks active
monomials, and map inversion lifts one degree at a time using a fixed linear solve.
Division directly solves `denominator*result=numerator` coefficient by coefficient;
square roots solve `result^2=input` using each symmetric monomial pair once.
`muladd(scale,p,q)` fuses scalar multiplication and polynomial addition into one
output buffer. `add!`, `scale!`, `weighted_sum!` and `mul!` reuse existing buffers.

`init(order, variables; table_bytes=32*1024^2)` bounds the combined product and
recurrence tables. Larger bases use combinatorial index calculation and Horner
composition. The budget excludes basis metadata and polynomial buffers. Dense
storage can be costly for very large, sparsely populated bases; performance
depends on order, dimension, sparsity and coefficient type.

Reuse buffers to avoid allocations in Float32/Float64 kernels:

```julia
using LinearAlgebra
DifferentialAlgebra.init(5, 2)
x, y = DifferentialAlgebra.identity()
p, q = sin(x)*exp(y), x+y*y
product = zero(p)
mul!(product, p, q)            # Allocates capacity on first use; then reuses it.

map = DifferentialAlgebra.compile([p, q])
args, result = [0.1, 0.2], zeros(2)
work = zeros(DifferentialAlgebra.getOrd(map)+1)
DifferentialAlgebra.evaluate!(result, map, args, work)
```

`mul!` allows aliases, using a temporary when needed. Numeric evaluation requires
disjoint inputs, output and workspace. BigFloat scalar arithmetic still allocates.
Compile once when evaluating a map repeatedly. The DiffEqBase extension loads
only when DiffEqBase is used, so ordinary calculations do not load the ODE stack.

## Substitution and persistence

`plug(p,i,value)` partially evaluates one coordinate. `replaceVariable(p,i,j,a)`
substitutes `a*x[j]` for `x[i]`; `scaleVariable(p,i,a)` and
`translateVariable(p,i,a,c)` substitute `a*x[i]` and `a*x[i]+c`. These operations
preserve the number of independent variables. Like the C core, `plug` respects
the current truncation order, while scaling, translation and replacement retain
terms through the initialized maximum order.

`evalMonomials(p,q)` takes a coefficient dot product. `filterMonomials(p,mask)`
retains monomials present in the mask. `nterms(p)` counts nonzero coefficients;
Julia's `size(p)` still describes a scalar. `estimNorm(p; errors=true)` returns
estimates and positive fit residuals; `convRadius(p,tolerance)` extrapolates the
first omitted order. These estimates are heuristic, not certified error bounds.

`fromString(toString(p), T)` and `parse(DA{T}, text)` round-trip text at the chosen
coefficient precision. Float64 output uses the C core's DACE text layout.
`write(io,p)` and `read(io,DA)` interoperate with packed DACE binary blobs for
Float64. The reader accepts both byte orders, validates records, and discards
terms outside the current initialization. Legacy blobs have 32-bit exponent
indices and Float64 coefficients; use text for larger bases or other precisions.

See [Functionality audit](../api-coverage.md) for the correspondence with the
previous Julia wrapper and upstream interfaces. The reproducible orbit comparison
with TaylorSeries is in `benchmark/integration.jl`, with setup and results in
`benchmark/README.md`.
