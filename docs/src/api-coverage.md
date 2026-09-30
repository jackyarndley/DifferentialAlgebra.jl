# API coverage and DACE compatibility

The engine was checked against the previous
[DACE.jl wrapper](https://github.com/UoA-AstroGroup/DACE.jl/tree/d390ea1), the
[C++/Julia interface](https://github.com/UoA-AstroGroup/dace/tree/julia-interface),
and the upstream C/C++ mathematical interfaces. The retained functionality is
listed below; C++ spelling is mapped to Julia conventions where appropriate.

When migrating Julia code, replace `using DACE` with `using DifferentialAlgebra`
and qualify functions with `DifferentialAlgebra`. The exception type is `DifferentialAlgebra.DAError`;
the scalar polynomial type remains `DA`. DACE-format persistence is unchanged.

| DACE functionality | Julia API | Verification |
|---|---|---|
| Initialization, epsilon, truncation stack | `init`, `isInitialized`, `getMax*`, `setEps`, `getEps`, `getEpsMac`, `setTO`, `getTO`, `pushTO`, `popTO` | Existing validation and API regression tests |
| Constants, variables, monomials, filled and random polynomials | `DA{T}`, `variable`, `identity`, `monomial`, `filled`, `random` | Three coefficient types, coefficient access and independent multiplication oracle |
| Coefficient extraction and filtering | `cons`, `linear`, `getCoefficient`, `setCoefficient!`, `getMonomial(s)`, `getExponents`, `trim`, `filterMonomials`, `nterms` | Sparse and dense coefficients, bounds, copying and stale-context tests |
| Arithmetic, real/imaginary parts and comparisons | Julia operators; `sqr`, `powi`, `powd`, `multiplyMonomials`; `add!`, `scale!`, `weighted_sum!`, `mul!`, `muladd` | All coefficient types, aliases, table and table-free kernels; comparisons intentionally use constants |
| Derivatives and integrals | `deriv`, `integ`/`integrate`, `divide`, `gradient`, `jacobian`, `hessian`, `hess_stack` | Mixed partials, repeated derivatives, variable division and polynomial matrix tests |
| Roots, powers, exp/log, circular and hyperbolic functions | Julia elementary functions; `root`, `isrt`, `icrt`; `atan(y,x)` corresponds to `atan2` | Identities, domains, constant endpoints, typed expansions and C-core coefficient comparisons |
| Error, gamma, psi and integer-order Bessel functions | `erf`, `erfc`, `gamma`, `loggamma`, `PsiFunction(p,n)`, `digamma`, `polygamma`, `besselj/y/i/k`, `besselix/kx` | Float32/64 and 128/256-bit BigFloat; recurrence identities, scaled Wronskians and independent high-precision reference values |
| Coefficient norms, estimates and bounds | `norm`, `maxNorm`, `orderNorm`, `estimNorm(; errors=true)`, `bound`, `convRadius` | Exact geometric-series estimates, residuals and interval bounds |
| Evaluation, compilation and substitution | `evaluate`, `evalScalar`, `compile`, `compiledDA`, `getDim/Ord/Vars/Terms`, `evaluate!`, `plug`, `replaceVariable`, `scaleVariable`, `translateVariable`, `evalMonomials` | Numeric/polynomial evaluation, missing coordinates, mixed precision, aliases and affine substitutions |
| Inverse polynomial maps | `invert` | Full and partial maps, nonsingular linear systems, coefficient residuals |
| Vectors and matrices | `AlgebraicVector`, `AlgebraicMatrix`, Julia arrays, elementwise vector functions, `vnorm`, `normalize`, `extract`, `concat`, Julia `dot`, `cross`, `inv`, `det`, `transpose`, `frobenius`, `eigh` | Previous eigen/partial tests, typed residuals, copy ownership, empty arrays and unchanged global truncation order |
| Raw/central moments and index generation | `getMultiIndices`, `getRawMoments`, `getCentralMoments` | Independent stored moment data and typed factorial factors |
| Text and binary persistence | `toString`, `fromString`, `parse`, `write(io,p)`, `read(io,DA)` | Exact typed text round trips, C blob fixture, C-core interoperability, malformed/truncated records |

`evaluate` intentionally replaces the old polynomial `eval` name. Julia's `inv`,
`^`, `log(base,p)` and `length`/`nterms` replace C++ `minv`, `pow`, `logb` and
container/monomial size spellings. C++ containers can be replaced by ordinary
Julia arrays; elementwise arithmetic on `AlgebraicVector` retains its old behavior.
The old wrapper's commented-out `hess_vec` was never an exposed function.

CxxWrap reference/allocated types, Eigen conversions, raw C allocation routines,
thread initialization, native version/ABI checks, memory dumps and C error codes
are intentionally absent. Julia owns storage and reports errors as exceptions.
`copy(p)` copies storage; `close(p)` invalidates it. This is a mathematical API
replacement, not an implementation of the removed C ABI.

Binary interoperability follows the legacy Float64 blob format; arbitrary
precision is supported by text IO. The parser accepts DACE records and
whitespace-separated COSY-style records, including Fortran `D` exponents; it does
not claim support for every historical COSY fixed-column variant.

Taylor expansions require an analytic center. General repeated symmetric
eigenvalues do not define a unique eigenbasis, and coefficient-norm extrapolation
does not provide a rigorous convergence certificate. These mathematical limits
are separate from API availability.
