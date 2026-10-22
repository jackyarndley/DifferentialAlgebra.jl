# Tests

Prepare the test environment once with normal compilation. This also prepares
dependencies that generate native SIMD code during initialization:

```sh
julia --startup-file=no --project=test --check-bounds=yes -O0 \
  -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

Then run the numerical and integration suite in a fresh process:

```sh
julia --startup-file=no --project=test --check-bounds=yes --compile=min -O0 \
  --threads=4 test/runtests.jl
```

CI uses latest stable Julia on Linux. Reduced compilation avoids spending minutes
optimizing one-use test callbacks. This runs every numerical test, including SciML,
ForwardDiff, intervals, BigFloat, ownership and concurrency. `Pkg.test()` runs the
same suite with normal compilation; it takes longer. Optimized kernel performance
and allocation measurements belong in `benchmark/`.

Prepare dependencies with the same bounds-checking flag: Julia requires an
[exact match when loading package images](https://docs.julialang.org/en/v1/devdocs/pkgimg/).

Aqua is an independent CI quality gate, with all optional integrations loaded:

```sh
julia --startup-file=no --project=test test/quality.jl
```

The documentation workflow executes the standalone Literate tutorials and checks
doctests, catalog registration, links and referenced assets. Documentation builds,
dependency installation/precompilation and Aqua are timed separately from the
numerical test process. No numerical tests are hidden behind a slow-suite flag.

## Assessment and coverage

At `3d755ac`, tests comprised 35 dispatched files, roughly 7,300 lines and 46,125
assertions. Latest-Julia Linux CI took 12m43s, including setup. Expensive files
repeated type × precision × order × backend grids and compiled many almost identical
callbacks. Assertion counts were misleading: 2,456 interval-kernel checks took
under a second locally. Other checks compared two library operations or subtracted
polynomials before comparing, allowing working-order/tolerance filtering to hide
errors.

The replacement is organized by contract:

| Subject | Independent evidence retained or added |
|---|---|
| Algebra | Seeded dictionary convolution, both products/backends, exact affine binomials, checked axes/counts, scalar types and extremes |
| Functions | Closed coefficients/derivatives, final-only filtering, high-degree Bessel fixtures, precision, singularities and NaNs |
| Maps | Exact derivatives/inverses/moments, dense evaluation, aliasing, lifetime and actual weak-owner collection |
| Models | Raw polynomial plus absolute remainder checked without certificates, rational discarded-product oracle, metadata, normalization and unsafe-rounding rejection |
| Adaptive maps | Shared construction, refinement/status, split trees and closed faces, certified box/polygon/continuous queries, endpoint jets and flow checkpoints |
| Time/interop | Exact recorded forcings, eager series, solver rejection/mutation/interpolation/discontinuities, Dual nesting and eigenpair residuals |
| Contexts | Scoped factories, inherited policies, exception restoration, concurrent initialization and a real Distributed worker |

Loading checks run first in the fresh test process, before optional imports;
BigFloat/Rational diagonal eigenpairs are exercised without GenericLinearAlgebra.
SafeTestsets isolate names; algebra-mutating tests remain sequential. Shared RHS
fixtures avoid recompiling equivalent closures. Exact allocation assertions and
large numerical sweeps were removed rather than copied into an extended suite.

The assessment also exposed a missing huge-center `hypot` curvature regression.
The separate fix reuses the existing triangular root kernel with balanced scaling;
both backends are checked against a high-precision derivative formula.

## Reference data

`support/special_references.jl` stores selected coefficients from independent
1024-bit scalar origin series and K integrals/derivative identities. Generate it
offline with `julia --project test/generate_references.jl`. Tests do not spend time
recomputing reference quadrature. Decimal fixtures retain 512 bits, comfortably
beyond the tested precision; BigFloat constant references additionally retain their
mpmath provenance beside the assertions.

Measure fresh-process wall time with dependencies already installed, and state the
compiler flags. Do not equate a warmed benchmark, test-summary sum or CI setup time
with that measurement. At `36b9621`, all 475 numerical checks passed on Julia 1.13.2:
local fresh-process runs took 17.7–32.6 seconds with the flags above. The Linux CI
test step took 9 seconds, after 54 seconds of dependency preparation. Aqua and
documentation are separate gates; the PR records their results.
