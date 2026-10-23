# Tests

## Run

Prepare dependencies once with normal compilation, including packages that
generate native SIMD code during initialization:

```sh
julia --startup-file=no --project=test --check-bounds=yes -O0 \
  -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

Run the complete numerical/integration suite in a fresh process:

```sh
julia --startup-file=no --project=test --check-bounds=yes --compile=min -O0 \
  --threads=4 test/runtests.jl
```

`Pkg.test()` runs the same suite with normal compilation. The reduced-compilation
command avoids optimizing one-use test callbacks; no tests are hidden behind a
slow-suite switch. Prepare dependencies with matching bounds-check flags because
Julia requires an [exact package-image match](https://docs.julialang.org/en/v1/devdocs/pkgimg/).

Run Aqua with all optional integrations loaded:

```sh
julia --startup-file=no --project=test test/quality.jl
```

CI has one Linux job using Julia `1` (latest stable), with numerical/integration,
Aqua and Runic checks as steps. ForwardDiff 1 is required and tested within the
complete suite. The documentation workflow executes all tutorials and checks
doctests, registration, local links and figures.

## Coverage by contract

| Subject | Independent evidence |
|:--|:--|
| Algebra | Dictionary convolution for ordinary/aliased products, both multiplication backends, exact affine binomials, axes/counts and scalar extremes |
| Functions | Closed coefficients/derivatives under both backends, final-only filtering, high-degree Bessel fixtures, precision and singularities |
| Maps | Exact repeated calculus, shifted inverses and moments; mixed-owner/workspace rejection, snapshot lifetime and owner collection |
| Models | Independent norm/tail oracles, raw nonlinear remainders, failed-initialization ownership, metadata, normalization and unsafe-rounding rejection |
| ADS | Explicit estimator formulas, exact blend weights, closed-face queries, refinement/status and nonlinear flow restart/final-time contracts |
| Time/interop | Closed-form eager/recorded coefficients, independent uncertainty degrees, local epochs, backward dense output, solver rejection/mutation, nested duals and eigenpair residuals |
| Context | Scoped factories, inherited policies, restoration, concurrent initialization and a Distributed worker |

Loading checks verify all four optional extensions are absent before imports.
SafeTestsets isolates names; algebra-mutating subjects run sequentially. Shared RHS
fixtures avoid recompiling equivalent closures. Add focused independent regressions;
avoid Cartesian grids, filtered polynomial differences that hide errors, and exact
allocation counts.
Optimized performance belongs in `benchmark/`.

## Reference data

`support/special_references.jl` contains selected coefficients from independent
1024-bit scalar origin series and K integrals/derivative identities. Regenerate
with `julia --project test/generate_references.jl`. Fixtures retain 512 bits;
additional BigFloat constants identify their mpmath provenance beside assertions.
Tests use stored data rather than repeating reference quadrature.

Measure fresh-process wall time with dependencies prepared and record compiler
flags. Keep setup, precompilation, Aqua and documentation timings separate.
Historical suite counts and benchmark comparisons belong in the PR/Git history,
not in assertions about the current checkout.
