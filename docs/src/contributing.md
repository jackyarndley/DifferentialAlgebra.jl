# Contributing

## Run the checks

From the repository root:

```sh
julia --project=test -e 'using Pkg; Pkg.instantiate()'
julia --project=test --check-bounds=yes test/runtests.jl
julia --project=test test/quality.jl
```

CI has one **Linux / Julia `1` (latest stable)** test job. Its steps run the
numerical/integration suite, Aqua and Runic; documentation has its own workflow.
The ForwardDiff integration requires version 1 and runs within the complete
suite. For the faster fresh-process command and dependency preparation,
see the repository's [test guide](https://github.com/jackyarndley/DifferentialAlgebra.jl/blob/main/test/README.md).
Julia 1.12 is the minimum supported version, not an additional CI job.

Prefer a small independent oracle over repeated type/order grids or two paths
through the same implementation. Closed forms, exact rationals and externally
derived coefficients are useful. Interval checks must inspect the raw polynomial
plus remainder, separately from stored range certificates; sampled agreement
alone cannot establish a uniform enclosure. Kernel timings and allocation
measurements belong in `benchmark/`, rather than compiler-sensitive test assertions.

SafeTestsets isolates file namespaces; tests run sequentially because default
algebra initialization invalidates earlier polynomials. Loading tests run before
optional imports. The only generated numerical fixture is documented beside its
regeneration script in `test/`.

## Source responsibilities

| Files | Responsibility |
|:--|:--|
| `basis.jl`, `polynomial.jl` | Monomial indexing, storage, algebra ownership and policies |
| `arithmetic.jl`, `functions.jl` | Polynomial arithmetic and elementary kernels |
| `special_functions.jl`, `special_scalars.jl` | Polynomial special functions and scalar precision support |
| `coefficients.jl`, `calculus.jl`, `substitution.jl` | Coefficient queries, derivatives/integrals, coordinate substitution |
| `evaluation.jl`, `linear_algebra.jl`, `statistics.jl` | Compiled maps, eigenpairs, moment operations |
| `ads_driver.jl`, `ads_estimators.jl` | Shared refinement lifecycle, lookup and error indicators |
| `domain_splitting.jl`, `polygon_domains.jl`, `polygon_ads.jl` | Box fitting and exact polygon geometry/adapters |
| `adaptive_flow.jl`, `continuous_ads.jl` | Propagation monitoring and overlapping smooth fits |
| `taylor_models.jl`, `validated_ads.jl` | Dependency-free model/snapshot storage and interfaces |
| `ext/interval_*.jl` | Interval kernels, model arithmetic and adapters to the shared ADS driver |
| `time_series.jl`, `time_kernels.jl`, `time_recurrence.jl` | Public time series, shared coefficient recurrences, recorded RHS graph |
| `ext/DifferentialAlgebra*Ext.jl` | Optional scalar and SciML integration |
| `display.jl`, `precompile.jl` | Readable polynomial output and the small compile workload |

Keep one owner for each algorithm. Thin dispatch methods express type-specific
contracts; they are not duplicated algorithms. Share a kernel only when its
mathematical invariant and ownership rules agree. Rigorous interval bounds and
ordinary coefficient estimates intentionally have different contracts.

The [function inventory](https://github.com/jackyarndley/DifferentialAlgebra.jl/tree/main/dev)
lists named methods, constructors and generated templates across source,
extensions, tests, tutorials and benchmarks. Regenerate it after source changes
with `julia dev/inventory.jl`; `dev/README.md` defines its scope and limitations.

## Readable code and comments

Follow [SciMLStyle](https://docs.sciml.ai/SciMLStyle/dev/) and the surrounding code.
Prefer ordinary arrays, broadcasting and suitable `AbstractArray` methods; test
views when indexing behavior matters. Keep public docstrings next to definitions.
Describe inputs, result, mutation/ownership and numerical restrictions once.
Comments should explain a recurrence, invariant or non-obvious choice; avoid
repeating an assignment or recording patch history in production code.

The [API reference](api.md) identifies exported and supported qualified operations.
Other helpers, including `coefficient_dot`, `filter_terms` and `convergence_radius`,
are internal even when examples/tests use qualified names.

Format with the version pinned in CI:

```sh
julia --project=@runic -e 'using Pkg; Pkg.add(PackageSpec(name="Runic", version=v"1.11.1"))'
julia --project=@runic -e 'using Runic; exit(Runic.main(ARGS))' -- --inplace --docstrings src ext test examples benchmark docs dev
```

Replace `--inplace` with `--check` to verify. Keep formatter versions and shared
compatibility bounds aligned across the root, test, docs, example and benchmark
environments. IntervalArithmetic stays within the 1.0 patch series because its
rounding/metadata hooks require review before widening compatibility.

## Documentation and examples

```sh
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/test_literate.jl
julia --project=docs docs/make.jl
```

Open `docs/build/index.html`. The build executes every standalone example through
Literate, checks public docstrings/doctests, and verifies navigation and local
assets. Register each script once under a topic in `docs/src/examples.md`; that
catalog also drives the sidebar. Generated pages and local manifests stay untracked.

Start an example with the question, then show the calculation and a meaningful
result check. Keep each script runnable from the checkout; put shared physical
models in `examples/support/` so companion tutorials cannot drift apart. Explain
their equations and assumptions in the relevant tutorial. Use `##` for comments
inside Julia blocks so Literate does not split an incomplete function or loop.
End figure blocks with `fig`; separate printed reports from figures with prose
or `#-`. Label sampled errors, display floors and rigorous bounds accurately.
Keep numerical reference provenance and citations beside the relevant calculation.

The documentation workflow builds PRs and deploys successful `main` builds via
GitHub Actions Pages. Source changes do not require a separate deployment key.
