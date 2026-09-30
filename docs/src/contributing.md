# Contributing

## Tests

From a repository checkout:

```sh
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Tests run on Linux with Julia `1.10` and `1` (the latest stable release).
Examples and benchmarks have separate dependency environments.

## Benchmarks

Run `julia benchmark/setup.jl`, then
`julia --project=benchmark benchmark/orbits.jl orbits.toml` from the repository root.
This checks and times adaptive orbit propagation, first-order state transition
matrices, and higher-order Kepler maps. See the
[benchmark report](https://github.com/jackyarndley/DifferentialAlgebra.jl/blob/main/benchmark/orbit-results.md)
for solver settings, accuracy checks and interpretation of the results.

## Code style and organization

Use [Runic.jl](https://github.com/fredrikekre/Runic.jl) for Julia source formatting.
Install it in a separate environment, then format the source and executable examples:

```sh
julia --project=@runic -e 'using Pkg; Pkg.add(PackageSpec(name="Runic", version=v"1.11.1"))'
julia --project=@runic -e 'using Runic; exit(Runic.main(ARGS))' -- --inplace --docstrings src test ext examples benchmark docs/make.jl
```

Replace `--inplace` with `--check` to verify formatting without modifying files.
CI runs this check on the latest stable Julia release.
Keep the Runic version in these instructions and the CI workflow synchronized;
formatter upgrades may change the output.

The source is organized by responsibility:

- `basis.jl` defines the monomial basis; `polynomial.jl` defines polynomial storage and algebra configuration.
- `arithmetic.jl`, `functions.jl` and `special_scalars.jl` implement numerical kernels.
- `coefficients.jl` and `substitution.jl` implement coefficient access, calculus and substitutions.
- `evaluation.jl`, `linear_algebra.jl` and `statistics.jl` implement maps and derived operations.
- `display.jl` handles polynomial output. Public docstrings live immediately beside the definitions they document.

Use standard Julia arrays and broadcasting. Add methods for `AbstractArray`
interfaces when appropriate, and test views as well as dense arrays.

## Build the documentation

From the repository root:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Open `docs/build/index.html` after a local build. Documenter.jl builds the site,
Literate.jl executes the example scripts, DocumenterCitations.jl renders references,
and DocumenterCodeBlocks.jl adds linked and highlighted code blocks.

Keep examples executable and add references to `docs/src/references.bib`.
Generated pages, build output and local manifests are excluded from version control.

## Publish documentation

The Documentation workflow checks pushes and pull requests. Successful builds
on `main` deploy to [GitHub Pages](https://jackyarndley.github.io/DifferentialAlgebra.jl/).
Pull requests build the site without deploying it.

In the repository's **Settings → Pages**, the source must be **GitHub Actions**.
The workflow uses the repository token and the `github-pages` environment;
no separate deployment key is required.
