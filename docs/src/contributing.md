# Contributing

## Tests

From a repository checkout:

```sh
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

Tests run on Linux with Julia `1.10` and `1` (the latest stable release).
Examples and benchmarks have separate dependency environments.

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
