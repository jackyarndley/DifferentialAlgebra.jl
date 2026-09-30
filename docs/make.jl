using Documenter, DocumenterCitations, DocumenterCodeBlocks, Literate
using DifferentialAlgebra

examples = [
    "Elementary functions" => "sine",
    "Differentiation" => "gradient",
    "Map inversion" => "polynomial_inversion",
    "Orbit integration" => "ode_integration",
]

for (_, name) in examples
    Literate.markdown(joinpath(@__DIR__, "..", "examples", name * ".jl"),
                      joinpath(@__DIR__, "src", "generated");
                      flavor=Literate.DocumenterFlavor(), execute=true)
end

DocMeta.setdocmeta!(DifferentialAlgebra, :DocTestSetup,
                   :(using DifferentialAlgebra); recursive=true)

makedocs(
    root=@__DIR__,
    sitename="DifferentialAlgebra.jl",
    modules=[DifferentialAlgebra],
    checkdocs=:exports,
    plugins=[
        CitationBibliography(joinpath(@__DIR__, "src", "references.bib"); style=:authoryear),
        CodeBlocks(),
    ],
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jackyarndley.github.io/DifferentialAlgebra.jl/",
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Getting started" => "manual/getting-started.md",
            "Polynomial maps" => "manual/maps.md",
            "Coefficient types" => "manual/coefficient-types.md",
            "Configuration and storage" => "manual/configuration.md",
        ],
        "Examples" => [title => "generated/$name.md" for (title, name) in examples],
        "API reference" => "api.md",
        "Mathematical background" => "background.md",
        "References" => "references.md",
        "Contributing" => "contributing.md",
    ],
)
