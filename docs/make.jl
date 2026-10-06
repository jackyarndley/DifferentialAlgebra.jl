using Documenter, DocumenterCitations, DocumenterCodeBlocks, Literate
using DifferentialAlgebra

include("literate.jl")

# The overview is the single ordered catalog for both navigation and execution.
overview = read(joinpath(@__DIR__, "src", "examples.md"), String)
groups = example_groups(overview)
examples = reduce(vcat, last.(groups))

# Every standalone script must be executed and linked in the documentation.
example_directory = joinpath(@__DIR__, "..", "examples")
registered = Set(name * ".jl" for (_, name) in examples)
available = Set(filter(name -> endswith(name, ".jl"), readdir(example_directory)))
registered == available || error("Example registration does not match examples/")
length(examples) == length(registered) || error("Each example must be registered exactly once")
destination = joinpath(@__DIR__, "src", "generated")
# This directory contains only generated Literate pages and display assets.
isdir(destination) && rm(destination; recursive = true)

for (_, name) in examples
    Literate.markdown(
        joinpath(@__DIR__, "..", "examples", name * ".jl"),
        destination;
        flavor = Literate.DocumenterFlavor(), execute = true,
        postprocess = markdown -> example_images(markdown, destination, name)
    )
end

DocMeta.setdocmeta!(
    DifferentialAlgebra, :DocTestSetup,
    :(using DifferentialAlgebra); recursive = true
)

makedocs(
    root = @__DIR__,
    sitename = "DifferentialAlgebra.jl",
    modules = [DifferentialAlgebra],
    checkdocs = :exports,
    plugins = [
        CitationBibliography(joinpath(@__DIR__, "src", "references.bib"); style = :authoryear),
        CodeBlocks(),
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        # Keep the same tutorial topics and links visible on every page.
        collapselevel = 3,
        canonical = "https://jackyarndley.github.io/DifferentialAlgebra.jl/",
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Getting started" => "manual/getting-started.md",
            "Polynomial maps" => "manual/maps.md",
            "Time expansions" => "manual/time-expansions.md",
            "Automatic domain splitting" => "manual/domain-splitting.md",
            "Continuous ADS surrogates" => "manual/continuous-ads.md",
            "Coefficient types" => "manual/coefficient-types.md",
            "Intervals and Taylor models" => "manual/interval-models.md",
            "Configuration and storage" => "manual/configuration.md",
        ],
        "Tutorials" => [
            "Overview" => "examples.md",
            [topic => [title => "generated/$name.md" for (title, name) in entries] for (topic, entries) in groups]...,
        ],
        "API reference" => "api.md",
        "Mathematical background" => "background.md",
        "References" => "references.md",
        "Contributing" => "contributing.md",
    ],
)

include(joinpath(@__DIR__, "test_navigation.jl"))
