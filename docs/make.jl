using Documenter, DocumenterCitations, DocumenterCodeBlocks, Literate
using DifferentialAlgebra

include("literate.jl")

examples = [
    "Elementary functions" => "ex01_basic_da",
    "Trigonometric identity" => "ex02_trig_identity",
    "Rational expansion" => "ex03_rational_expansion",
    "Truncation" => "ex04_trig_polynomial",
    "Calculus" => "ex05_diff_integral",
    "Gaussian integral" => "ex06_gaussian_integral",
    "Sombrero at the origin" => "ex07_sombrero_origin",
    "Sombrero gradient" => "ex08_sombrero_gradient",
    "Inverse sine" => "ex09_sinx_inversion",
    "Inverse maps" => "ex10_direct_inverse_map",
    "Matrices and vectors" => "ex11_matrix_vector",
    "Linear algebra" => "ex12_linearalgebra",
    "Legendre basis" => "ex13_legendre_basis",
    "Gaussian ADS" => "ex14_ads_gaussian",
    "Sombrero ADS" => "ex15_ads_sombrero",
    "Implicit Kepler map" => "ex16_kepler_eq_coupled",
    "Pseudo-arclength continuation" => "ex17_kepler_eq_pseudo_arc_length",
    "Newton and fixed point" => "ex18_kepler_eq_fixedpoint",
    "Adaptive Runge–Kutta" => "ex19_kepler_flow_da",
    "OrdinaryDiffEq flow" => "ex20_kepler_flow_ordinarydiffeq",
    "Integrated Kepler ADS" => "ex21_kepler_flow_ads",
    "Picard time expansion" => "ex22_kepler_flow_picard",
    "Serialization" => "ex23_serialization",
    "State transition matrix" => "ode_integration",
    "Kepler Monte Carlo" => "damc_kepler",
    "Analytic Kepler ADS" => "ads_kepler",
]

# Every standalone script must be executed and linked in the documentation.
example_directory = joinpath(@__DIR__, "..", "examples")
registered = Set(name * ".jl" for (_, name) in examples)
available = Set(filter(name -> endswith(name, ".jl"), readdir(example_directory)))
registered == available || error("Example registration does not match examples/")
destination = joinpath(@__DIR__, "src", "generated")

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
        canonical = "https://jackyarndley.github.io/DifferentialAlgebra.jl/",
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Getting started" => "manual/getting-started.md",
            "Polynomial maps" => "manual/maps.md",
            "Automatic domain splitting" => "manual/domain-splitting.md",
            "Coefficient types" => "manual/coefficient-types.md",
            "Configuration and storage" => "manual/configuration.md",
            "Moving from DACEjl" => "manual/dacejl.md",
        ],
        "Examples" => [
            "Overview" => "examples.md",
            [title => "generated/$name.md" for (title, name) in examples]...,
        ],
        "API reference" => "api.md",
        "Mathematical background" => "background.md",
        "References" => "references.md",
        "Contributing" => "contributing.md",
    ],
)
