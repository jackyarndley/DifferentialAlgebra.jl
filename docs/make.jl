using Documenter, Literate, DifferentialAlgebra

examples = [
    "Sine function" => "sine",
    "Gradient" => "gradient",
    "Polynomial inversion" => "polynomial_inversion",
    "ODE integration" => "ode_integration",
]

for (_, name) in examples
    Literate.markdown(joinpath(@__DIR__, "..", "examples", name * ".jl"),
                      joinpath(@__DIR__, "src", "generated");
                      flavor=Literate.DocumenterFlavor(), execute=true)
end

makedocs(
    root=@__DIR__,
    sitename="DifferentialAlgebra.jl",
    format=Documenter.HTML(prettyurls=get(ENV, "CI", "false") == "true"),
    pages=[
        "Home" => "index.md",
        "Engine and precision" => "tutorials/native-julia.md",
        "Development" => "tutorials/setting-up-your-development-environment.md",
        "Examples" => [title => "generated/$name.md" for (title, name) in examples],
        "API" => "api.md",
        "API coverage" => "api-coverage.md",
    ],
)
