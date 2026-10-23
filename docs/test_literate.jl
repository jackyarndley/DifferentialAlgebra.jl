using Test
using Literate
using DifferentialAlgebra
import IntervalArithmetic as IA
include("literate.jl")
include("../examples/support/kepler.jl")
include("../examples/support/moments.jl")

@testset "Shared tutorial numerics" begin
    for eccentricity in (0.0, 0.3, 0.999999), mean in (-0.6, 1.0e-4)
        anomaly = kepler_root(mean, eccentricity)
        @test anomaly - eccentricity * sin(anomaly) ≈ mean atol = 2.0e-15
    end
    @test_throws DomainError kepler_root(0.5, 1.0)

    # A polar circular orbit fixes both the rotation and the force convention.
    state = element_state(4.0, 0.0, π / 2, π / 2, 0.0, 0.0, 1.0)
    @test state ≈ [0, 4, 0, 0, 0, 0.5] atol = 1.0e-15
    expected = [0, 0, 0.5, 0, -1 / 16, 0]
    @test two_body(state, 1.0, 0.0) ≈ expected atol = 1.0e-15
    derivative = similar(state)
    two_body!(derivative, state, 1.0, 0.0)
    @test derivative ≈ expected atol = 1.0e-15

    # The same orbital helper must retain guarantees for native intervals.
    enclosed = element_state(IA.interval.((4, 0, 0, 0, 0, 0, 1))...)
    @test all(IA.isguaranteed, enclosed)
    @test all(IA.inf(value) <= expected <= IA.sup(value) for (value, expected) in zip(enclosed, (4, 0, 0, 0, 0.5, 0)))

    # Moments must use all stored terms, independently of arithmetic policy.
    x, = variables(1; order = 2)
    p = 1 + 2x + 3x^2
    previous_tolerance = set_coefficient_tolerance!(10)
    try
        with_order(0) do
            for (moment, mean, variance) in ((uniform_moment, 2, 32 / 15), (normal_moment, 4, 22))
                μ, covariance = polynomial_statistics([p], moment)
                @test only(μ) ≈ mean
                @test only(covariance) ≈ variance
            end
        end
    finally
        set_coefficient_tolerance!(previous_tolerance)
    end
end

@testset "Tutorial topic navigation and registration" begin
    overview = read(joinpath(@__DIR__, "src/examples.md"), String)
    groups = example_groups(overview)
    names = [name for (_, entries) in groups for (_, name) in entries]
    @test !isempty(groups) && all(group -> !isempty(last(group)), groups)
    @test length(names) == length(unique(names))
    @test Set(name * ".jl" for name in names) == Set(filter(n -> endswith(n, ".jl"), readdir(joinpath(@__DIR__, "../examples"))))
    @test_throws ErrorException example_groups("[Uncategorized](generated/missing.md)")
end

@testset "Literate numerical reports and rich figures" begin
    mktempdir() do directory
        source = joinpath(directory, "rendering.jl")
        write(
            source, """
            # # Rendering example
            using CairoMakie
            println("First report: ", 42)
            #-
            fig = Figure(size = (300, 200))
            ax = Axis(fig[1, 1])
            lines!(ax, [0, 1], [0, 1])
            fig
            # A second report and figure must also survive.
            println("Second report: ", 17)
            #-
            fig = Figure(size = (300, 200))
            ax = Axis(fig[1, 1])
            scatter!(ax, [0, 1], [1, 0])
            fig
            """
        )
        output = Literate.markdown(
            source, directory; execute = true,
            flavor = Literate.DocumenterFlavor(),
            postprocess = text -> example_images(text, directory, "rendering")
        )
        markdown = read(output, String)
        @test occursin("First report: 42", markdown)
        @test occursin("Second report: 17", markdown)
        @test !occursin("data:image", markdown)
        images = collect(eachmatch(r"!\[\]\(([^)]+\.png)\)", markdown))
        @test length(images) == 2
        for match in images
            bytes = read(joinpath(directory, only(match.captures)))
            @test bytes[1:8] == UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]
        end
    end
end
