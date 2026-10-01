using Test
using Literate
include("literate.jl")

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
