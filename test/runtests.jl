using DifferentialAlgebra
using Test

include("utils.jl")

@testset verbose = true "DifferentialAlgebra tests" begin
    @testset verbose = true "Tutorials" begin
        include("tutorial_tests.jl")
    end

    @testset verbose = true "Validation tests" begin
        include("validation_1.jl")
        include("validation_2.jl")
    end

    @testset verbose = true "Operators" begin
        include("comparison_operators.jl")
    end

    @testset verbose = true "Special Functions" begin
        include("special_functions.jl")
    end

    @testset verbose = true "Linear Algebra" begin
        include("linear_algebra.jl")
    end

    @testset verbose = true "Hessians" begin
        include("hessians.jl")
    end

    @testset verbose = true "Statistics" begin
        include("statistics.jl")
    end

    @testset verbose = true "Variables" begin
        include("variables.jl")
    end

    include("polynomial.jl")
    include("coefficient_types.jl")
    include("regressions.jl")
    include("constructors.jl")
    include("public_api.jl")
    include("arrays.jl")
    include("display.jl")
    include("staticarrays.jl")
end
