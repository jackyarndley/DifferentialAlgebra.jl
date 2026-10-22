using Test, DifferentialAlgebra, ForwardDiff

include("support/polygon_queries.jl")
using .PolygonQueryTestSupport: exact_evaluate

@testset "Promoted polygon queries and coefficient derivatives" begin
    polygon = ConvexPolygon([(-4, -4), (4, -4), (4, 4), (-4, 4)])
    f(v) = 2 + v[1]^2 + v[1] * v[2] / 4 - v[2]
    fit = adaptive_map(
        f, polygon; directions = [1 1; -1 1], order = 2,
    )
    for point in ([0.25, -0.5], [0.25, 0.25], [4.0, 4.0])
        expected = [2point[1] + point[2] / 4, point[1] / 4 - 1]
        @test ForwardDiff.gradient(fit, point) ≈ expected
        @test ForwardDiff.hessian(fit, point) ≈ [2 1 / 4; 1 / 4 0]
    end
    setprecision(256) do
        point = [BigFloat(1) / 3, BigFloat(1) / 7]
        @test fit(point) isa BigFloat
        @test fit(Tuple(point)) isa BigFloat
        @test fit(point) ≈ f(point) rtol = 32eps(BigFloat)
    end
    @test fit((π, 0)) ≈ f((π, 0))
    @test fit([π, 0]) ≈ f([π, 0])
    vector_fit = adaptive_map(
        v -> [f(v), 3v[1] - v[2]], polygon;
        directions = :axes, order = 2
    )
    @test ForwardDiff.jacobian(vector_fit, [0.25, -0.5]) ≈ [3 / 8 -15 / 16; 3 -1]
    @test vector_fit(BigFloat[0.25, -0.5]) isa Vector{BigFloat}
end

@testset "Filtered routing preserves promoted polynomial and dual evaluation" begin
    DA = DifferentialAlgebra
    f(v) = exp(v[1] + 2v[2])
    polygon = ConvexPolygon([(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)])
    fit = adaptive_map(
        f, polygon; directions = [1 1 // 3; -1 // 3 1],
        order = 3, atol = 1 // 16
    )
    @test length(fit.patches) > 1
    for point in ([0.123, 0.234], [-0.321, 0.432])
        oracle(x) = exact_evaluate(fit, x)
        @test fit(point) == oracle(point)
        @test ForwardDiff.gradient(fit, point) ≈ ForwardDiff.gradient(oracle, point)
        @test ForwardDiff.hessian(fit, point) ≈ ForwardDiff.hessian(oracle, point)
    end
    setprecision(512) do
        point = BigFloat[1 // 7, 2 // 9]
        oracle(x) = exact_evaluate(fit, x)
        @test fit(point) isa BigFloat
        @test fit(point) == oracle(point)
        @test ForwardDiff.gradient(fit, point) ≈ ForwardDiff.gradient(oracle, point)
        @test ForwardDiff.hessian(fit, point) ≈ ForwardDiff.hessian(oracle, point)
    end
end
