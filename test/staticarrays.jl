using StaticArrays, LinearAlgebra

@testset "Static array interoperability" begin
    x, y = variables(2; order = 4)
    f = SVector(x + y^2, y)
    @test evaluate(f, SVector(0.1, 0.2)) ≈ [0.14, 0.2]
    @test maximum(coefficient_norm, evaluate(f, invert(f)) - [x, y]) < 1.0e-14
    @test sin.(f) isa SVector
    A = SMatrix{2, 2}(2 + x, y, y, 3 - x)
    @test coefficient_norm(det(A) - ((2 + x) * (3 - x) - y^2)) < 1.0e-14
    @test maximum(coefficient_norm, A * inv(A) - I) < 1.0e-14
    @test coefficient(x * y, SVector(1, 1)) == 1
    @test constant_term(f) == SVector(0.0, 0.0)
    @test constant_term(jacobian(f)) == [1 0;0 1]
    map = CompiledMap(f)
    out, work = MVector{2, Float64}(undef), MVector{3, Float64}(undef)
    @test evaluate!(out, map, SVector(0.1, 0.2), work) === out
    @test out ≈ [0.14, 0.2]
    @test maximum(coefficient_norm, differentiate(f, 2) - SVector(2y, one(y))) == 0
end
