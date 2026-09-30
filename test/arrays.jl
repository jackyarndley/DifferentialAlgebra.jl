using LinearAlgebra

@testset "Native arrays and broadcasting" begin
    for T in (Float32, Float64, BigFloat)
        x, y = variables(T, 2; order = 4)
        tolerance = 1000eps(T)
        v = [1 + x, 2 + y]
        A = [1 + x x * y; y 2 + y]
        @test v isa Vector{DA{T}}
        @test A isa Matrix{DA{T}}
        @test all(iszero, v .* v - [1 + 2x + x^2, 4 + 4y + y^2])
        @test maximum(DifferentialAlgebra.norm, sin.(v) .^ 2 + cos.(v) .^ 2 .- 1) < tolerance
        @test maximum(DifferentialAlgebra.norm, v .* inv.(v) .- 1) < tolerance
        @test DifferentialAlgebra.norm(norm(normalize(v)) - 1) < tolerance
        normalized = copy(v)
        @test normalize!(normalized) === normalized
        @test maximum(DifferentialAlgebra.norm, normalized - normalize(v)) < tolerance
        @test constant_term(v) == T[1, 2]
        @test DifferentialAlgebra.norm(norm(normalize(v, 1), 1) - 1) < tolerance
        @test_throws DifferentialAlgebra.DAError normalize([x, y])
        @test constant_term(A) == T[1 0; 0 2]
        @test all(iszero, differentiate(A, 1) - [one(x) y; zero(x) zero(x)])
        @test all(iszero, differentiate(integrate(A, 2), 2) - A)
        @test all(iszero, DifferentialAlgebra.trim(A, 1) - [x x * y; y y])
        @test all(iszero, DifferentialAlgebra.plug(A, 1, 0) - [one(x) zero(x); y 2 + y])
        @test maximum(DifferentialAlgebra.norm, A * inv(A) - I) < tolerance

        # Views use the same array methods and keep their dimensions.
        column = @view A[:, 1]
        @test constant_term(column) == T[1, 0]
        @test all(iszero, differentiate(column, 1) - [one(x), zero(x)])
        @test evaluate(column, T[1, 2]) == T[2, 2]
        @test CompiledMap(column)(T[1, 2]) == T[2, 2]
        @test all(iszero, jacobian(column) - [one(x) zero(x); zero(x) one(x)])
        @test sin.(column) isa Vector{DA{T}}

        # Native array copies share entries; explicit element copies are independent.
        independent = copy.(A)
        DifferentialAlgebra.setCoefficient!(independent[1, 1], [0, 0], 9)
        @test constant_term(A[1, 1]) == 1
        @test constant_term(independent[1, 1]) == 9
        empty = Matrix{DA{T}}(undef, 0, 2)
        @test size(constant_term(empty)) == (0, 2)
        @test size(differentiate(empty, 1)) == (0, 2)
        @test size(normalize(empty)) == (0, 2)
        @test normalize!(empty) === empty
    end
end
