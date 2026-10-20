using Test, DifferentialAlgebra
using LinearAlgebra, GenericLinearAlgebra, ForwardDiff

include("support/processes.jl")
include("support/coefficients.jl")

import ForwardDiff: Dual

@testset "Diagonal eigenpairs without optional scalar solvers" begin
    script = "include(" *
        repr(joinpath(@__DIR__, "support", "diagonal_eigenpairs.jl")) * ")"
    @test fresh_julia(script)
end

@testset "Eigenpair validation and scalar result types" begin
    for T in (Int, Float32, Float64, BigFloat)
        variables(T, 1; order = 2)
        # A one-element problem also supplies the working type for solvers that
        # cannot handle empty systems (including GenericLinearAlgebra BigFloat).
        expected = eigen(Symmetric(ones(T, 1, 1)))
        values, vectors = DifferentialAlgebra.eigenpairs(Matrix{TaylorPolynomial{T}}(undef, 0, 0))
        @test coefficient_type(eltype(values)) === eltype(expected.values)
        @test coefficient_type(eltype(vectors)) === eltype(expected.vectors)
        @test isempty(values) && size(vectors) == (0, 0)
        for A in (T[2 1; 1 4], T[2 0; 0 4])
            expected = eigen(Symmetric(A))
            values, vectors = DifferentialAlgebra.eigenpairs(TaylorPolynomial{T}.(A))
            @test coefficient_type(eltype(values)) === eltype(expected.values)
            @test coefficient_type(eltype(vectors)) === eltype(expected.vectors)
            @test constant_term(values) ≈ expected.values
            @test A * constant_term(vectors) ≈ constant_term(vectors) * Diagonal(expected.values)
        end
    end
    for A in ([1 2; 0 3],)
        @test_throws ArgumentError DifferentialAlgebra.eigenpairs(A)
        @test_throws ArgumentError DifferentialAlgebra.eigenpairs(TaylorPolynomial.(A))
    end
    x, = variables(Int, 1; order = 2)
    A = [2 + x one(x); one(x) 4 - x]
    values, vectors = DifferentialAlgebra.eigenpairs(A)
    @test eltype(values) === eltype(vectors) === TaylorPolynomial{Float64}
    @test maximum(coefficient_norm, A * vectors - vectors * Diagonal(values)) < 1.0e-13
    @test_throws DimensionMismatch DifferentialAlgebra.eigenpairs(ones(2, 3))
    x, = variables(1; order = 4)
    A = [2 + x one(x) + x^4; one(x) 4 + x]
    B = [2 + x one(x) + 1.0e-8x; one(x) 4 + x]
    with_order(2) do
        set_coefficient_tolerance!(1.0e-4)
        values, vectors = DifferentialAlgebra.eigenpairs(A)
        @test constant_term(values) ≈ [3 - sqrt(2), 3 + sqrt(2)]
        @test all(p -> coefficient(p, [1]) ≈ 1, values)
        @test all(p -> iszero(coefficient(p, [2])), values)
        @test maximum(degree, vectors) == 0
        @test_throws ArgumentError DifferentialAlgebra.eigenpairs(B)
        @test truncation_order() == 2
    end
end

@testset verbose = true "Eigendecomposition" begin

    DifferentialAlgebra.initialize!(1, 6)

    A = [4.0 2.0 1.0; 2.0 3.0 1.0; 1.0 1.0 2.0]
    Ada = Matrix{TaylorPolynomial}(undef, 3, 3)
    Afd = Matrix{Dual{Nothing, Float64, 6}}(undef, 3, 3)

    k = 1
    p = zeros(6)
    for i in 1:3
        for j in 1:i
            p[k] = 1.0
            Ada[i, j] = A[i, j] + TaylorPolynomial(k, 1.0)
            Afd[i, j] = Dual(A[i, j], p...)
            p[k] = 0.0
            k += 1
        end
    end

    for i in 1:3
        for j in 1:(i - 1)
            Ada[j, i] = Ada[i, j]
            Afd[j, i] = Afd[i, j]
        end
    end

    Fda = LinearAlgebra.eigen(Hermitian(Ada))
    Ffd = LinearAlgebra.eigen(Hermitian(Afd))

    λda, Vda = Fda.values, Fda.vectors          # GenericLinearAlgebra + TaylorPolynomial
    λcc, Vcc = DifferentialAlgebra.eigenpairs(Ada)                   # Taylor eigenpairs in Julia
    λfd, Vfd = Ffd.values, Ffd.vectors          # GenericLinearAlgebra + ForwardDiff

    cλda = DifferentialAlgebra.constant_term.(λda)
    cλcc = DifferentialAlgebra.constant_term.(λcc)
    cλfd = [λ.value for λ in λfd]

    @testset verbose = true "Eigenvalues constant part" begin
        @test all(isapprox.(cλda, cλcc, atol = 1.0e-14, rtol = 1.0e-14))
        @test all(isapprox.(cλda, cλfd, atol = 1.0e-14, rtol = 1.0e-14))
    end

    Jλda = vcat([DifferentialAlgebra.linear_part.(λda)...]'...)
    Jλcc = vcat([DifferentialAlgebra.linear_part.(λcc)...]'...)
    Jλfd = vcat([[λ.partials...] for λ in λfd]'...)

    @testset verbose = true "Eigenvalues partials" begin
        @test all(isapprox.(Jλda, Jλcc, atol = 1.0e-13, rtol = 1.0e-13))
        @test all(isapprox.(Jλda, Jλfd, atol = 1.0e-14, rtol = 1.0e-14))
    end

    cVda = DifferentialAlgebra.constant_term.(Vda)
    cVcc = DifferentialAlgebra.constant_term.(Vcc)
    cVfd = [V.value for V in Vfd]

    scc = [sign.(cVda[1, i]) == sign.(cVcc[1, i]) ? 1 : -1 for i in 1:3]'
    sfd = [sign.(cVda[1, i]) == sign.(cVfd[1, i]) ? 1 : -1 for i in 1:3]'

    @testset verbose = true "Eigenvectors constant part" begin
        @test all(isapprox.(cVda, cVcc .* scc, atol = 1.0e-14, rtol = 1.0e-14))
        @test all(isapprox.(cVda, cVfd .* sfd, atol = 1.0e-14, rtol = 1.0e-14))
    end

    JVda = [[DifferentialAlgebra.linear_part(V)...] for V in Vda]
    JVcc = [[DifferentialAlgebra.linear_part(V)...] for V in Vcc]
    JVfd = [[V.partials...] for V in Vfd]

    @testset verbose = true "Eigenvectors partials" begin
        for i in 1:3
            for j in 1:3
                @test all(isapprox.(JVda[i, j], JVcc[i, j] .* scc[j], atol = 1.0e-13, rtol = 1.0e-13))
                @test all(isapprox.(JVda[i, j], JVfd[i, j] .* sfd[j], atol = 1.0e-13, rtol = 1.0e-13))
            end
        end
    end
end

@testset "Eigenpairs use the working order" begin
    for budget in (0, 32 * 1024^2)
        x, = variables(1; order = 5, table_bytes = budget)
        diagonal = [2 + x + x^4, 4 - x + x^5]
        # Above-order off-diagonal entries must not disable the diagonal path.
        A = [diagonal[1] x^5; zero(x) diagonal[2]]
        values, vectors = with_order(() -> DifferentialAlgebra.eigenpairs(A), 2)
        @test constant_term(values) == [2, 4]
        @test coefficient.(values, Ref([1])) == [1, -1]
        @test all(p -> degree(p) <= 2, values)
        @test constant_term(vectors) == Matrix{Float64}(I, 2, 2)
        @test values[1] !== diagonal[1]
        set_coefficient!(values[1], [1], 7)
        @test coefficient(diagonal[1], [1]) == 1
        @test coefficient(diagonal[1], [4]) == 1
        @test truncation_order() == 5
    end
end

# Every matrix access checks the public order, including accesses inside the
# eigenpair lifting loop. This deterministically detects hidden set_truncation_order! calls.
struct OrderCheckedMatrix{T} <: AbstractMatrix{T}
    data::Matrix{T}
    order::Int
end
Base.size(A::OrderCheckedMatrix) = size(A.data)
function Base.getindex(A::OrderCheckedMatrix, i::Int, j::Int)
    truncation_order() == A.order || error("eigenpairs changed the working order")
    return A.data[i, j]
end

@testset "Eigenpair precision and working-order scope" begin
    for T in (Float32, Float64, BigFloat), budget in (0, 32 * 1024^2)
        x, y, z = variables(T, 3; order = 6, table_bytes = budget)
        tolerance = T === Float32 ? T(2.0e-5) :
            T === Float64 ? T(2.0e-13) : T(2)^(-200)
        A = OrderCheckedMatrix([2 + x / 10 one(x) / 5; one(x) / 5 4 + y / 10], 6)
        values, vectors = DifferentialAlgebra.eigenpairs(A)
        @test eltype(values) === TaylorPolynomial{T}
        @test coefficientwise_error(A * vectors, vectors * Diagonal(values)) < 10tolerance
        @test coefficientwise_error(vectors' * vectors, Matrix{T}(I, 2, 2)) < 10tolerance
        @test truncation_order() == 6
    end
end

@testset "Higher-order symmetric eigenpairs" begin
    x, y = variables(2; order = 4)
    set_coefficient_tolerance!(1.0e-24)
    A = [2 + x 0.3 + y 0.2x; 0.3 + y 4 - y 0.1 + x * y; 0.2x 0.1 + x * y 7 + x + y]
    values, vectors = DifferentialAlgebra.eigenpairs(A)
    @test coefficientwise_error(A * vectors, vectors * Diagonal(values)) < 2.0e-12
    @test coefficientwise_error(
        transpose(vectors) * vectors,
        TaylorPolynomial.(Matrix{Float64}(I, 3, 3))
    ) < 2.0e-12
    @test DifferentialAlgebra.truncation_order() == 4
    # Evaluate the Taylor eigenvalues against independent numeric LAPACK.
    for point in ([1.0e-3, 2.0e-3], [-2.0e-3, 1.0e-3])
        numeric = [DifferentialAlgebra.evaluate(a, point) for a in A]
        @test evaluate(values, point) ≈ eigvals(Symmetric(numeric)) atol = 2.0e-13
    end
    repeated, basis = DifferentialAlgebra.eigenpairs(
        [1 + x TaylorPolynomial(); TaylorPolynomial() 1 + y]
    )
    @test coefficientwise_error(repeated, [1 + x, 1 + y]) == 0
    @test coefficientwise_error(basis, TaylorPolynomial.(Matrix{Float64}(I, 2, 2))) == 0
    @test_throws ArgumentError DifferentialAlgebra.eigenpairs([1 + x y; y 1 - x])
    @test_throws ArgumentError DifferentialAlgebra.eigenpairs([x y; x y])
end
