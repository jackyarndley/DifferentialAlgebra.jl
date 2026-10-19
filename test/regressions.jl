using Test, DifferentialAlgebra
using LinearAlgebra, SpecialFunctions, ForwardDiff, GenericLinearAlgebra

@testset "Mixed precision reusable buffers" begin
    initialize!(1, 1)
    a, b = TaylorPolynomial(16777217.0), TaylorPolynomial(-16777216.0)
    out = TaylorPolynomial{Float32}(0)
    # Preserve operand precision before the final assignment to the buffer.
    # Narrowing either operand first loses this exactly representable result.
    DifferentialAlgebra.weighted_sum!(out, a, 1.0, b, 1.0)
    @test constant_term(out) == 1.0f0
    DifferentialAlgebra.scale!(out, TaylorPolynomial(1.0e40), 1.0e-40)
    @test constant_term(out) == 1.0f0
end

@testset "Working-order arithmetic" begin
    x, = variables(1; order = 3)
    p = x^3
    with_order(1) do
        @test (degree(p + 1), degree(p + TaylorPolynomial(1.0)), degree(1 * p)) == (0, 0, 0)
        for q in (1 + p, p - 1, 1 - p, +p, -p, p^1, p / 1, muladd(1, p, one(p)))
            @test degree(q) <= 1
        end
        @test degree(copy(p)) == 3
        @test degree(TaylorPolynomial{Float32}(p)) == 3
        @test degree(p) == 3
    end
    @test truncation_order() == 3
    @test degree(p + 1) == 3
end

# Every matrix access checks the public order, including accesses inside the
# eigenpair lifting loop. This deterministically detects hidden set_truncation_order! calls.
struct OrderCheckedMatrix{T} <: AbstractMatrix{T}
    data::Matrix{T}
    order::Int
end
Base.size(A::OrderCheckedMatrix) = size(A.data)
function Base.getindex(A::OrderCheckedMatrix, i::Int, j::Int)
    DifferentialAlgebra.truncation_order() == A.order || error("eigenpairs changed the global truncation order")
    return A.data[i, j]
end

@testset "Substitution, analytic domains and eigenpairs" begin
    for T in (Float32, Float64, BigFloat), budget in (0, 32 * 1024^2)
        DifferentialAlgebra.initialize!(6, 3; table_bytes = budget)
        x, y, z = variable.(1:DifferentialAlgebra.nvariables(), T)
        p = 3 + x + 2y + x^2 * y + 3x * z^2 + x^5
        tolerance = T === Float32 ? T(2.0e-5) : T === Float64 ? T(2.0e-13) : T(2)^(-200)
        @test DifferentialAlgebra.coefficient_norm(DifferentialAlgebra.substitute(p, 1, 2) - (5 + 6y + 6z^2 + 32)) == 0
        @test DifferentialAlgebra.coefficient_norm(DifferentialAlgebra.substitute(p, 1) - (3 + 2y)) == 0
        @test DifferentialAlgebra.coefficient_norm(DifferentialAlgebra.replace_variable(p, 1, 2, 2) - (3 + 4y + 4y^3 + 6y * z^2 + 32y^5)) == 0
        for v in 1:3
            args = [x, y, z]
            args[v] = 2args[v] + 1
            @test DifferentialAlgebra.coefficient_norm(DifferentialAlgebra.translate_variable(p, v, 2, 1) - DifferentialAlgebra.evaluate(p, args)) < tolerance
            args[v] = -[x, y, z][v]
            @test DifferentialAlgebra.coefficient_norm(DifferentialAlgebra.scale_variable(p, v, -1) - DifferentialAlgebra.evaluate(p, args)) < tolerance
        end
        @test_throws ArgumentError DifferentialAlgebra.substitute(p, 0)
        @test_throws ArgumentError DifferentialAlgebra.replace_variable(p, 1, 4)
        @test DifferentialAlgebra.coefficient_dot(p, 2x + 3y + x^2 * y) == T(9)
        @test DifferentialAlgebra.coefficient_norm(DifferentialAlgebra.filter_terms(p, x + y) - (x + 2y)) == 0
        @test DifferentialAlgebra.nterms(p) == 6
        @test DifferentialAlgebra.coefficient_norm(p) == 3
        @test DifferentialAlgebra.coefficient_norm(DifferentialAlgebra.monomial([2, 1, 0], T(3)) - 3x^2 * y) == 0
        q = inv(1 - x / 4)
        estimates, errors = DifferentialAlgebra.estimate_norms(q, 0, 0, 8; errors = true)
        @test maximum(abs.(estimates - [T(1) / 4^i for i in 0:8])) < tolerance
        @test all(>=(0), errors) && length(errors) == 7
        @test DifferentialAlgebra.convergence_radius(q, T(1) / 2^14) ≈ one(T) rtol = tolerance
        @test_throws ArgumentError DifferentialAlgebra.convergence_radius(q, -1)
        @test DifferentialAlgebra.isinitialized()
        DifferentialAlgebra.set_truncation_order!(2)
        with_order(max_order()) do
            @test truncation_order() == 6
        end
        @test DifferentialAlgebra.truncation_order() == 2
        DifferentialAlgebra.set_truncation_order!()
        @test DifferentialAlgebra.truncation_order() == 6
        @test DifferentialAlgebra.coefficient_norm((p / q) * q - p) < tolerance
        @test DifferentialAlgebra.coefficient_norm(muladd(T(2), p, q) - (2p + q)) < tolerance
        @test DifferentialAlgebra.coefficient_norm(muladd(p, T(2), q) - (2p + q)) < tolerance
        @test DifferentialAlgebra.coefficient_norm(muladd(p, q, p) - (p * q + p)) < tolerance
        @test DifferentialAlgebra.coefficient_norm(sqrt(2 + x + y * z)^2 - (2 + x + y * z)) < tolerance
        @test asin(TaylorPolynomial{T}(1)) isa TaylorPolynomial{T}
        @test DifferentialAlgebra.constant_term(acos(TaylorPolynomial{T}(-1))) == acos(-one(T))
        @test iszero(acosh(TaylorPolynomial{T}(1)))
        @test_throws DifferentialAlgebra.TaylorError asin(1 + x)
        A = OrderCheckedMatrix([2 + x / 10 one(x) / 5;one(x) / 5 4 + y / 10], 6)
        values, vectors = DifferentialAlgebra.eigenpairs(A)
        @test maximum(DifferentialAlgebra.coefficient_norm, A * vectors - vectors * Diagonal(values)) < 20tolerance
        @test maximum(DifferentialAlgebra.coefficient_norm, vectors' * vectors - I) < 20tolerance
        @test isempty(first(DifferentialAlgebra.eigenpairs(Matrix{TaylorPolynomial{T}}(undef, 0, 0))))
        @test size(DifferentialAlgebra.linear_part(TaylorPolynomial{T}[])) == (0, 3)
        @test size(DifferentialAlgebra.hessian_tensor(TaylorPolynomial{T}[])) == (3, 3, 0)
    end

    DifferentialAlgebra.initialize!(5, 2)
    x, y = variable.(1:DifferentialAlgebra.nvariables())
    @test DifferentialAlgebra.constant_term(hypot(TaylorPolynomial(3.0e200), TaylorPolynomial(4.0e200))) ≈ 5.0e200
    @test DifferentialAlgebra.constant_term(hypot(TaylorPolynomial(3.0e-200), TaylorPolynomial(4.0e-200))) ≈ 5.0e-200
    @test DifferentialAlgebra.coefficient_norm(mod(-1.5 + x, 2) - (0.5 + x)) == 0
    setprecision(256) do
        DifferentialAlgebra.set_coefficient_tolerance!(BigFloat(1) - BigFloat(2)^(-100))
        @test DifferentialAlgebra.constant_term(TaylorPolynomial(1.0)) == 1
        @test DifferentialAlgebra.constant_term(TaylorPolynomial(1.0f0)) == 1
        DifferentialAlgebra.set_coefficient_tolerance!(0)
    end
    @test DifferentialAlgebra.coefficient((1.0e308x + 1.0e-308y)^2, [1, 1]) ≈ 2
    @test_throws ArgumentError DifferentialAlgebra.multiindices(1000, 1000)
end
