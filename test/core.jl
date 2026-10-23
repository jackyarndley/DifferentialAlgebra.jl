using Test, DifferentialAlgebra, Random, LinearAlgebra
using ..Fixtures
const DA = DifferentialAlgebra

@testset "Exact products, square path and working order" begin
    rng = MersenneTwister(7123)
    # A Cartesian dictionary convolution is independent of both basis-index kernels.
    indices = [(i, j) for i in 0:3 for j in 0:(3 - i)]
    a, b = rand(rng, -4:4, length(indices)), rand(rng, -4:4, length(indices))
    for T in (Float64, Rational{BigInt}), budget in (0, 32 * 1024^2)
        x, y = variables(T, 2; order = 3, table_bytes = budget)
        p = sum(c * x^i * y^j for (c, (i, j)) in zip(a, indices))
        q = sum(c * x^i * y^j for (c, (i, j)) in zip(b, indices))
        for (right, weights) in ((q, b), (p, a)), order in (2, 3)
            oracle = Dict(alpha => big(0) // 1 for alpha in indices)
            for (ac, alpha) in zip(a, indices), (bc, beta) in zip(weights, indices)
                gamma = alpha .+ beta
                sum(gamma) <= order && (oracle[gamma] += ac * bc)
            end
            with_order(order) do
                out = copy(p)
                mul!(out, out, right === p ? out : right)
                for result in (p * right, out)
                    @test all(coefficient(result, collect(k)) == v for (k, v) in oracle)
                end
            end
        end
    end
    x, y = variables(2; order = 2)
    @test coefficient((1.0e308x + 1.0e-308y)^2, [1, 1]) ≈ 2
end

@testset "Validation before indexing and conversion" begin
    x, y = variables(2; order = 3)
    for bad in ([1], [1, 0, 0]), f in (coefficient, differentiate, integrate)
        @test_throws DimensionMismatch f(x, bad)
    end
    for bad in ([-1, 0], [big(typemax(Int)) + 1, 0], ShiftedVector([1, 0]))
        for f in (coefficient, differentiate, integrate)
            @test_throws ArgumentError f(x, bad)
        end
        @test_throws ArgumentError set_coefficient!(copy(x), bad, 2)
    end
    @test_throws ArgumentError x(ShiftedVector([1, 2]))
    @test_throws ArgumentError TaylorPolynomial{Real}(1)
    @test_throws ArgumentError TimeSeries(Real[1, 2])
    @test_throws InexactError Float64(1 + x)
    @test_throws InexactError Bool(1 + x)
    @test_throws MethodError round(1 + x)
    @test_throws MethodError trunc(1 + x)
    @test (2 + x + 3y)([4]) == 6
    @test (2 + x + 3y)([4, 5, 99]) == 21
    @test_throws ArgumentError variables(0; order = 3)
    @test x([2]) == 2 # Failed setup leaves the owning algebra valid.
end

@testset "Scalar types, stable norms and exact affine factors" begin
    for T in (Int, Float32, Float64, BigFloat)
        x, y = variables(T, 2; order = 2)
        p = 3x + 4y
        @test coefficient_norm(p, 2) isa typeof(float(one(T)))
        @test coefficient_norm(p, 2) == 5
        @test DA.degree_norms(p, 0, 2) == [0, 5, 0]
        @test typeof(float(p)) == float(typeof(p))
        @test Monomial{T}(2, [1, 2]).coefficient == 2
    end
    x, y = variables(2; order = 2)
    for s in (1.0e-300, 1.0e300)
        @test coefficient_norm(3s * x + 4s * y, 2) / s ≈ 5
    end
    x, = variables(BigInt, 1; order = 60)
    translated = DA.translate_variable(x^60, 1, big(2), big(3))
    @test all(
        coefficient(translated, [k]) == binomial(big(60), k) *
            big(2)^k * big(3)^(60 - k) for k in (0, 1, 30, 59, 60)
    )
end

@testset "Stored identity, copy and inspection" begin
    x, y = variables(2; order = 4)
    p = 2 + x + 3x * y
    @test isequal(constant_term(Real[p, 3]), [2.0, 3.0])
    @test p == 2 && !isequal(p, 2)
    @test length(Set([p, copy(p), 2 + y])) == 2
    q = copy(p)
    set_coefficient!(q, [1, 0], 99)
    @test coefficient(p, [1, 0]) == 1
    @test exponents(Monomial(3, [2, 1])) == [2, 1]
    @test coefficient(DA.filter_terms(p, x + y), [1, 0]) == 1
    @test DA.coefficient_dot(p, x + 2x * y) == 7
    @test occursin("x₁x₂", string(x * y))
    @test_throws ArgumentError coefficient_norm(p, -1)
    @test_throws ArgumentError Monomial{Int}(1, [-1])
    @test_throws ArgumentError Monomial{Int}(1, [big(typemax(UInt32)) + 1])
    exps = [1, 2]
    m = Monomial(2, exps)
    exps[1] = 99
    @test exponents(m) == [1, 2]
    # Narrow the final result only, not the operands of a reusable buffer.
    out = TaylorPolynomial{Float32}(0)
    DA.weighted_sum!(
        out, TaylorPolynomial(16777217.0), 1.0,
        TaylorPolynomial(-16777216.0), 1.0
    )
    @test constant_term(out) == 1.0f0
    setprecision(512) do
        tolerance = big"1e-400"
        set_coefficient_tolerance!(tolerance)
        @test coefficient_tolerance() == tolerance
    end
end
