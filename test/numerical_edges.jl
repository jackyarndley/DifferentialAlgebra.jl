using Test, DifferentialAlgebra, IntervalArithmetic
include("support/intervals.jl")
using .IntervalTestSupport: IA, sameinterval, interval_contains

@testset "Exact affine binomial factors" begin
    for T in (BigInt, Rational{BigInt}), budget in (0, 32 * 1024^2)
        x, = variables(T, 1; order = 60, table_bytes = budget)
        # Degree 60 has binomial coefficients beyond exact Float64 integers.
        p = DifferentialAlgebra.translate_variable(x^60, 1, T(2), T(3))
        @test coefficient_type(p) === T
        for k in 0:60
            expected = binomial(big(60), k) * big(2)^k * big(3)^(60 - k)
            @test coefficient(p, [k]) == expected
        end
        @test coefficient(DifferentialAlgebra.scale_variable(x^60, 1, T(2)), [60]) == big(2)^60
    end
    for T in (Float64, BigFloat)
        x, = variables(Interval{T}, 1; order = 8)
        p = DifferentialAlgebra.translate_variable(x^8, 1, 2, 3)
        for k in 0:8
            c = coefficient(p, [k])
            @test interval_contains(c, binomial(8, k) * 2^k * 3^(8 - k))
            @test IA.isguaranteed(c) && IA.decoration(c) == IA.com
        end
        uncertain = IA.interval(T, 2, 3) * x
        c = coefficient(DifferentialAlgebra.translate_variable(uncertain, 1, 1, 1), [0])
        @test sameinterval(c, IA.interval(T, 2, 3))
    end
end

@testset "Requested central moment order" begin
    x, = variables(1; order = 6)
    # N(2,3) has centered moments 0,3,0,27,0,405.
    mgf = exp(2x + 3x^2 / 2)
    with_order(2) do
        @test_throws ArgumentError DifferentialAlgebra.central_moments(mgf, 6)
        @test truncation_order() == 2
        indices, moments = DifferentialAlgebra.central_moments(mgf, 2)
        @test only(last(indices)) == 2
        @test moments ≈ [1, 0, 3]
    end
    indices, moments = DifferentialAlgebra.central_moments(mgf, 6)
    @test length(indices) == 7
    @test moments ≈ [1, 0, 3, 0, 27, 0, 405]
    @test truncation_order() == 6
    @test_throws ArgumentError DifferentialAlgebra.central_moments(2mgf, 2)
    @test_throws ArgumentError DifferentialAlgebra.central_moments(mgf, -1)
end
