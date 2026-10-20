using Test, DifferentialAlgebra
using IntervalArithmetic

include("support/intervals.jl")
using .IntervalTestSupport: IA, sameinterval, interval_contains

@testset "Affine substitutions retain compact low-degree storage" begin
    for budget in (0, 32 * 1024^2), T in (Float64, Rational{BigInt})
        x = variables(T, 6; order = 10, table_bytes = budget)
        p = 1 + x[1] * x[6]
        translated = DifferentialAlgebra.translate_variable(p, 6, 2, 3)
        @test iszero(translated - (1 + 3x[1] + 2x[1] * x[6]))
        @test length(translated.coeffs) <= binomial(6 + 2, 2)
        substituted = DifferentialAlgebra.substitute(p, 6, 3)
        @test iszero(substituted - (1 + 3x[1]))
        @test length(substituted.coeffs) <= binomial(6 + 2, 2)
    end
end

@testset "Exact affine binomial factors" begin
    for T in (Int, Rational{Int})
        x, = variables(T, 1; order = 60)
        p = DifferentialAlgebra.translate_variable(x^60, 1, 1, 1)
        @test all(k -> coefficient(p, [k]) == binomial(big(60), k), 0:60)
    end
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

@testset "Composition filters only completed outputs" begin
    for budget in (0, 32 * 1024^2)
        x, y, z = variables(3; order = 3, table_bytes = budget)
        small = [0.006x, 0.006x, zero(x)]
        quadratic = [0.02x, zero(x), zero(x)]
        parts = (0.02x, -0.014x, 0.006x)
        p, q = x + y + z, 100x^2
        compiled = compile([p, q])
        set_coefficient_tolerance!(0.01)
        @test coefficient(evaluate(p, small), [1, 0, 0]) ≈ 0.012
        @test coefficient(evaluate(q, quadratic), [2, 0, 0]) ≈ 0.04
        @test coefficient(compiled(quadratic)[2], [2, 0, 0]) ≈ 0.04
        for (i, j, k) in (
                (1, 2, 3), (1, 3, 2), (2, 1, 3),
                (2, 3, 1), (3, 1, 2), (3, 2, 1),
            )
            result = evaluate(p, [parts[i], parts[j], parts[k]])
            @test coefficient(result, [1, 0, 0]) ≈ 0.012
        end
        result = evaluate(p, small)
        set_coefficient!(result, [1, 0, 0], 7)
        @test coefficient(small[1], [1, 0, 0]) == 0.006
        @test coefficient_tolerance() == 0.01
        @test truncation_order() == 3
    end
end

@testset "Substitution across coefficient types and table budgets" begin
    for T in (Float32, Float64, BigFloat), budget in (0, 32 * 1024^2)
        DifferentialAlgebra.initialize!(6, 3; table_bytes = budget)
        x, y, z = variable.(1:DifferentialAlgebra.nvariables(), T)
        p = 3 + x + 2y + x^2 * y + 3x * z^2 + x^5
        tolerance = T === Float32 ? T(2.0e-5) :
            T === Float64 ? T(2.0e-13) : T(2)^(-200)
        @test coefficient_norm(
            DifferentialAlgebra.substitute(p, 1, 2) - (5 + 6y + 6z^2 + 32)
        ) == 0
        @test coefficient_norm(DifferentialAlgebra.substitute(p, 1) - (3 + 2y)) == 0
        @test coefficient_norm(
            DifferentialAlgebra.replace_variable(p, 1, 2, 2) -
                (3 + 4y + 4y^3 + 6y * z^2 + 32y^5)
        ) == 0
        for v in 1:3
            args = [x, y, z]
            args[v] = 2args[v] + 1
            @test coefficient_norm(
                DifferentialAlgebra.translate_variable(p, v, 2, 1) - evaluate(p, args)
            ) < tolerance
            args[v] = -[x, y, z][v]
            @test coefficient_norm(
                DifferentialAlgebra.scale_variable(p, v, -1) - evaluate(p, args)
            ) < tolerance
        end
        @test_throws ArgumentError DifferentialAlgebra.substitute(p, 0)
        @test_throws ArgumentError DifferentialAlgebra.replace_variable(p, 1, 4)
    end
end
