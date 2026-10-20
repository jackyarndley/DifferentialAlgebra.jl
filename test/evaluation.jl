using Test, DifferentialAlgebra
using ForwardDiff, IntervalArithmetic

include("support/intervals.jl")
using .IntervalTestSupport: IA, interval_contains

@testset "Inverse construction uses the input's algebra" begin
    for T in (Float64, BigFloat, Rational{BigInt}), budget in (0, 32 * 1024^2)
        root, = variables(1; order = 2)
        DifferentialAlgebra.with_algebra(5, 2; table_bytes = budget) do _
            x, y = variable(1, T), variable(2, T)
            full, partial = [2x + y^2, y], [2x + x^2]
            result = fetch(
                Threads.@spawn begin
                    complete = invert(full)
                    incomplete = only(invert(partial))
                    (; complete, incomplete, dimensions = nvariables(), order = max_order())
                end
            )
            @test result.dimensions == 1 && result.order == 2
            @test all(p -> p.algebra === x.algebra, result.complete)
            @test result.incomplete.algebra === x.algebra
            @test isequal(result.complete, [x / 2 - y^2 / 2, y])
            # Inverse of 2x+x^2 is sqrt(1+u)-1; use exact binomial coefficients.
            for k in 0:5
                expected = k == 0 ? 0 : (1 // 2, -1 // 8, 1 // 16, -5 // 128, 7 // 256)[k]
                @test coefficient(result.incomplete, [k, 0]) == expected
            end
            @test coefficient(partial[1], [1, 0]) == 2
        end
        @test variable(1).algebra === root.algebra && max_order() == 2
    end
end

@testset "Compact polynomial linear transforms" begin
    for budget in (0, 32 * 1024^2)
        x = variables(6; order = 10, table_bytes = budget)
        result = DifferentialAlgebra._linear_transform([2 3; 0 4], [1 + x[1], x[6]])
        @test iszero(result[1] - (2 + 2x[1] + 3x[6]))
        @test iszero(result[2] - 4x[6])
        @test all(p -> length(p.coeffs) <= 7, result)
    end
end

@testset "Dense evaluation uses the stored basis traversal" begin
    for budget in (0, 32 * 1024^2)
        x = variables(6; order = 10, table_bytes = budget)
        p = (1 + sum(x) / 8)^10
        point = fill(0.1, 6)
        @test @inferred(p(point)) ≈ (1 + sum(point) / 8)^10 rtol = 80eps(Float64)
        @test p(Float64[]) == 1
        @test p([1.0]) ≈ (1 + 1 / 8)^10
        @test p([point; 100]) == p(point)
    end
end

@testset "Direct numeric evaluation" begin
    for T in (Float32, Float64, BigFloat, Rational{BigInt}), budget in (0, 32 * 1024^2)
        x, y, z = variables(T, 3; order = 4, table_bytes = budget)
        p = (1 + x / 2 + y / 4 + z / 8)^4
        for point in (T[], T[1 // 4], T[1 // 4, -1 // 2, 1], T[1, 2, 3, 100])
            padded = vcat(point, zeros(T, 3))
            expected = (1 + padded[1] / 2 + padded[2] / 4 + padded[3] / 8)^4
            @test @inferred(p(point)) ≈ expected
        end
        @test p(1 // 4) ≈ (1 + T(1 // 8))^4
        @test coefficient_type(p) === T
        with_order(1) do
            @test p(T[1, 2, 3]) ≈ (1 + T(1) / 2 + T(2) / 4 + T(3) / 8)^4
        end
        set_coefficient_tolerance!(0.1)
        @test p(T[1, 2, 3]) ≈ (1 + T(1) / 2 + T(2) / 4 + T(3) / 8)^4
        setprecision(256) do
            point = [BigFloat(1) / 3, BigFloat(1) / 7, BigFloat(1) / 9]
            expected = (1 + point[1] / 2 + point[2] / 4 + point[3] / 8)^4
            @test p(point) isa BigFloat
            @test p(point) ≈ expected rtol = 32eps(BigFloat)
        end
    end
    x, y = variables(2; order = 3)
    p = 2 + x^2 + x * y / 4 - y
    @test ForwardDiff.gradient(p, [0.25, -0.5]) == [3 / 8, -15 / 16]
    @test ForwardDiff.hessian(p, [0.25, -0.5]) == [2 1 / 4; 1 / 4 0]
    @test p(Real[0.25, -0.5]) == p([0.25, -0.5])
    @test TaylorPolynomial(7)(BigFloat[]) isa BigFloat
    @test (1 + x)([2.0, Inf]) == 3
    @test (1 + x * y)([Inf]) == 1
end

@testset "Direct evaluation preserves interval bounds and metadata" begin
    for T in (Float64, BigFloat), budget in (0, 32 * 1024^2)
        I = Interval{T}
        x, y = variables(I, 2; order = 3, table_bytes = budget)
        p = (1 + x / 2 + y / 4)^3
        box = [IA.interval(T, 1 // 4, 1 // 2), IA.interval(T, -1 // 2, -1 // 4)]
        result = p(box)
        # Monotonicity of the positive affine base gives this independent range.
        @test interval_contains(result, 1)
        @test interval_contains(result, (1 + 1 // 4 - 1 // 16)^3)
        @test IA.isguaranteed(result) && IA.decoration(result) == IA.com
        degraded = IA.setdecoration(box[1], IA.def)
        result = p([degraded, box[2]])
        @test IA.decoration(result) == IA.def && IA.isguaranteed(result)
        zero_coefficient = IA.setdecoration(IA.interval(T, 0), IA.def)
        owned = TaylorPolynomial{I}([zero_coefficient, IA.interval(T, 1)], 2, x.algebra)
        @test IA.decoration(owned([IA.interval(T, 1)])) == IA.def
        nonguaranteed = convert(I, T(0))
        owned = TaylorPolynomial{I}([nonguaranteed, IA.interval(T, 1)], 2, x.algebra)
        @test !IA.isguaranteed(owned([IA.interval(T, 1)]))
    end
end
