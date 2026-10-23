using Test, DifferentialAlgebra, SpecialFunctions
using ..Fixtures
const DA = DifferentialAlgebra

@testset "Analytic coefficient oracles" begin
    for budget in (0, 32 * 1024^2)
        x, = variables(1; order = 4, table_bytes = budget)
        cases = (
            (exp(x), [1, 1, 1 / 2, 1 / 6, 1 / 24]),
            (inv(1 + x), [1, -1, 1, -1, 1]),
            (log(1 + x), [0, 1, -1 / 2, 1 / 3, -1 / 4]),
            (sin(x), [0, 1, 0, -1 / 6, 0]),
            (cos(x), [1, 0, -1 / 2, 0, 1 / 24]),
            (tan(x), [0, 1, 0, 1 / 3, 0]),
            (sinh(x), [0, 1, 0, 1 / 6, 0]),
            (cosh(x), [1, 0, 1 / 2, 0, 1 / 24]),
            (tanh(x), [0, 1, 0, -1 / 3, 0]),
            (asin(x), [0, 1, 0, 1 / 6, 0]),
            (atan(x), [0, 1, 0, -1 / 3, 0]),
            (asinh(x), [0, 1, 0, -1 / 6, 0]),
            (atanh(x), [0, 1, 0, 1 / 3, 0]),
            (sqrt(1 + x), [1, 1 / 2, -1 / 8, 1 / 16, -5 / 128]),
        )
        for (result, expected) in cases
            @test coefficients(result) ≈ expected atol = 2.0e-15
        end
        @test coefficient(acos(x), [3]) ≈ -1 / 6
        @test coefficient(acosh(2 + x), [1]) ≈ 1 / sqrt(3)
        @test coefficient(cbrt(1 + x), [2]) ≈ -1 / 9
        @test coefficient(DA.nthroot(1 + x, 5), [2]) ≈ -2 / 25
    end
end

@testset "Final-only filtering, including compound operations" begin
    for budget in (0, 32 * 1024^2)
        x, = variables(1; order = 3, table_bytes = budget)
        set_coefficient_tolerance!(1.0e-6)
        @test coefficients(1.0e8 / (1.0e8 + 1.0e6x))[1:2] ≈ [1, -0.01]
        x, = variables(1; order = 10, table_bytes = budget)
        set_coefficient_tolerance!(1.0e-4)
        @test coefficient(exp(-5 + 100x), [10]) ≈ exp(-5) * 100.0^10 / factorial(10)
        x, = variables(1; order = 2, table_bytes = budget)
        set_coefficient_tolerance!(1.0e-4)
        @test coefficient(erf(4 + 100x), [2]) ≈ -0.005079293868746624
        @test coefficient(sqrt(1.0e16 + 1.0e11x), [1]) ≈ 500
        set_coefficient_tolerance!(0.01)
        @test coefficient(log2(2 + 0.018x), [1]) ≈ 0.018 / (2log(2))
        @test coefficient(log10(2 + 0.18x), [1]) ≈ 0.18 / (2log(10))
        @test coefficient(log(1.1, 2 + 0.018x), [1]) ≈ 0.018 / (2log(1.1))
        @test coefficient((10 + 0.09x)^3, [2]) ≈ 0.243
        @test coefficient((0.5 + 0.09x)^(-3), [2]) ≈ 1.5552
        @test coefficient(hypot(1 + 100x, 0.005), [2]) ≈
            100^2 * 0.005^2 / (2hypot(1, 0.005)^3)
        # Scalar logarithmic derivatives, independent of polynomial composition.
        for c in (10.0, -2.5)
            expected = gamma(c) * (digamma(c)^2 + trigamma(c)) * 0.09^2 / 2
            @test coefficient(gamma(c + 0.09x), [2]) ≈ expected
        end
        a, b, q, r = 10.0, 0.09, 2.0, 0.018
        d1, d2 = r * log(a) + q * b / a, 2r * b / a - q * (b / a)^2
        @test coefficient((a + b * x)^(q + r * x), [2]) ≈ a^q * (d1^2 + d2) / 2
        @test coefficient_tolerance() == 0.01 && truncation_order() == 2
    end
end

@testset "Precision, extremes and domains" begin
    for T in (UInt8, UInt64, UInt128, Rational{UInt64})
        x, = variables(T, 1; order = 3)
        @test coefficients(inv(3 + x)) ≈ [1 / 3, -1 / 9, 1 / 27, -1 / 81]
    end
    setprecision(256) do
        x, = variables(1; order = 3)
        power = big(3) / 7
        result = (3 + x)^power
        expected = big(3)^power * power * (power - 1) * (power - 2) / (6big(3)^3)
        @test coefficient_type(result) === BigFloat
        @test coefficient(result, [3]) ≈ expected rtol = 100eps(BigFloat)
    end
    for budget in (0, 32 * 1024^2), scale in (1.0e-170, 1.0e170)
        x, = variables(1; order = 2, table_bytes = budget)
        result = hypot(scale + x, scale + x / 2)
        expected = setprecision(256) do
            s = BigFloat(scale)
            (s / 2)^2 / (2hypot(s, s)^3)
        end
        @test coefficient(result, [2]) ≈ Float64(expected) rtol = 2.0e-14
    end
    x, = variables(1; order = 2)
    p = 2.0^30 + 8x
    s, c = sincos(p)
    @test coefficient(s, [2]) ≈ -32sin(2.0^30)
    @test coefficient(c, [1]) ≈ -8sin(2.0^30)
    for f in (exp, log, sqrt, gamma, p -> besselj(1, p))
        @test isnan(constant_term(f(TaylorPolynomial(NaN))))
    end
    for (f, p) in ((log, -1 + x), (asin, 1 + x), (sqrt, x), (gamma, -2 + x))
        @test_throws DA.TaylorError f(p)
    end
    @test constant_term(asin(one(x))) == pi / 2
    @test iszero(sqrt(zero(x)))
    wrapped = 5 + 2x + 3x^2
    @test all(
        coefficients(p) == [2, 2, 3] for p in
            (mod(wrapped, 3), mod(wrapped, TaylorPolynomial(3)))
    )
    @test_throws ArgumentError mod(wrapped, 3 + x)
    @test_throws ArgumentError mod(3, wrapped)
    with_order(1) do
        @test coefficients(mod(wrapped, 3)) == [2, 2, 0]
    end
end
