using Test, DifferentialAlgebra

@testset "Polynomial differentiation and integration" begin
    x, y = variables(2; order = 5)
    p = 2 + x + 0.2y
    @test coefficient_norm(differentiate(p^3, [2, 1]) - 1.2) < 1.0e-13
    @test coefficient_norm(differentiate(integrate(p, 1), 1) - p) < 1.0e-14
    @test coefficient_norm(differentiate(integrate(p, [1, 1]), [1, 1]) - p) < 1.0e-14
    @test iszero(DifferentialAlgebra.divide_variable(x^2 * y, 1, 2) - y)
end

@testset "Calculus storage follows the retained input degree" begin
    for budget in (0, 32 * 1024^2), T in (Float64, Rational{BigInt})
        x = variables(T, 6; order = 10, table_bytes = budget)
        p = x[1]^3 + 2x[1] * x[6]
        d, integral = differentiate(p, 1), integrate(p, 1)
        @test iszero(d - (3x[1]^2 + 2x[6]))
        @test iszero(integral - (x[1]^4 / 4 + x[1]^2 * x[6]))
        @test length(d.coeffs) <= binomial(6 + 2, 2)
        @test length(integral.coeffs) <= binomial(6 + 4, 4)
        @test length(differentiate(x[1], 1).coeffs) == 1
        p = x[1]^4
        with_order(2) do
            derivative = differentiate(p, [2, 0, 0, 0, 0, 0])
            @test iszero(derivative - 12x[1]^2)
            @test length(derivative.coeffs) <= binomial(6 + 2, 2)
        end
    end
end

@testset "Gaussian probability integral" begin
    x, = variables(1; order = 24)
    density = exp(-x^2 / 2) / sqrt(2pi)
    primitive = integrate(density, 1)
    # Independent standard normal probability between -1 and 1.
    @test primitive(1.0) - primitive(-1.0) ≈ 0.6826894921370859 atol = 1.0e-12
end
