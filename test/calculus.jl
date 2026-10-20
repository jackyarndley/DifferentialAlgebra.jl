using Test, DifferentialAlgebra

include("support/coefficients.jl")

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

@testset "Gradient components and owning contexts" begin
    for T in (Float32, Float64, BigFloat, Rational{BigInt}), budget in (0, 32 * 1024^2)
        x, y, z = variables(T, 3; order = 4, table_bytes = budget)
        p = 3x^2 * y + 2x * y^2 + 5z
        expected = [6x * y + 2y^2, 3x^2 + 4x * y, 5one(z)]
        # A different current factory must not change an object's variable count.
        result = DifferentialAlgebra.with_algebra(2, 1) do _
            gradient(p)
        end
        @test size(DifferentialAlgebra.linear_part(TaylorPolynomial{T}[])) == (0, 3)
        @test size(DifferentialAlgebra.hessian_tensor(TaylorPolynomial{T}[])) == (3, 3, 0)
        @test length(result) == 3
        @test all(iszero, result - expected)
        @test all(d -> d.algebra === p.algebra, result)
        @test all(i -> result[i].coeffs !== p.coeffs, eachindex(result))
        @test result[1].coeffs !== result[2].coeffs
        snapshot = copy(p.coeffs)
        set_coefficient!(result[1], [0, 0, 0], 11)
        @test p.coeffs == snapshot
        with_order(1) do
            g = gradient(p)
            @test all(iszero, g[1:2])
            @test iszero(g[3] - 5)
        end
    end
end

@testset "Map Jacobian and Hessian shapes" begin
    x, y = variables(2; order = 4)
    @test size(DifferentialAlgebra.hessian_tensor([x * x + y, x * y])) == (2, 2, 2)
    @test coefficientwise_error(DifferentialAlgebra.hessian(x * x + y)[1, 1], 2) == 0
    @test coefficientwise_error(
        DifferentialAlgebra.jacobian([x * x + y, x * y])[1, 2],
        1
    ) == 0
end
