using Test, DifferentialAlgebra

@testset "Polynomial differentiation and integration" begin
    x, y = variables(2; order = 5)
    p = 2 + x + 0.2y
    @test coefficient_norm(differentiate(p^3, [2, 1]) - 1.2) < 1.0e-13
    @test coefficient_norm(differentiate(integrate(p, 1), 1) - p) < 1.0e-14
    @test coefficient_norm(differentiate(integrate(p, [1, 1]), [1, 1]) - p) < 1.0e-14
    @test iszero(DifferentialAlgebra.divide_variable(x^2 * y, 1, 2) - y)
end

@testset "Gaussian probability integral" begin
    x, = variables(1; order = 24)
    density = exp(-x^2 / 2) / sqrt(2pi)
    primitive = integrate(density, 1)
    # Independent standard normal probability between -1 and 1.
    @test primitive(1.0) - primitive(-1.0) ≈ 0.6826894921370859 atol = 1.0e-12
end
