using SpecialFunctions

@testset "Test erf" begin
    DifferentialAlgebra.initialize!(1, 3)

    x = DifferentialAlgebra.random_polynomial()
    @test isapprox(DifferentialAlgebra.constant_term(erf(x)), erf(DifferentialAlgebra.constant_term(x)), atol = 1.0e-15, rtol = 1.0e-15)

end

@testset "Test erfc" begin
    DifferentialAlgebra.initialize!(1, 3)

    x = DifferentialAlgebra.random_polynomial()
    @test isapprox(DifferentialAlgebra.constant_term(erfc(x)), erfc(DifferentialAlgebra.constant_term(x)), atol = 1.0e-15, rtol = 1.0e-15)

end
