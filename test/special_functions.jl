using SpecialFunctions

@testset "Test erf" begin
    DifferentialAlgebra.init(1, 3)

    x = DifferentialAlgebra.random(-1)
    @test isapprox(DifferentialAlgebra.cons(erf(x)), erf(DifferentialAlgebra.cons(x)), atol = 1.0e-15, rtol = 1.0e-15)

end

@testset "Test erfc" begin
    DifferentialAlgebra.init(1, 3)

    x = DifferentialAlgebra.random(-1)
    @test isapprox(DifferentialAlgebra.cons(erfc(x)), erfc(DifferentialAlgebra.cons(x)), atol = 1.0e-15, rtol = 1.0e-15)

end
