@testset "Hessian" begin

    DifferentialAlgebra.initialize!(4, 4)

    v = [DifferentialAlgebra.random_polynomial() for i in 1:4]
    a = @view v[:]

    H1 = DifferentialAlgebra.hessian(a)
    H2 = DifferentialAlgebra.hessian(v)

    H3 = DifferentialAlgebra.hessian_tensor(a)
    H4 = DifferentialAlgebra.hessian_tensor(v)

    @test all(iszero, stack(H1, dims = 3) - stack(H2, dims = 3))
    @test all(iszero, H3 - H4)
    @test all(iszero, stack(H1, dims = 3) - H3)
end
