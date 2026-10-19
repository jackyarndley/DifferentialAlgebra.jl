using Test, DifferentialAlgebra
@testset "Hessians from independent polynomial derivatives" begin
    x, y, z = variables(3; order = 4)
    p, q = x^2 * y + 3x * y * z^2, y^2 + x * z
    expected_p = [2y 2x + 3z^2 6y * z; 2x + 3z^2 zero(x) 6x * z; 6y * z 6x * z 6x * y]
    expected_q = [zero(x) zero(x) one(x); zero(x) 2one(x) zero(x); one(x) zero(x) zero(x)]
    @test all(iszero, hessian(p) - expected_p)
    map = [p, q]
    result = hessian(view(map, :))
    @test all(iszero, result[1] - expected_p)
    @test all(iszero, result[2] - expected_q)
    tensor = DifferentialAlgebra.hessian_tensor(view(map, :))
    @test all(iszero, tensor[:, :, 1] - expected_p)
    @test all(iszero, tensor[:, :, 2] - expected_q)
    # Symmetric entries and source polynomials own independent coefficients.
    set_coefficient!(result[1][1, 2], [0, 0, 0], 99)
    @test iszero(result[1][2, 1] - expected_p[2, 1])
    @test iszero(p - (x^2 * y + 3x * y * z^2))
end
