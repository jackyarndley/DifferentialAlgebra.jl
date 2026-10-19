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

@testset "Stored degrees above the working order" begin
    for budget in (0, 32 * 1024^2)
        x, = variables(1; order = 4, table_bytes = budget)
        p = x^4
        positive, negative = 1 + p, -1 - p
        with_order(2) do
            @test coefficient(differentiate(p, [2]), [2]) == 12
            @test coefficient(hessian(p)[1, 1], [2]) == 12
            @test degree(abs(positive)) == degree(abs(negative)) == 0
        end
        x, y = variables(2; order = 4, table_bytes = budget)
        p = x^4 + 2x^2 * y^2
        with_order(2) do
            d = differentiate(p, [2, 0])
            H = hessian(p)
            @test coefficient(d, [2, 0]) == coefficient(H[1, 1], [2, 0]) == 12
            @test coefficient(d, [0, 2]) == coefficient(H[1, 1], [0, 2]) == 4
            @test coefficient(differentiate(p, [1, 1]), [1, 1]) == 8
            @test coefficient(H[1, 2], [1, 1]) == coefficient(H[2, 1], [1, 1]) == 8
            @test truncation_order() == 2
            @test degree(p) == 4
            @test H[1, 2].coeffs !== H[2, 1].coeffs
            @test d.coeffs !== H[1, 1].coeffs
        end
        @test truncation_order() == 4
        # An intermediate derivative below tolerance grows on differentiation.
        p = 1.0e-5x^4
        set_coefficient_tolerance!(1.0e-4)
        @test coefficient(differentiate(p, [2, 0]), [2, 0]) ≈ 0.00012
        @test coefficient(hessian(p)[1, 1], [2, 0]) ≈ 0.00012
        @test coefficient_tolerance() == 1.0e-4
    end
end
