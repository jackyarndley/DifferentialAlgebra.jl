using Test, DifferentialAlgebra, LinearAlgebra
using ..Fixtures
const DA = DifferentialAlgebra

Base.@noinline function detached_snapshot()
    x, y = variables(2; order = 2)
    map = compile([1 + x * y])
    initialize!(2, 1)
    return map
end

@testset "Compiled ownership releases obsolete algebras" begin
    map = detached_snapshot()
    GC.gc(true)
    @test map.algebra.value === nothing
    @test map([2, 3]) == [7]
end

@testset "Calculus through stored higher degrees" begin
    for budget in (0, 32 * 1024^2)
        x, y = variables(2; order = 4, table_bytes = budget)
        p = x^4 + 2x^2 * y^2
        positive, negative = 1 + p, -1 - p
        with_order(2) do
            H = hessian(p)
            @test coefficient(differentiate(p, [2, 0]), [2, 0]) == 12
            @test samecoefficients(H[1, 1], 12x^2 + 4y^2)
            @test samecoefficients(H[1, 2], 8x * y)
            set_coefficient!(H[1, 2], [0, 0], 9)
            @test coefficient(H[2, 1], [0, 0]) == 0
            @test degree(abs(positive)) == degree(abs(negative)) == 0
        end
        @test degree(p) == 4 && truncation_order() == 4
        @test samecoefficients(differentiate(integrate(x + y, 1), 1), x + y)
        p = 1.0e-5x^4
        set_coefficient_tolerance!(1.0e-4)
        @test coefficient(differentiate(p, [2, 0]), [2, 0]) ≈ 0.00012
        @test coefficient(hessian(p)[1, 1], [2, 0]) ≈ 0.00012
    end
    x, y = variables(2; order = 3)
    p = 2 + x + 3x * y + y^2
    @test all(samecoefficients.(gradient(p), [1 + 3y, 3x + 2y]))
    @test constant_term.(jacobian([p, x * y])) == [1 0; 0 0]
    @test size(DA.hessian_tensor([p, x * y])) == (2, 2, 2)
end

@testset "Evaluation, composition and inverse maps" begin
    x, y = variables(2; order = 4)
    p = 2 + x + 3x * y + y^2
    map = CompiledMap([p, x + y])
    point = [2.0, 3.0]
    out, work = zeros(2), zeros(nnodes(map))
    evaluate!(out, map, point, work)
    @test out == [31, 5] && p(point) == 31
    @test samecoefficients(only(compile([x + y])([x, x])), 2x)
    inverse = invert([2x + y^2, y])
    @test all(samecoefficients.(inverse, [x / 2 - y^2 / 2, y]))
    @test_throws DA.TaylorError invert([x, 2x])
    @test_throws DimensionMismatch evaluate!(zeros(1), map, point, work)
    # A snapshot owns its coefficients and survives global reinitialization.
    set_coefficient!(p, [0, 0], 99)
    @test map(point) == [31, 5]
    initialize!(2, 1)
    @test map(point) == [31, 5]
    @test_throws ArgumentError p(point)
    @test_throws ArgumentError map([variable(1)])
end

@testset "Gaussian moments from exact closed values" begin
    x, y = variables(2; order = 4)
    mgf = exp(2x + y + x^2 / 2 + 2y^2)
    indices, raw = DA.raw_moments(mgf, 4)
    moments = Dict(Tuple(a) => v for (a, v) in zip(indices, raw))
    @test moments[(2, 0)] ≈ 5
    @test moments[(0, 4)] ≈ 73
    @test moments[(2, 2)] ≈ 25
    _, central = DA.central_moments(mgf, 2)
    @test central ≈ [1, 0, 0, 1, 0, 4]
    with_order(1) do
        @test_throws ArgumentError DA.central_moments(mgf, 2)
        @test truncation_order() == 1
    end
end
