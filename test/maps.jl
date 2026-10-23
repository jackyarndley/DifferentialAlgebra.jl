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

@testset "Map operations require one owning algebra" begin
    x, y = variables(2; order = 2)
    map = compile([1 + x + 3000y])
    set_coefficient_tolerance!(0.01)
    DA.with_algebra(2, 2) do _
        z = variable(1)
        mixed = [x, z]
        @test_throws ArgumentError linear_part(mixed)
        @test_throws ArgumentError jacobian(mixed)
        @test_throws ArgumentError DA.hessian_tensor(mixed)
        @test_throws ArgumentError DA.eigenpairs([x zero(x); zero(z) z])
        # Zero weights must not bypass ownership checks in the shared transform.
        @test_throws ArgumentError DA._linear_transform(zeros(2, 2), mixed)
        @test coefficients(only(map(Real[x, 0.005]))) == [16, 1, 0, 0, 0, 0]
    end
    @test x([0.5]) == 0.5
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
        @test coefficients(integrate(x + y, [1, 1])) ==
            [0, 0, 0, 0, 0, 0, 0, 1 / 2, 1 / 2, 0, zeros(5)...]
        p = 1.0e-5x^4
        set_coefficient_tolerance!(1.0e-4)
        @test coefficient(differentiate(p, [2, 0]), [2, 0]) ≈ 0.00012
        @test coefficient(hessian(p)[1, 1], [2, 0]) ≈ 0.00012
    end
    x, y = variables(2; order = 3)
    p = 2 + x + 3x * y + y^2
    @test all(samecoefficients.(gradient(p), [1 + 3y, 3x + 2y]))
    @test constant_term.(jacobian([p, x * y])) == [1 0; 0 0]
    tensor = DA.hessian_tensor([p, x * y])
    @test all(iszero, degree.(tensor)) &&
        constant_term.(tensor) == cat([0 3; 3 2], [0 1; 1 0]; dims = 3)
end

@testset "Evaluation, composition and inverse maps" begin
    x, y = variables(2; order = 4)
    p = 2 + x + 3x * y + y^2
    map = CompiledMap([p, x + y])
    @test_throws ArgumentError map(ShiftedVector([2, 3]))
    point = [2.0, 3.0]
    out, work = zeros(2), zeros(nnodes(map))
    evaluate!(out, map, point, work)
    @test out == [31, 5] && p(point) == 31
    @test coefficients.(map([1 + x, 2 + y])) ==
        [[13, 7, 7, 0, 3, 1, zeros(9)...], [3, 1, 1, zeros(12)...]]
    inverse = invert([2 + 2x + y^2, -1 + y])
    @test coefficients.(inverse) ==
        [[-1.5, 0.5, -1, 0, 0, -0.5, zeros(9)...], [1, 0, 1, zeros(12)...]]
    @test_throws DA.TaylorError invert([x, 2x])
    @test_throws DimensionMismatch evaluate!(zeros(1), map, point, work)
    @test_throws DimensionMismatch evaluate!(out, map, point, zeros(degree(map)))
    for (result, coordinates, scratch) in ((out, out, work), (out, point, out), (out, point, point))
        @test_throws ArgumentError evaluate!(result, map, coordinates, scratch)
    end
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

@testset "Public map and moment calculations filter completed results" begin
    x, y = variables(2; order = 2)
    f = [x - y + 0.06x^2, y + 0.06x^2]
    set_coefficient_tolerance!(0.1)
    inverse = invert(f)
    # Quadratic inverse: z - A^-1*q(A^-1*z), independent of DA composition.
    @test coefficient.(inverse, Ref([1, 1])) ≈ [-0.24, -0.12]
    @test [coefficient(inverse[1], a) for a in ([2, 0], [0, 2])] ≈ [-0.12, -0.12]
    x, = variables(1; order = 3)
    mgf = exp(0.01x + 5000x^2)
    set_coefficient_tolerance!(0.1)
    @test last(DA.central_moments(mgf, 3)) ≈ [1, 0, 10000, 0]
end

@testset "Moment scaling avoids factorial overflow" begin
    x, = variables(1; order = 200, table_bytes = 0)
    @test last(last(DA.raw_moments(one(x), 200))) == 0
    p = monomial([200], 1.0e-300)
    expected = Float64(BigFloat(1.0e-300) * factorial(big(200)))
    @test last(last(DA.raw_moments(p, 200))) ≈ expected rtol = 1.0e-13
    for T in (BigInt, Rational{BigInt})
        x, = variables(T, 1; order = 30, table_bytes = 0)
        @test last(last(DA.raw_moments(x^30, 30))) == factorial(big(30))
    end
    x, = variables(Int, 1; order = 21, table_bytes = 0)
    @test_throws OverflowError DA.raw_moments(x^21, 21)
end
