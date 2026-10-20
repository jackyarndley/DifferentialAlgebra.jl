using Test, DifferentialAlgebra
using LinearAlgebra

@testset verbose = true "Comparison operators" begin
    @testset "Test ==" begin
        DifferentialAlgebra.initialize!(10, 4)

        x = 1 + TaylorPolynomial(1, 1.0)
        @test 1 == x

        y = TaylorPolynomial(1)
        @test 1 == y

        @test x == y
    end

    @testset "Test >=" begin
        DifferentialAlgebra.initialize!(10, 4)

        x = 1 + TaylorPolynomial(1, 1)
        y = TaylorPolynomial(2)
        @test y >= x
        @test y >= 1
        @test y >= 2

    end

    @testset "Test isnan" begin
        DifferentialAlgebra.initialize!(2, 2)

        x = 1.0 + TaylorPolynomial(1, 1.0)
        y = 1.0 + TaylorPolynomial(2, Inf)
        z = 1.0 + TaylorPolynomial(2, NaN)
        @test isnan(x) == false
        @test isnan(y) == false
        @test isnan(z) == true
    end

    @testset "Test isinf" begin
        DifferentialAlgebra.initialize!(2, 2)

        x = 1.0 + TaylorPolynomial(1, 1.0)
        y = 1.0 + TaylorPolynomial(2, Inf)
        z = 1.0 + TaylorPolynomial(2, NaN)
        @test isinf(x) == false
        @test isinf(y) == true
        @test isinf(z) == false
    end

    @testset "Test isfinite" begin
        DifferentialAlgebra.initialize!(2, 2)

        x = 1.0 + TaylorPolynomial(1, 1.0)
        y = 1.0 + TaylorPolynomial(2, Inf)
        z = 1.0 + TaylorPolynomial(2, NaN)
        @test isfinite(x) == true
        @test isfinite(y) == false
        @test isfinite(z) == false
    end
end

@testset "Predicates inspect the complete polynomial" begin
    x, y = variables(2; order = 3)
    @test !isone(1 + x)
    @test isone(one(x))
    @test !isinteger(1 + x)
    @test isinteger(TaylorPolynomial(2))
    @test !isinteger(TaylorPolynomial(1.5))
    for scale! in (lmul!, rmul!)
        v = [x, y]
        scale! === lmul! ? scale!(1 + x, v) : scale!(v, 1 + x)
        @test coefficient(v[1], [2, 0]) == 1
        @test coefficient(v[2], [1, 1]) == 1
    end
    A = reshape([x, y], 2, 1)
    B = reshape([one(x)], 1, 1)
    C = [zero(x), zero(x)]
    mul!(reshape(C, 2, 1), A, B, 1 + x, false)
    @test coefficient(C[1], [2, 0]) == 1
    @test coefficient(C[2], [1, 1]) == 1
end
