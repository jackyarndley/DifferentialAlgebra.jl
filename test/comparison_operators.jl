@testset verbose = true "Comparison operators" begin
    @testset "Test ==" begin
        DifferentialAlgebra.init(10, 4)

        x = 1 + DA(1, 1.0)
        @test 1 == x

        y = DA(1)
        @test 1 == y

        @test x == y
    end

    @testset "Test >=" begin
        DifferentialAlgebra.init(10, 4)

        x = 1 + DA(1, 1)
        y = DA(2)
        @test y >= x
        @test y >= 1
        @test y >= 2

    end

    @testset "Test isnan" begin
        DifferentialAlgebra.init(2, 2)

        x = 1.0 + DA(1, 1.0)
        y = 1.0 + DA(2, Inf)
        z = 1.0 + DA(2, NaN)
        @test isnan(x) == false
        @test isnan(y) == false
        @test isnan(z) == true
    end

    @testset "Test isinf" begin
        DifferentialAlgebra.init(2, 2)

        x = 1.0 + DA(1, 1.0)
        y = 1.0 + DA(2, Inf)
        z = 1.0 + DA(2, NaN)
        @test isinf(x) == false
        @test isinf(y) == true
        @test isinf(z) == false
    end

    @testset "Test isfinite" begin
        DifferentialAlgebra.init(2, 2)

        x = 1.0 + DA(1, 1.0)
        y = 1.0 + DA(2, Inf)
        z = 1.0 + DA(2, NaN)
        @test isfinite(x) == true
        @test isfinite(y) == false
        @test isfinite(z) == false
    end
end
