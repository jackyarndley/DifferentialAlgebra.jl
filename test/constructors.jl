using DifferentialAlgebra, Test

@testset "Monomial and compiled map constructors" begin
    for T in (Float32, Float64, BigFloat)
        DifferentialAlgebra.initialize!(4, 2)
        m = Monomial{T}()
        @test DifferentialAlgebra.coefficient(m) isa T && iszero(DifferentialAlgebra.coefficient(m))
        @test DifferentialAlgebra.exponents(m) == UInt32[0, 0]
        @test DifferentialAlgebra.degree(Monomial()) == 0
        x, y = variable.(1:DifferentialAlgebra.nvariables(), T)
        map = DifferentialAlgebra.compile([1 + x, 2 + y])
        for other in (copy(map), DifferentialAlgebra.CompiledMap(map))
            @test other.coefficients !== map.coefficients
            @test other.levels !== map.levels && other.indices !== map.indices
            @test DifferentialAlgebra.evaluate(other, T[2, 3]) == T[3, 5]
        end
        other = copy(map)
        map.coefficients[1, 1] = T(10)
        @test DifferentialAlgebra.evaluate(other, T[2, 3]) == T[3, 5]
        DifferentialAlgebra.initialize!(2, 1)
        @test DifferentialAlgebra.evaluate(other, T[2, 3]) == T[3, 5]
    end
end
