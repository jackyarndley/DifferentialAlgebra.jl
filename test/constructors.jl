using DifferentialAlgebra,Test

@testset "Legacy default and copy constructors" begin
    for T in (Float32,Float64,BigFloat)
        DifferentialAlgebra.init(4,2)
        m = Monomial{T}()
        @test DifferentialAlgebra.getCoefficient(m) isa T && iszero(DifferentialAlgebra.getCoefficient(m))
        @test DifferentialAlgebra.getExponents(m) == UInt32[0,0]
        @test DifferentialAlgebra.order(Monomial()) == 0
        @test size(AlgebraicMatrix{DA{T}}()) == (0,0)
        x,y = DifferentialAlgebra.identity(T)
        map = DifferentialAlgebra.compile([1+x,2+y])
        for other in (copy(map),DifferentialAlgebra.compiledDA(map))
            @test other.coefficients !== map.coefficients
            @test other.levels !== map.levels && other.indices !== map.indices
            @test DifferentialAlgebra.evaluate(other,T[2,3]) == T[3,5]
        end
        other = copy(map)
        map.coefficients[1,1] = T(10)
        @test DifferentialAlgebra.evaluate(other,T[2,3]) == T[3,5]
        DifferentialAlgebra.init(2,1)
        @test DifferentialAlgebra.evaluate(other,T[2,3]) == T[3,5]
    end
end
