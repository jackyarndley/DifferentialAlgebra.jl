using Test, DifferentialAlgebra

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

@testset "Converting monomials" begin
    for T in (Int, Float32, Float64, BigFloat, Rational{Int})
        m = Monomial{T}(2, [1, 3])
        @test coefficient(m) isa T
        @test coefficient(m) == T(2)
        @test exponents(m) == UInt32[1, 3]
        @test coefficient(Monomial(T(2), view([1, 3], :))) isa T
        @test_throws ArgumentError Monomial{T}(2, [-1])
        @test_throws ArgumentError Monomial{T}(2, [big(typemax(UInt32)) + 1])
    end
    @test_throws ArgumentError Monomial{Real}(2, [1])
    @test_throws InexactError Monomial{Int}(1.5, [1])
end

@testset "Monomial identity, degree and exponent ownership" begin
    x, = variables(1; order = 3)
    @test monomials(1 + x + x^2) == monomials(1 + x + x^2)
    powers = UInt32[1, 2]
    m = Monomial(3.0, powers)
    powers[1] = 9
    @test exponents(m) == UInt32[1, 2]
    @test m == Monomial(3, [1, 2])
    @test isequal(m, Monomial(3, [1, 2]))
    @test hash(m) == hash(Monomial(3, [1, 2]))
    @test length(Set([m, Monomial(3, [1, 2])])) == 1
    @test degree(m) === 3
    @test degree(Monomial(1.0, [0, 0])) - 1 === -1
    @test degree(Monomial(1.0, UInt32[])) === 0
end
