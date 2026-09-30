@testset "Public interface and generic coefficients" begin
    x,y = variables(2;order=4)
    @test [x,y] isa Vector{DA{Float64}}
    p = 2+x+3x*y+y*y
    @test coefficient(p,[1,1]) == 3
    @test_throws DimensionMismatch coefficient(p,[1])
    @test constant_term(p) == 2
    @test constant_term([p,x]) == [2,0]
    @test p([2,3]) == 31
    @test DifferentialAlgebra.norm(p([y,x])-(2+y+3x*y+x*x)) == 0
    map = CompiledMap([p,x+y])
    @test map([2,3]) == [31,5]
    @test DifferentialAlgebra.norm(differentiate(p,1)-(1+3y)) == 0
    @test constant_term(differentiate(p,[1,1])) == 3
    @test DifferentialAlgebra.norm(differentiate(integrate(p,1),1)-p) == 0
    @test_throws ArgumentError variables(AbstractFloat,2;order=4)
    @test_throws ArgumentError variables(0;order=4)
    @test_throws ArgumentError variables(2;order=0)
    @test p([2,3]) == 31 # Invalid setup leaves the current algebra intact.
    variables(1;order=2)
    @test_throws ArgumentError p([2,3])
    @test map([2,3]) == [31,5] # Numeric compiled maps own their coefficients.

    # An additional floating-point type exercises the generic scalar kernels.
    for budget in (0,32*1024^2)
        x,y = variables(Float16,2;order=3,table_bytes=budget)
        p = Float16(2)+x/2+y/4
        @test p isa DA{Float16}
        @test p(Float16[1,2]) == Float16(3)
        @test coefficient(p*p,[1,1]) == Float16(1)/4
        @test coefficient(differentiate(p,1),[0,0]) == Float16(1)/2
        @test coefficient(integrate(x,1),[2,0]) == Float16(1)/2
        @test exp(x) isa DA{Float16}
        @test coefficient(exp(x),[2,0]) ≈ Float16(1)/2

        # Exact coefficients are Real but not AbstractFloat.
        x,y = variables(Rational{BigInt},2;order=4,table_bytes=budget)
        p = (1+x+y)^3
        @test p isa DA{Rational{BigInt}}
        @test coefficient(p,[1,1]) == 6
        @test p([1//3,1//4]) == (1+1//3+1//4)^3
        @test inv(1-x) isa DA{Rational{BigInt}}
        @test coefficient(inv(1-x),[4,0]) == 1
        @test inv(DA{Rational{BigInt}}(3)) isa DA{Rational{BigInt}}
        @test constant_term(inv(DA{Rational{BigInt}}(3))) == 1//3
        @test DifferentialAlgebra.norm(differentiate(integrate(x*y,1),1)-x*y) == 0
        @test CompiledMap([p,x*y])([1//3,1//4]) == [p([1//3,1//4]),1//12]
    end

    # Differentiation through the coefficients requires a wider bound than AbstractFloat.
    T = ForwardDiff.Dual{Nothing,Float64,1}
    x, = variables(T,1;order=3)
    a = ForwardDiff.Dual(2.0,1.0)
    p = exp(DA{T}(a)*x)
    c = coefficient(p,[2])
    @test p isa DA{T}
    @test ForwardDiff.value(c) == 2
    @test only(ForwardDiff.partials(c)) == 2
    value = p([0.1])
    @test ForwardDiff.value(value) ≈ 1+0.2+0.2^2/2+0.2^3/6
    @test only(ForwardDiff.partials(value)) ≈ 0.1+2*0.1^2+2*0.1^3
end
