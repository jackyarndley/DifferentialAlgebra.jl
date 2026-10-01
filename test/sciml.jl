import DiffEqBase

@testset "SciML scalar values and types" begin
    for T in (Float32, Float64, BigFloat, Rational{Int})
        x = only(variables(T, 1; order = 3))
        p = T(2) + x + x^2
        @test @inferred(DiffEqBase.value(p)) === constant_term(p)
        @test @inferred(DiffEqBase.value(typeof(p))) === T
    end
end
