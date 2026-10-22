using Test, DifferentialAlgebra

import DiffEqBase
using ForwardDiff: Dual

@testset "SciML scalar values and types" begin
    for T in (Float32, Float64, BigFloat, Rational{Int})
        x = only(variables(T, 1; order = 3))
        p = T(2) + x + x^2
        @test @inferred(DiffEqBase.value(p)) === constant_term(p)
        @test @inferred(DiffEqBase.value(typeof(p))) === T
    end
end

@testset "SciML scalar values recursively unwrap coefficients" begin
    variables(1; order = 2)
    for value in (Dual(2.0, 1.0), Dual(Dual(2.0, 1.0), Dual(1.0, 0.0)))
        p = TaylorPolynomial{typeof(value)}(value)
        @test @inferred(DiffEqBase.value(p)) === 2.0
        @test @inferred(DiffEqBase.value(typeof(p))) === Float64
    end
end
