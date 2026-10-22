using Test, DifferentialAlgebra, ForwardDiff, LinearAlgebra, GenericLinearAlgebra
using ..Fixtures
const DA = DifferentialAlgebra

@testset "Dual coefficients and independent scalar convolution" begin
    x, = variables(1; order = 3)
    p, q, d = 2 + x, 3 + 2x, ForwardDiff.Dual(2.0, 3.0)
    @test ForwardDiff.derivative(s -> constant_term(p^s), 2.0) == 4log(2)
    @test ForwardDiff.derivative(s -> constant_term(s^p), 2.0) == 4
    @test ForwardDiff.derivative(s -> constant_term(log(s, p)), 2.0) ≈ -1 / (2log(2))
    @test promote_type(TaylorPolynomial, typeof(d)) === TaylorPolynomial{typeof(d)}
    @test only(
        ForwardDiff.partials(coefficient(ForwardDiff.Dual(0.0, 1.0) * x^2, [2]))
    ) == 1
    cs(a) = a isa TaylorPolynomial ? coefficients(a) : [a, zero(a), zero(a), zero(a)]
    for (a, b, c) in (
            (d, p, q), (p, d, q), (p, q, d),
            (d, p, 1 // 2), (1 // 2, d, p), (p, 1 // 2, d),
        )
        ac, bc, cc = cs(a), cs(b), cs(c)
        expected = [sum(ac[j + 1] * bc[k - j + 1] for j in 0:k) + cc[k + 1] for k in 0:3]
        result = coefficients(muladd(a, b, c))
        @test ForwardDiff.value.(result) == ForwardDiff.value.(expected)
        @test ForwardDiff.partials.(result) == ForwardDiff.partials.(expected)
    end
    map = CompiledMap([p^2])
    @test ForwardDiff.gradient(v -> only(map(v)), [0.25]) ≈ [4.5]
    @test ForwardDiff.hessian(v -> only(map(v)), [0.25]) ≈ reshape([2], 1, 1)
end

@testset "Scalar eigensolver types, symmetry and polynomial residuals" begin
    x, = variables(BigFloat, 1; order = 2)
    A = [2 + x one(x); one(x) 4 - x]
    values, vectors = DA.eigenpairs(A)
    @test coefficient_type(first(values)) === BigFloat
    @test maximum(coefficient_norm, A * vectors - vectors * Diagonal(values)) < big"1e-65"
    @test_throws ArgumentError DA.eigenpairs([one(x) x; zero(x) 2one(x)])
    @test_throws DimensionMismatch DA.eigenpairs(ones(2, 3))
    @test DA.coefficient_norm(x) isa BigFloat
end
