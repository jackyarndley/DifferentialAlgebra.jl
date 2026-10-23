using Test, DifferentialAlgebra, ForwardDiff, LinearAlgebra, GenericLinearAlgebra
using SpecialFunctions: gamma, besselj, besseli, besselix, besselk
using ..Fixtures
const DA = DifferentialAlgebra

@testset "Dual coefficients and independent scalar convolution" begin
    x, = variables(1; order = 3)
    p, q, d = 2 + x, 3 + 2x, ForwardDiff.Dual(2.0, 3.0)
    @test ForwardDiff.derivative(s -> constant_term(p^s), 2.0) == 4log(2)
    @test ForwardDiff.derivative(s -> constant_term(s^p), 2.0) == 4
    @test ForwardDiff.derivative(s -> constant_term(log(s, p)), 2.0) ≈ -1 / (2log(2))
    @test promote_type(TaylorPolynomial, typeof(d)) === TaylorPolynomial{typeof(d)}
    dual_zero = ForwardDiff.Dual(0.0, 1.0)
    nested_zero = ForwardDiff.Dual{:outer}(dual_zero, zero(dual_zero))
    for tangent in (dual_zero, nested_zero)
        jet = tangent * x^2
        @test coefficient(jet, [2]) == tangent && !iszero(jet)
        @test !isone(TaylorPolynomial(1 + tangent))
        @test_throws DA.TaylorError inv(TaylorPolynomial(tangent) + x)
        @test_throws DA.TaylorError p / tangent
        @test_throws DA.TaylorError sqrt(TaylorPolynomial(tangent))
        @test iszero(sqrt(TaylorPolynomial(zero(tangent))))
        @test iszero(hypot(TaylorPolynomial(zero(tangent)), zero(x)))
    end
    singular = TaylorPolynomial(nested_zero)
    @test_throws DA.TaylorError (singular + x)^0.5
    @test_throws DA.TaylorError hypot(singular, x)
    @test_throws DA.TaylorError asin(1 + singular)
    @test_throws DA.TaylorError gamma(singular)
    @test_throws DA.TaylorError besselk(1, singular)
    @test_throws DA.TaylorError besselix(1, singular)
    for f in (besselj, besseli)
        c = constant_term(f(1, TaylorPolynomial(dual_zero) + x))
        @test ForwardDiff.value(c) == 0 && only(ForwardDiff.partials(c)) ≈ 0.5
    end
    cs(a) = a isa TaylorPolynomial ? coefficients(a) : [a, zero(a), zero(a), zero(a)]
    for (a, b, c) in (
            (d, p, q), (p, d, q), (p, q, d),
            (d, p, 1 // 2), (1 // 2, d, p), (p, 1 // 2, d),
            (d, d, p), (d, p, d), (p, d, d),
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
    r = sqrt(big(2))
    @test hcat(coefficients.(values)...) ≈
        [3 - r 3 + r; 1 / r -1 / r; -1 / (4r) 1 / (4r)] atol = big"1e-65" rtol = 0
    # Scalar coefficient convolutions cannot hide errors through DA truncation.
    ac = ([2 1; 1 4], [1 0; 0 -1], zeros(Int, 2, 2))
    vc = [coefficient.(vectors, Ref([k])) for k in 0:2]
    lc = [Diagonal(coefficient.(values, Ref([k]))) for k in 0:2]
    for k in 0:2
        residual = sum(ac[j + 1] * vc[k - j + 1] - vc[j + 1] * lc[k - j + 1] for j in 0:k)
        norm_coefficient = sum(transpose(vc[j + 1]) * vc[k - j + 1] for j in 0:k)
        @test norm(residual) < big"1e-65"
        @test norm(norm_coefficient - (k == 0 ? I : zeros(2, 2))) < big"1e-65"
    end
    @test_throws ArgumentError DA.eigenpairs([one(x) x; zero(x) 2one(x)])
    @test_throws DimensionMismatch DA.eigenpairs(ones(2, 3))
end

@testset "Range-stable dual recurrences and final eigenpair filtering" begin
    x, = variables(1; order = 2)
    @test_throws DA.TaylorError ForwardDiff.Dual(0.0, 1.0)^(2 + x)
    d = ForwardDiff.Dual(1.0e308, 1.0)
    result = (d + d * x)^0.25
    linear = coefficient(result, [1])
    @test ForwardDiff.value(linear) ≈ 2.5e76
    @test only(ForwardDiff.partials(linear)) ≈ 6.25e-233
    nested = ForwardDiff.Dual{:outer}(d, one(d))
    q = coefficient((nested + nested * x)^0.25, [1])
    @test ForwardDiff.value(only(ForwardDiff.partials(q))) ≈ 6.25e-233
    A = [2 + 0.018x one(x); one(x) 4 - 0.018x]
    set_coefficient_tolerance!(0.01)
    values, _ = DA.eigenpairs(A)
    @test coefficient.(values, Ref([1])) ≈ [0.018 / sqrt(2), -0.018 / sqrt(2)]
end

@testset "Rational exponents stay constant with dual coefficients" begin
    x, = variables(typeof(ForwardDiff.Dual(0.0, 1.0)), 1; order = 2)
    @test coefficients(x^(2 // 1)) == coefficients(x^2)
    @test coefficient((1 + x)^(1 // 2), [2]) == ForwardDiff.Dual(-1 / 8, 0.0)
end
