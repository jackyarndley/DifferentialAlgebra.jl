using Test, DifferentialAlgebra
using ForwardDiff: ForwardDiff, Dual

@testset "Zero-valued coefficients retain nonzero Dual partials" begin
    for T in (Float32, Float64, BigFloat), budget in (0, 32 * 1024^2)
        x, = variables(T, 1; order = 3, table_bytes = budget)
        nilpotent = TaylorPolynomial{typeof(Dual(zero(T), one(T)))}(
            Dual(zero(T), one(T))
        )
        @test !iszero(nilpotent)
        @test !isone(1 + nilpotent)
        @test_throws DifferentialAlgebra.TaylorError inv(nilpotent + x)
        derivative = ForwardDiff.derivative(s -> coefficient(exp(s * x), [1]), zero(T))
        @test derivative == one(T)
    end
end

@testset "Dual powers and logarithms retain coefficient nesting" begin
    for T in (Float32, Float64, BigFloat), budget in (0, 32 * 1024^2)
        x, = variables(T, 1; order = 3, table_bytes = budget)
        d = Dual(T(2), T(1))
        @test promote_type(TaylorPolynomial, typeof(d)) === TaylorPolynomial{typeof(d)}
        @test promote_type(typeof(d), TaylorPolynomial) === TaylorPolynomial{typeof(d)}
        derivative = ForwardDiff.derivative(s -> constant_term((1 + x)^s), T(2))
        @test derivative isa T
        @test iszero(derivative)
        second(s) = coefficient((1 + x)^s, [2])
        @test ForwardDiff.derivative(s -> ForwardDiff.derivative(second, s), T(2)) ≈ one(T)

        result = (1 + x)^d
        @test result isa TaylorPolynomial{typeof(d)}
        if result isa TaylorPolynomial
            for (k, value, partial) in ((0, 1, 0), (1, 2, 1), (2, 1, 3 // 2), (3, 0, 1 // 3))
                c = coefficient(result, [k])
                @test ForwardDiff.value(c) ≈ T(value) atol = 16eps(T) rtol = 32eps(T)
                @test only(ForwardDiff.partials(c)) ≈ T(partial) atol = 16eps(T) rtol = 32eps(T)
            end
        end

        power = d^(1 + x)
        logarithm = log(d, 4 + x)
        @test power isa TaylorPolynomial{typeof(d)}
        @test logarithm isa TaylorPolynomial{typeof(d)}
        logtwo = log(T(2))
        for k in 0:3
            if power isa TaylorPolynomial
                c = coefficient(power, [k])
                value = 2logtwo^k / factorial(k)
                partial = logtwo^k / factorial(k) +
                    (k == 0 ? zero(T) : logtwo^(k - 1) / factorial(k - 1))
                @test ForwardDiff.value(c) ≈ value
                @test only(ForwardDiff.partials(c)) ≈ partial
            end
            if logarithm isa TaylorPolynomial
                c = coefficient(logarithm, [k])
                value = k == 0 ? log(T(4)) : T((-1)^(k - 1)) / (k * T(4)^k)
                @test ForwardDiff.value(c) ≈ value / logtwo
                @test only(ForwardDiff.partials(c)) ≈ -value / (2logtwo^2)
            end
        end
    end
end

@testset "Dual muladd keeps a polynomial of scalar Dual coefficients" begin
    x, = variables(1; order = 3)
    p, q, d = 4 + x, 1 - x, Dual(2.0, 3.0)
    coefficients(a) = a isa TaylorPolynomial ? [coefficient(a, [k]) for k in 0:3] :
        [a, zero(a), zero(a), zero(a)]
    cases = Tuple[(d, p, q), (p, d, q), (p, q, d), (d, d, p), (d, p, d), (p, d, d)]
    for scalar in (1, 1.0, 1 // 2, ℯ, BigFloat(1))
        append!(
            cases, (
                (d, p, scalar), (p, d, scalar), (d, scalar, p),
                (scalar, d, p), (p, scalar, d), (scalar, p, d),
            )
        )
    end
    for (a, b, c) in cases
        ac, bc, cc = coefficients(a), coefficients(b), coefficients(c)
        # Independent scalar convolution, without polynomial multiplication.
        expected = [
            sum(ac[j + 1] * bc[k - j + 1] for j in 0:k) + cc[k + 1] for k in 0:3
        ]
        result = muladd(a, b, c)
        @test result isa TaylorPolynomial{typeof(first(expected))}
        if result isa TaylorPolynomial
            @test all(0:3) do k
                actual, reference = coefficient(result, [k]), expected[k + 1]
                isequal(ForwardDiff.value(actual), ForwardDiff.value(reference)) &&
                    isequal(ForwardDiff.partials(actual), ForwardDiff.partials(reference))
            end
        end
    end
    @test ForwardDiff.derivative(s -> constant_term(muladd(s, p, q)), 2.0) == 4
    @test ForwardDiff.derivative(2.0) do s
        ForwardDiff.derivative(t -> constant_term(muladd(s, t, p)), 3.0)
    end == 1
end
