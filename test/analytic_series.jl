using Test, DifferentialAlgebra
using ForwardDiff, SpecialFunctions

include("support/series.jl")

@testset "Filter completed analytic results" begin
    for budget in (0, 32 * 1024^2)
        x, = variables(1; order = 3, table_bytes = budget)
        set_coefficient_tolerance!(1.0e-6)
        q = 1.0e8 / (1.0e8 + 1.0e6x)
        @test constant_term(q) == 1
        @test coefficient(q, [1]) ≈ -0.01
        @test iszero(inv(1.0e8 + 1.0e6x)) # Public reciprocal obeys tolerance.
        @test coefficient_tolerance() == 1.0e-6

        x, = variables(1; order = 10, table_bytes = budget)
        set_coefficient_tolerance!(1.0e-4)
        # Scalar differentiation gives exp(-5)*100^k/k!, independently of DA.
        p = exp(-5 + 100x)
        @test coefficient(p, [10]) ≈ exp(-5) * 100.0^10 / factorial(10)
        @test coefficient_tolerance() == 1.0e-4

        x, = variables(1; order = 2, table_bytes = budget)
        set_coefficient_tolerance!(1.0e-4)
        # erf''(4)*100^2/2 = -8exp(-16)*10000/sqrt(pi).
        p = erf(4 + 100x)
        @test coefficient(p, [2]) ≈ -0.005079293868746624
        @test coefficient_tolerance() == 1.0e-4

        x, = variables(1; order = 2, table_bytes = budget)
        set_coefficient_tolerance!(1.0e-4)
        # A normalized slope below tolerance can yield a retained coefficient
        # after multiplication by the square root of a large expansion point.
        @test coefficient(sqrt(1.0e16 + 1.0e11x), [1]) ≈ 500
        @test coefficient(log(1.0e16 + 1.0e11x), [1]) == 0
    end
end

@testset "Promoted analytic intermediates" begin
    for budget in (0, 32 * 1024^2)
        for T in (
                UInt8, UInt64, UInt128, Rational{UInt8}, Rational{UInt64},
                Rational{UInt128},
            )
            x, = variables(T, 1; order = 4, table_bytes = budget)
            denominator = T(3) * one(x) + x
            for result in (inv(denominator), one(T) / denominator)
                for k in 0:4
                    # The geometric series is independent of either DA kernel.
                    expected = (-big(1))^k // big(3)^(k + 1)
                    @test coefficient(result, [k]) ≈ expected
                end
            end
            @test constant_term(-1 / TaylorPolynomial{T}(T(3))) ≈ -1 // 3
            @test coefficient(denominator, [1]) == 1
        end
        setprecision(256) do
            x, = variables(1; order = 6, table_bytes = budget)
            denominator = 3 + x
            exponent = BigFloat(3) / 7
            for (result, power, scale) in (
                    (BigFloat(2) / denominator, -one(BigFloat), BigFloat(2)),
                    (denominator^exponent, exponent, one(BigFloat)),
                )
                @test coefficient_type(result) === BigFloat
                expected = setprecision(512) do
                    p = BigFloat(power)
                    [
                        BigFloat(scale) * big(3)^p *
                            prod(p - j for j in 0:(k - 1); init = one(p)) /
                            (big(3)^k * factorial(big(k))) for k in 0:6
                    ]
                end
                for k in 0:6
                    @test coefficient(result, [k]) ≈ expected[k + 1] rtol = 64eps(BigFloat)
                end
            end
        end
    end
end

@testset "Univariate analytic coefficients" begin
    order = 10
    # Exact scalar formulas supply the oracle, including Bernoulli numbers for
    # tangent and rational binomial coefficients for fractional powers.
    arcsine(k) = iseven(k) ? 0 : factorial(big(k - 1)) //
        (big(4)^((k - 1) ÷ 2) * factorial(big((k - 1) ÷ 2))^2 * k)
    tangent(k) = iseven(k) ? 0 : (-1)^((k - 1) ÷ 2) * big(2)^(k + 1) *
        (big(2)^(k + 1) - 1) * bernoulli(k + 1) / factorial(big(k + 1))
    cases = (
        (inv, k -> (-1)^k),
        (exp, k -> 1 // factorial(big(k))),
        (log, k -> k == 0 ? 0 : (-1)^(k + 1) // k),
        (sin, k -> iseven(k) ? 0 : (-1)^((k - 1) ÷ 2) // factorial(big(k))),
        (cos, k -> isodd(k) ? 0 : (-1)^(k ÷ 2) // factorial(big(k))),
        (tan, tangent), (asin, arcsine),
        (acos, k -> k == 0 ? pi / 2 : -arcsine(k)),
        (atan, k -> iseven(k) ? 0 : (-1)^((k - 1) ÷ 2) // k),
        (sinh, k -> iseven(k) ? 0 : 1 // factorial(big(k))),
        (cosh, k -> isodd(k) ? 0 : 1 // factorial(big(k))),
        (tanh, k -> iseven(k) ? 0 : (-1)^((k - 1) ÷ 2) * tangent(k)),
        (asinh, k -> iseven(k) ? 0 : (-1)^((k - 1) ÷ 2) * arcsine(k)),
        (atanh, k -> iseven(k) ? 0 : 1 // k),
    )
    for budget in (0, 32 * 1024^2)
        x, = variables(1; order, table_bytes = budget)
        for (f, oracle) in cases
            @testset "$f, table_bytes=$budget" begin
                p = f in (inv, log) ? f(1 + x) : f(x)
                actual = [coefficient(p, [k]) for k in 0:order]
                expected = Float64[oracle(k) for k in 0:order]
                @test all(isapprox.(actual, expected; atol = 1.0e-15, rtol = 0))
            end
        end
        for n in (2, 5, 10)
            p = (1 + x)^n
            @test [coefficient(p, [k]) for k in 0:order] == [binomial(n, k) for k in 0:order]
        end
        for n in 2:5
            p = DifferentialAlgebra.nthroot(1 + x, n)
            expected = [prod(1 // n - j for j in 0:(k - 1); init = 1) / factorial(big(k)) for k in 0:order]
            @test all(isapprox.([coefficient(p, [k]) for k in 0:order], Float64.(expected); atol = 1.0e-15, rtol = 0))
        end
    end
end

# Scalar automatic differentiation supplies f^(k)(c), independently of the DA
# coefficient recurrences and composition engine. These are numerical regression
# expectations, not uniform interval inclusion oracles.
function scalar_derivative(f, x, n)
    return n == 0 ? f(x) : ForwardDiff.derivative(t -> scalar_derivative(f, t, n - 1), x)
end

@testset "Filter completed compound results" begin
    cases = (
        t -> log2(2 + 0.018t), t -> log10(2 + 0.18t),
        t -> log(1.1, 2 + 0.018t), t -> (10 + 0.09t)^3,
        t -> (10 + 0.09t)^5, t -> (0.5 + 0.09t)^(-3),
        t -> hypot(10 + 0.09t, 2 + 0.018t),
        t -> gamma(10 + 0.09t), t -> gamma(-2.5 + 0.09t),
        t -> (10 + 0.09t)^(2 + 0.018t), t -> 1.1^(100 + 0.018t),
        t -> atan(2 + 0.018t, 10 + 0.09t),
        t -> atan(10 + 0.09t, 2 + 0.018t),
        t -> muladd(10 + 0.09t, 10 + 0.09t, 0.018t^2),
    )
    scalar_cases = (
        t -> hypot(1 + 100t, 0.005), t -> hypot(0.005, 1 + 100t),
        t -> atan(1 + 100t, 0.005), t -> atan(0.005, 1 + 100t),
    )
    for budget in (0, 32 * 1024^2)
        x, = variables(1; order = 4, table_bytes = budget)
        stored = 10 + 0.009x + x^4
        set_coefficient_tolerance!(0.01)
        # BigFloat derivatives avoid cancellation in the small scalar-input
        # corrections to hypot; the polynomial calculation still uses Float64.
        for (functions, point) in ((cases, 0.0), (scalar_cases, BigFloat(0)))
            for f in functions
                result = f(x)
                for k in 0:4
                    expected = scalar_derivative(f, point, k) / factorial(k)
                    expected = abs(expected) <= 0.01 ? 0.0 : expected
                    @test coefficient(result, [k]) ≈ expected rtol = 1.0e-10 atol = 1.0e-12
                end
                @test coefficient_tolerance() == 0.01
                @test truncation_order() == 4
            end
        end
        @test coefficient(log2(2 + 0.018x), [1]) ≈ 0.01298425536800067
        @test coefficient((10 + 0.09x)^3, [2]) ≈ 0.243
        # The public result is owned, including identity and negative powers.
        p = 10 + 0.09x
        snapshot = copy(p.coeffs)
        for n in (0, 1, 2, 3, -1, -3)
            result = p^n
            result.coeffs[1] = 7
            @test p.coeffs == snapshot
        end
        with_order(2) do
            @test degree(p^5) <= 2
            @test coefficient(stored^3, [1]) ≈ 2.7
            @test truncation_order() == 2
        end
        @test truncation_order() == 4
    end
end

@testset "Compound coefficient precision" begin
    for bits in (128, 256)
        setprecision(BigFloat, bits) do
            x, = variables(BigFloat, 1; order = 3)
            slope = BigFloat("0.018")
            p = 2 + slope * x
            set_coefficient_tolerance!(0.01)
            for (result, expected) in (
                    (log2(p), slope / (2log(big(2)))),
                    (log(1 // 2, p), slope / (2log(big(1) // 2))),
                    (3^p, 9log(big(3)) * slope),
                )
                @test coefficient_type(result) === BigFloat
                actual = coefficient(result, [1])
                @test abs(expected) <= 0.01 ? iszero(actual) :
                    isapprox(actual, expected; rtol = 20eps(BigFloat))
            end
        end
    end
    x, = variables(1; order = 2)
    @test constant_term(asin(one(x))) == pi / 2
    @test constant_term(acos(one(x))) == 0
end

# Expand f(c+h), h = a*x + b*y + c*x*z + d*y^2, directly by the multinomial
# theorem. Counts of x*z and x are fixed by the target x/z exponents; enumerate
# the possible y^2 counts. No polynomial multiplication or basis metadata is used.
function multinomial_composition_coefficient(series, weights, alpha)
    nxz = alpha[3]
    nx = alpha[1] - nxz
    nx < 0 && return zero(BigFloat)
    value = zero(BigFloat)
    for nyy in 0:(alpha[2] ÷ 2)
        ny = alpha[2] - 2nyy
        counts = (nx, ny, nxz, nyy)
        k = sum(counts)
        k < length(series) || continue
        multiplicity = factorial(big(k)) ÷ prod(n -> factorial(big(n)), counts)
        weight = prod(w^n for (w, n) in zip(weights, counts))
        value += BigFloat(series[k + 1]) * BigFloat(multiplicity * weight)
    end
    return value
end

@testset "Multivariate coefficients from scalar derivatives" begin
    setprecision(256) do
        @testset "Combinatorial oracle identities" begin
            series, weights = (1, 2, 3), (2 // 1, 3 // 1, 5 // 1, 7 // 1)
            @test multinomial_composition_coefficient(series, weights, (0, 0, 0)) == 1
            @test multinomial_composition_coefficient(series, weights, (0, 0, 1)) == 0
            @test multinomial_composition_coefficient(series, weights, (1, 0, 1)) == 10
            @test multinomial_composition_coefficient(series, weights, (0, 2, 0)) == 41
            @test multinomial_composition_coefficient(series, weights, (1, 1, 1)) == 90
            @test multinomial_composition_coefficient(series, weights, (2, 0, 2)) == 75
        end
        # Float64 inputs denote their stored binary values. Make those values
        # exact rationals for the combinatorial part of the independent oracle.
        weights = Rational{BigInt}.((0.1, 0.2, 0.05, 0.03))
        for budget in (0, 32 * 1024^2)
            x, y, z = variables(3; order = 6, table_bytes = budget)
            cases = (
                (sqrt, 2.0), (inv, 2.0), (exp, 2.0), (log, 2.0),
                (sin, 2.0), (cos, 2.0), (tan, 2.0), (atan, 2.0),
                (tanh, 2.0), (gamma, 2.0), (loggamma, 2.0), (erf, 2.0),
                (asin, 0.25), (acos, 0.25), (asinh, 0.25), (atanh, 0.25),
                (acosh, 2.0), (sinh, 0.25), (cosh, 0.25), (cbrt, 2.0),
                (log2, 2.0), (log10, 2.0), (erfc, 0.25), (digamma, 2.0),
                (Base.Fix1(besselj, 2), 2.0), (Base.Fix1(bessely, 2), 2.0),
                (Base.Fix1(besseli, 2), 2.0), (Base.Fix1(besselk, 2), 2.0),
                (Base.Fix1(besselix, 2), 2.0), (Base.Fix1(besselkx, 2), 2.0),
            )
            argument = 2 + 0.1x + 0.2y + 0.05x * z + 0.03y^2
            pair = sincos(argument)
            @test all(isequal.(pair, (sin(argument), cos(argument))))
            @test pair[1].coeffs !== pair[2].coeffs
            @test all(p -> p.coeffs !== argument.coeffs, pair)
            for (f, center) in cases
                @testset "$(f) at $center, table_bytes=$budget" begin
                    p = center + 0.1x + 0.2y + 0.05x * z + 0.03y^2
                    series = [scalar_derivative(f, center, k) / factorial(k) for k in 0:6]
                    result = f(p)
                    for i in 0:6, j in 0:(6 - i), k in 0:(6 - i - j)
                        alpha = (i, j, k)
                        expected = multinomial_composition_coefficient(series, weights, alpha)
                        @test coefficient(result, collect(alpha)) ≈ expected atol = 1.0e-12 rtol = 0
                        if f === sin || f === cos
                            member = f === sin ? pair[1] : pair[2]
                            actual = coefficient(member, collect(alpha))
                            @test actual ≈ expected atol = 1.0e-12 rtol = 0
                        end
                    end
                end
            end
        end
    end
end

@testset "Trigonometric kernels at the working order" begin
    references = setprecision(512) do
        s, c = sincos(big(1) / 2)
        sh, ch = sinh(big(1) / 2), cosh(big(1) / 2)
        map(
            (
                (sin, (s, c, -s, -c)), (cos, (c, -s, -c, s)),
                (sinh, (sh, ch, sh, ch)), (cosh, (ch, sh, ch, sh)),
            )
        ) do (f, cycle)
            # Compose the scalar Taylor series with 2x + 3x^5 explicitly.
            values = [
                cycle[mod(k, 4) + 1] * big(2)^k / factorial(big(k)) +
                    (k == 5 ? 3cycle[2] : 0) for k in 0:5
            ]
            (f, values)
        end
    end
    for T in (Float32, Float64, BigFloat), budget in (0, 32 * 1024^2)
        x, = variables(T, 1; order = 5, table_bytes = budget)
        argument = T(1 // 2) + 2x + 3x^5
        constant = T(1 // 2) * one(x)
        for order in (1, 2, 5)
            with_order(order) do
                for (f, reference) in references
                    value = @inferred f(argument)
                    @test value isa TaylorPolynomial{T}
                    @test degree(value) <= order
                    @test value.coeffs !== argument.coeffs
                    for k in 0:order
                        @test coefficient(value, [k]) ≈ T(reference[k + 1]) rtol = 64eps(T)
                    end
                    scalar = @inferred f(constant)
                    @test constant_term(scalar) == f(T(1 // 2))
                    @test scalar.coeffs !== constant.coeffs
                end
                pair = @inferred sincos(argument)
                for (value, (_, reference)) in zip(pair, references[1:2])
                    for k in 0:order
                        @test coefficient(value, [k]) ≈ T(reference[k + 1]) rtol = 64eps(T)
                    end
                end
                @test pair[1].coeffs !== pair[2].coeffs
                cosine = copy(pair[2].coeffs)
                set_coefficient!(pair[1], [1], T(9))
                @test pair[2].coeffs == cosine
                @test coefficient(argument, [1]) == T(2)
                @test coefficient(argument, [5]) == T(3)
            end
        end
    end
end

@testset "Filter completed trigonometric results" begin
    for budget in (0, 32 * 1024^2)
        x, y = variables(2; order = 3, table_bytes = budget)
        argument = 0.001x + 100y
        tiny_constant = TaylorPolynomial(0.005)
        set_coefficient_tolerance!(0.01)
        s, c = sincos(argument)
        # The mixed coefficients grow from an input slope below tolerance.
        for value in (s, sin(argument))
            @test coefficient(value, [1, 0]) == 0
            @test coefficient(value, [1, 2]) ≈ -5
        end
        for value in (c, cos(argument))
            @test coefficient(value, [1, 1]) ≈ -0.1
        end
        @test coefficient(argument, [1, 0]) == 0.001
        @test iszero(sin(tiny_constant))
        @test iszero(first(sincos(tiny_constant)))
        @test constant_term(tiny_constant) == 0.005
        @test coefficient_tolerance() == 0.01
        @test truncation_order() == 3
        set_coefficient_tolerance!(0)
    end
end

@testset "Inverse functions at boundary and large centers" begin
    for budget in (0, 32 * 1024^2)
        x, = variables(1; order = 4, table_bytes = budget)
        for (f, center) in (
                (asin, 1 - 1.0e-8), (acos, -1 + 1.0e-8),
                (atanh, 1 - 1.0e-8), (atanh, -1 + 1.0e-8), (acosh, 1 + 1.0e-8),
            )
            p = f(center + x)
            expected = setprecision(BigFloat, 256) do
                [scalar_derivative(f, BigFloat(center), k) / factorial(k) for k in 0:4]
            end
            for k in 0:4
                @test coefficient(p, [k]) ≈ Float64(expected[k + 1]) rtol = 5.0e-13
            end
        end
        for f in (asinh, acosh)
            p = f(1.0e200 + x)
            @test coefficient(p, [1]) ≈ 1.0e-200
            p = f(1.0e200 + 1.0e200x)
            for k in 1:4
                @test coefficient(p, [k]) ≈ (-1)^(k - 1) / k rtol = 5.0e-14
            end
        end
        p = asinh(-1.0e200 + 1.0e200x)
        for k in 1:4
            @test coefficient(p, [k]) ≈ 1 / k rtol = 5.0e-14
        end
    end
end

@testset "Hypotenuse at small and large centers" begin
    for budget in (0, 32 * 1024^2)
        x, = variables(1; order = 4, table_bytes = budget)
        for scale in (1.0e-170, 1.0e170), sign in (-1, 1)
            p = sign * scale + x
            for result in (hypot(p, 0), hypot(0, p))
                @test constant_term(result) == scale
                @test coefficient(result, [1]) == sign
                @test all(k -> iszero(coefficient(result, [k])), 2:4)
                @test isfinite(result)
            end
            # A uniformly scaled 3-4-5 line supplies an exact independent oracle.
            p, q = 3scale * (1 + x), 4scale * (1 + x)
            result = hypot(p, q)
            @test constant_term(result) ≈ 5scale
            @test coefficient(result, [1]) ≈ 5scale
            @test all(k -> abs(coefficient(result, [k])) <= 100eps(scale), 2:4)
        end
        @test iszero(hypot(zero(x), zero(x)))
        @test_throws DifferentialAlgebra.TaylorError hypot(x, zero(x))
        @test_throws DifferentialAlgebra.TaylorError hypot(x, x)
        x, = variables(1; order = 2, table_bytes = budget)
        scale = 1.0e-170
        for result in (hypot(scale + x, scale + x / 2), hypot(scale + x / 2, scale + x))
            @test isfinite(result)
            @test constant_term(result) ≈ sqrt(2) * scale
            @test coefficient(result, [1]) ≈ 3 / (2sqrt(2)) rtol = 5.0e-14
            # h''(0)/2 = (a'b-ab')²/(2hypot(a,b)³), evaluated without tiny squares.
            @test coefficient(result, [2]) ≈ inv(16sqrt(2) * scale) rtol = 5.0e-14
        end
    end
end

@testset "NaN constants propagate through analytic function families" begin
    for T in (Float32, Float64, BigFloat)
        variables(T, 1; order = 3)
        p = TaylorPolynomial{T}(T(NaN))
        for fn in (
                exp, log, sqrt, sin, cos, gamma, loggamma, digamma,
                x -> besselj(0, x), x -> besselix(2, x),
            )
            result = fn(p)
            @test coefficient_type(result) === T
            @test isnan(constant_term(result))
        end
    end
end

@testset "Real bases and analytic domains" begin
    for budget in (0, 32 * 1024^2), T in (Float32, Float64, BigFloat)
        x, = variables(T, 1; order = 3, table_bytes = budget)
        @test iszero(T(0)^(2 + x))
        @test isone(T(0)^zero(x))
        @test constant_term(T(-2)^(2 + zero(x))) == 4
        for exponent in (x, -x, -2 + x, -one(x))
            @test_throws DifferentialAlgebra.TaylorError T(0)^exponent
        end
        @test_throws DifferentialAlgebra.TaylorError T(-2)^x
        @test_throws DifferentialAlgebra.TaylorError (-one(x) + x)^T(0.5)
        @test_throws ArgumentError DifferentialAlgebra.nthroot(x, 0)
        @test_throws ArgumentError DifferentialAlgebra.nthroot(T(1), 0)
    end
end

@testset "Large integer exponents" begin
    for budget in (0, 32 * 1024^2)
        x, = variables(1; order = 2, table_bytes = budget)
        for n in (big(2)^31, big(10)^100, typemax(UInt))
            @test isone(one(x)^n)
            @test iszero(x^n)
            p = (1 + x)^n
            @test constant_term(p) == 1
            @test coefficient(p, [1]) ≈ Float64(n)
            @test coefficient(p, [2]) ≈ Float64(binomial(big(n), 2))
        end
        @test isone(one(x)^typemin(Int))
        @test isone((-one(x))^typemin(Int))
        @test_throws DifferentialAlgebra.TaylorError zero(x)^typemin(Int)
        @test_throws DifferentialAlgebra.TaylorError x^typemin(Int)
        p = 1 + x
        q = p^big(2)^31
        set_coefficient!(q, [1], 7)
        @test coefficient(p, [1]) == 1
    end
end

@testset "Integer storage with a constant negative exponent" begin
    x, = variables(Int, 1; order = 3)
    @test constant_term(2^(-one(x))) == 0.5
    @test constant_term((-2)^(-one(x))) == -0.5
    @test coefficient_type(2^one(x)) === Float64
    for base in (2, -2)
        exponent = TimeSeries([-1, 0, 0])
        @test coefficient(base^exponent, 0) == inv(base)
        @test coefficient_type(base^exponent) === Float64
    end
end
