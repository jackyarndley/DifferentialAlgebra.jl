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
                @test [coefficient(p, [k]) for k in 0:order] ≈
                    Float64[oracle(k) for k in 0:order] atol = 1.0e-14 rtol = 1.0e-14
            end
        end
        for n in (2, 5, 10)
            p = (1 + x)^n
            @test [coefficient(p, [k]) for k in 0:order] ≈ [binomial(n, k) for k in 0:order]
        end
        for n in 2:5
            p = DifferentialAlgebra.nthroot(1 + x, n)
            expected = [prod(1 // n - j for j in 0:(k - 1); init = 1) / factorial(big(k)) for k in 0:order]
            @test [coefficient(p, [k]) for k in 0:order] ≈ Float64.(expected) atol = 1.0e-15
        end
    end
end

# Scalar automatic differentiation supplies f^(k)(c), independently of the DA
# coefficient recurrences and composition engine. These are numerical regression
# expectations, not uniform interval inclusion oracles.
function scalar_derivative(f, x, n)
    return n == 0 ? f(x) : ForwardDiff.derivative(t -> scalar_derivative(f, t, n - 1), x)
end

# Expand f(2+h), h = a*x + b*y + c*x*z + d*y^2, directly by the multinomial
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
            p = 2 + 0.1x + 0.2y + 0.05x * z + 0.03y^2
            for f in (sqrt, inv, exp, log, sin, cos, tan, atan, tanh, gamma, loggamma, erf)
                @testset "$(f), table_bytes=$(budget)" begin
                    series = [scalar_derivative(f, 2.0, k) / factorial(k) for k in 0:6]
                    result = f(p)
                    for i in 0:6, j in 0:(6 - i), k in 0:(6 - i - j)
                        alpha = (i, j, k)
                        expected = multinomial_composition_coefficient(series, weights, alpha)
                        @test coefficient(result, collect(alpha)) ≈ expected atol = 1.0e-12 rtol = 0
                    end
                end
            end
        end
    end
end
