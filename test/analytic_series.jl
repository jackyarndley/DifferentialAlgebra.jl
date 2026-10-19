using Test, DifferentialAlgebra
using ForwardDiff, SpecialFunctions
include("support/series.jl")

@testset verbose = true "TEST 2: Single-variable functions" begin
    k = 10  # order
    eps = 1.0e-15  # tolerance for comparing real numbers

    @testset "2.1 Power function" begin
        for n in 1:10
            @testset "2.1 Power function (n=$(n))" begin
                DifferentialAlgebra.initialize!(k, 1)

                x = TaylorPolynomial(1, 1)
                f = (1 + x)^n

                jj = Vector{UInt32}(undef, 1)
                for i in 0:k
                    jj[1] = i
                    a_da = DifferentialAlgebra.coefficient(f, jj)
                    if i <= n
                        a_exact = binomial(n, i)
                    else
                        a_exact = 0.0
                    end
                    @test isapprox(a_da, a_exact, atol = eps)
                end
            end
        end
    end

    @testset "2.2 Division" begin
        DifferentialAlgebra.initialize!(k, 1)

        x = TaylorPolynomial(1, 1)
        f = 1 / (1 - x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.coefficient(f, jj)
            a_exact = 1.0
            @test isapprox(a_da, a_exact, atol = eps)
        end
    end

    @testset "2.3 nth root function" begin
        for n in 2:5
            @testset "2.3 nth root function (n=$(n))" begin
                DifferentialAlgebra.initialize!(k, 1)

                x = TaylorPolynomial(1, 1)
                f = DifferentialAlgebra.nthroot(1 + x, n)

                jj = Vector{UInt32}(undef, 1)
                args = Vector{Float64}(undef, 0)
                for i in 0:k
                    jj[1] = i
                    a_da = DifferentialAlgebra.coefficient(f, jj)
                    if i == 0
                        a_exact = 1.0
                    else
                        push!(args, 1.0 / n - (i - 1))
                        a_exact = 1.0 / factorial(i) * prod(args)
                    end
                    @test isapprox(a_da, a_exact, atol = eps)
                end
            end
        end
    end

    @testset "2.4 Exponential function" begin
        DifferentialAlgebra.initialize!(k, 1)

        x = TaylorPolynomial(1, 1)
        f = exp(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.coefficient(f, jj)
            a_exact = 1 / factorial(i)
            @test isapprox(a_da, a_exact, atol = eps)
        end
    end

    @testset "2.5 Logarithmic function" begin
        DifferentialAlgebra.initialize!(k, 1)

        x = TaylorPolynomial(1, 1)
        f = log(1 + x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.coefficient(f, jj)
            if i == 0
                a_exact = 0.0
            else
                a_exact = (-1)^(i + 1) * 1.0 / i
            end
            @test isapprox(a_da, a_exact, atol = eps)
        end
    end

    @testset "2.6 Sine function" begin
        DifferentialAlgebra.initialize!(k, 1)

        x = TaylorPolynomial(1, 1)
        f = sin(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.coefficient(f, jj)
            if iseven(i)
                a_exact = 0.0
            else
                a_exact = (-1)^((i - 1) / 2) / factorial(i)
            end
            @test isapprox(a_da, a_exact, atol = eps)
        end
    end

    @testset "2.7 Cosine function" begin
        DifferentialAlgebra.initialize!(k, 1)

        x = TaylorPolynomial(1, 1)
        f = cos(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.coefficient(f, jj)
            if iseven(i)
                a_exact = (-1)^(i / 2) / factorial(i)
            else
                a_exact = 0.0
            end
            @test isapprox(a_da, a_exact, atol = eps)
        end
    end

    @testset "2.8 Tangent function" begin
        DifferentialAlgebra.initialize!(k, 1)

        x = TaylorPolynomial(1, 1)
        f = tan(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.coefficient(f, jj)
            if iseven(i)
                a_exact = 0.0
            else
                a_exact = (-1)^((i - 1) / 2) * 2^(i + 1) * (2^(i + 1) - 1) * bernoulli(i + 1) / factorial(i + 1)
            end
            @test isapprox(a_da, a_exact, atol = eps)
        end
    end

    @testset "2.9 Arcsine function" begin
        DifferentialAlgebra.initialize!(k, 1)

        x = TaylorPolynomial(1, 1)
        f = asin(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.coefficient(f, jj)
            if iseven(i)
                a_exact = 0.0
            else
                a_exact = factorial(i - 1) / (4.0^((i - 1) / 2) * factorial(div(i - 1, 2))^2 * i)
            end
            @test isapprox(a_da, a_exact, atol = eps)
        end
    end

    @testset "2.10 Arccosine function" begin
        DifferentialAlgebra.initialize!(k, 1)

        x = TaylorPolynomial(1, 1)
        f = acos(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.coefficient(f, jj)
            if i == 0
                a_exact = pi / 2.0
            elseif iseven(i)
                a_exact = 0.0
            else
                a_exact = - factorial(i - 1) / (4.0^((i - 1) / 2) * factorial(div(i - 1, 2))^2 * i)
            end
            @test isapprox(a_da, a_exact, atol = eps)
        end
    end

    @testset "2.11 Arctangent function" begin
        DifferentialAlgebra.initialize!(k, 1)

        x = TaylorPolynomial(1, 1)
        f = atan(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.coefficient(f, jj)
            if iseven(i)
                a_exact = 0.0
            else
                a_exact = (-1)^((i - 1) / 2) / i
            end
            @test isapprox(a_da, a_exact, atol = eps)
        end
    end

    @testset "2.12 Hyperbolic sine function" begin
        DifferentialAlgebra.initialize!(k, 1)

        x = TaylorPolynomial(1, 1)
        f = sinh(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.coefficient(f, jj)
            if iseven(i)
                a_exact = 0.0
            else
                a_exact = 1.0 / factorial(i)
            end
            @test isapprox(a_da, a_exact, atol = eps)
        end
    end

    @testset "2.13 Hyperbolic cosine function" begin
        DifferentialAlgebra.initialize!(k, 1)

        x = TaylorPolynomial(1, 1)
        f = cosh(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.coefficient(f, jj)
            if iseven(i)
                a_exact = 1.0 / factorial(i)
            else
                a_exact = 0.0
            end
            @test isapprox(a_da, a_exact, atol = eps)
        end
    end

    @testset "2.14 Hyperbolic tangent function" begin
        DifferentialAlgebra.initialize!(k, 1)

        x = TaylorPolynomial(1, 1)
        f = tanh(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.coefficient(f, jj)
            if iseven(i)
                a_exact = 0.0
                @test isapprox(a_da, a_exact, atol = eps)
            else
                # tanh(x) = sum(2^(2n)*(2^(2n)-1)*B_(2n)/(2n)! * x^(2n-1)).
                a_exact = 2^(i + 1) * (2^(i + 1) - 1) * bernoulli(i + 1) / factorial(i + 1)
                @test isapprox(a_da, a_exact, atol = eps)
            end
        end
    end

    @testset "2.15 Hyperbolic arcsine function" begin
        DifferentialAlgebra.initialize!(k, 1)

        x = TaylorPolynomial(1, 1)
        f = asinh(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.coefficient(f, jj)
            if iseven(i)
                a_exact = 0.0
            else
                a_exact = (-1)^((i - 1) / 2) * factorial(i - 1) / (4^((i - 1) / 2) * factorial(div(i - 1, 2))^2 * i)
            end
            @test isapprox(a_da, a_exact, atol = eps)
        end
    end

    @testset "2.16 Hyperbolic arctangent function" begin
        DifferentialAlgebra.initialize!(k, 1)

        x = TaylorPolynomial(1, 1)
        f = atanh(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.coefficient(f, jj)
            if iseven(i)
                a_exact = 0.0
            else
                a_exact = 1 / i
            end
            # Compare the analytic 1/i coefficients with an absolute tolerance.
            @test isapprox(a_da, a_exact, atol = 1.0e-12)
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
