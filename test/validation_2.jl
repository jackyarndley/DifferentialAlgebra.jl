@testset verbose = true "TEST 2: Single-variable functions" begin
    k = 10  # order
    eps = 1e-15  # tolerance for comparing real numbers

    @testset "2.1 Power function" begin
        for n = 1:10
            @testset "2.1 Power function (n=$(n))" begin
                DifferentialAlgebra.init(k, 1)

                x = DA(1, 1)
                f = (1 + x)^n

                jj = Vector{UInt32}(undef, 1)
                for i in 0:k
                    jj[1] = i
                    a_da = DifferentialAlgebra.getCoefficient(f, jj)
                    if i <= n
                        a_exact = binomial(n, i)
                    else
                        a_exact = 0.0
                    end
                    @test isapprox(a_da, a_exact, atol=eps)
                end
            end
        end
    end

    @testset "2.2 Division" begin
        DifferentialAlgebra.init(k, 1)

        x = DA(1, 1)
        f = 1 / (1 - x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.getCoefficient(f, jj)
            a_exact = 1.0
            @test isapprox(a_da, a_exact, atol=eps)
        end
    end

    @testset "2.3 nth root function" begin
        for n = 2:5
            @testset "2.3 nth root function (n=$(n))" begin
                DifferentialAlgebra.init(k, 1)

                x = DA(1, 1)
                f = DifferentialAlgebra.root(1 + x, n)

                jj = Vector{UInt32}(undef, 1)
                args = Vector{Float64}(undef, 0)
                for i in 0:k
                    jj[1] = i
                    a_da = DifferentialAlgebra.getCoefficient(f, jj)
                    if i == 0
                        a_exact = 1.0
                    else
                        push!(args, 1.0 / n - (i - 1))
                        a_exact = 1.0 / factorial(i) * prod(args)
                    end
                    @test isapprox(a_da, a_exact, atol=eps)
                end
            end
        end
    end

    @testset "2.4 Exponential function" begin
        DifferentialAlgebra.init(k, 1)

        x = DA(1, 1)
        f = exp(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.getCoefficient(f, jj)
            a_exact = 1 / factorial(i)
            @test isapprox(a_da, a_exact, atol=eps)
        end
    end

    @testset "2.5 Logarithmic function" begin
        DifferentialAlgebra.init(k, 1)

        x = DA(1, 1)
        f = log(1 + x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.getCoefficient(f, jj)
            if i == 0
                a_exact = 0.0
            else
                a_exact = (-1)^(i+1) * 1.0 / i
            end
            @test isapprox(a_da, a_exact, atol=eps)
        end
    end

    @testset "2.6 Sine function" begin
        DifferentialAlgebra.init(k, 1)

        x = DA(1, 1)
        f = sin(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.getCoefficient(f, jj)
            if iseven(i)
                a_exact = 0.0
            else
                a_exact = (-1)^((i - 1) / 2) / factorial(i)
            end
            @test isapprox(a_da, a_exact, atol=eps)
        end
    end

    @testset "2.7 Cosine function" begin
        DifferentialAlgebra.init(k, 1)

        x = DA(1, 1)
        f = cos(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.getCoefficient(f, jj)
            if iseven(i)
                a_exact = (-1)^(i / 2) / factorial(i)
            else
                a_exact = 0.0
            end
            @test isapprox(a_da, a_exact, atol=eps)
        end
    end

    @testset "2.8 Tangent function" begin
        DifferentialAlgebra.init(k, 1)

        x = DA(1, 1)
        f = tan(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.getCoefficient(f, jj)
            if iseven(i)
                a_exact = 0.0
            else
                a_exact = (-1)^((i-1)/2) * 2^(i+1) * (2^(i+1) - 1) * bernoulli(i+1) / factorial(i+1)
            end
            @test isapprox(a_da, a_exact, atol=eps)
        end
    end

    @testset "2.9 Arcsine function" begin
        DifferentialAlgebra.init(k, 1)

        x = DA(1, 1)
        f = asin(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.getCoefficient(f, jj)
            if iseven(i)
                a_exact = 0.0
            else
                a_exact = factorial(i - 1) / (4.0^((i - 1) / 2) * factorial(div(i - 1, 2))^2 * i)
            end
            @test isapprox(a_da, a_exact, atol=eps)
        end
    end

    @testset "2.10 Arccosine function" begin
        DifferentialAlgebra.init(k, 1)

        x = DA(1, 1)
        f = acos(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.getCoefficient(f, jj)
            if i == 0
                a_exact = pi / 2.0
            elseif iseven(i)
                a_exact = 0.0
            else
                a_exact = - factorial(i - 1) / (4.0^((i - 1) / 2) * factorial(div(i - 1, 2))^2 * i)
            end
            @test isapprox(a_da, a_exact, atol=eps)
        end
    end

    @testset "2.11 Arctangent function" begin
        DifferentialAlgebra.init(k, 1)

        x = DA(1, 1)
        f = atan(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.getCoefficient(f, jj)
            if iseven(i)
                a_exact = 0.0
            else
                a_exact = (-1)^((i-1)/2) / i
            end
            @test isapprox(a_da, a_exact, atol=eps)
        end
    end

    @testset "2.12 Hyperbolic sine function" begin
        DifferentialAlgebra.init(k, 1)

        x = DA(1, 1)
        f = sinh(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.getCoefficient(f, jj)
            if iseven(i)
                a_exact = 0.0
            else
                a_exact = 1.0 / factorial(i)
            end
            @test isapprox(a_da, a_exact, atol=eps)
        end
    end

    @testset "2.13 Hyperbolic cosine function" begin
        DifferentialAlgebra.init(k, 1)

        x = DA(1, 1)
        f = cosh(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.getCoefficient(f, jj)
            if iseven(i)
                a_exact = 1.0 / factorial(i)
            else
                a_exact = 0.0
            end
            @test isapprox(a_da, a_exact, atol=eps)
        end
    end

    @testset "2.14 Hyperbolic tangent function" begin
        DifferentialAlgebra.init(k, 1)

        x = DA(1, 1)
        f = tanh(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.getCoefficient(f, jj)
            if iseven(i)
                a_exact = 0.0
                @test isapprox(a_da, a_exact, atol=eps)
            else
                # tanh(x) = sum(2^(2n)*(2^(2n)-1)*B_(2n)/(2n)! * x^(2n-1)).
                a_exact = 2^(i+1) * (2^(i+1) - 1) * bernoulli(i+1) / factorial(i+1)
                @test isapprox(a_da, a_exact, atol=eps)
            end
        end
    end

    @testset "2.15 Hyperbolic arcsine function" begin
        DifferentialAlgebra.init(k, 1)

        x = DA(1, 1)
        f = asinh(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.getCoefficient(f, jj)
            if iseven(i)
                a_exact = 0.0
            else
                a_exact = (-1)^((i-1)/2) * factorial(i-1) / (4^((i-1)/2) * factorial(div(i-1, 2))^2 * i)
            end
            @test isapprox(a_da, a_exact, atol=eps)
        end
    end

    @testset "2.16 Hyperbolic arctangent function" begin
        DifferentialAlgebra.init(k, 1)

        x = DA(1, 1)
        f = atanh(x)

        jj = Vector{UInt32}(undef, 1)
        for i in 0:k
            jj[1] = i
            a_da = DifferentialAlgebra.getCoefficient(f, jj)
            if iseven(i)
                a_exact = 0.0
            else
                a_exact = 1 / i
            end
            # Preserve the reference validation's 1e-12 coefficient tolerance.
            @test isapprox(a_da, a_exact, atol=1e-12)
        end
    end
end
