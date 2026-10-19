using Test, DifferentialAlgebra
using SpecialFunctions

special_error(p, c) = coefficient_norm(p - c, 0)

@testset "Special function derivatives" begin
    DifferentialAlgebra.initialize!(4, 1)
    DifferentialAlgebra.set_coefficient_tolerance!(0.0)
    x = 2.3 + TaylorPolynomial(1, 1.0)
    for f in (erf, erfc, gamma, loggamma, digamma)
        @test DifferentialAlgebra.constant_term(f(x)) ≈ f(2.3) rtol = 2.0e-13
    end
    for f in (besselj, bessely, besseli, besselk, besselix, besselkx), n in (0, 1, 2)
        @test DifferentialAlgebra.constant_term(f(n, x)) ≈ f(n, 2.3) rtol = 2.0e-13
    end
    @test special_error(DifferentialAlgebra.trim(DifferentialAlgebra.differentiate(loggamma(x), 1) - digamma(x), 0, 3), 0) < 1.0e-12
    @test special_error(DifferentialAlgebra.trim(DifferentialAlgebra.differentiate(besselj(1, x), 1) - (besselj(0, x) - besselj(2, x)) / 2, 0, 3), 0) < 1.0e-12
end

@testset "Complete special functions at coefficient precision" begin
    for bits in (128, 256)
        setprecision(bits) do
            DifferentialAlgebra.initialize!(4, 1)
            x = DifferentialAlgebra.variable(1, BigFloat)
            tol = BigFloat(2)^(-bits + 30)
            for n in 1:5
                expected = (-1)^(n + 1) * factorial(big(n)) * (zeta(BigFloat(n + 1)) - 1)
                @test DifferentialAlgebra.constant_term(polygamma(n, 2 + zero(x))) ≈ expected rtol = tol
                a = -BigFloat(1) / 2 + x / 10
                @test DifferentialAlgebra.coefficient_norm(polygamma(n, a + 1) - polygamma(n, a) - (-1)^n * factorial(big(n)) * a^(-n - 1)) < 100tol
            end
            # Independent 100-decimal-digit mpmath 1.4.1 reference values.
            for (f, expected) in (
                    (p -> besseli(2, p), "0.809196566007602118254159830538643220568413967382672328252350891607223341316661635156423216783493928"),
                    (p -> besselk(2, p), "0.2096093039998351971497158542154000893596098104942351105661592886876975266704431100028063620394724142"),
                )
                @test DifferentialAlgebra.constant_term(f(BigFloat(17) / 8 + zero(x))) ≈ parse(BigFloat, expected) rtol = tol
            end
            expected = "344.6756950613542964686513088383632754343010336756080598637557812072204830205585200518156277649550717"
            @test DifferentialAlgebra.constant_term(polygamma(3, -BigFloat(5) / 8 + zero(x))) ≈ parse(BigFloat, expected) rtol = tol
            for center in (BigFloat(2), -BigFloat(1) / 2)
                p = center + x / 10
                @test gamma(p) isa TaylorPolynomial{BigFloat}
                @test DifferentialAlgebra.coefficient_norm(gamma(p + 1) - p * gamma(p)) < 10tol
                @test abs(DifferentialAlgebra.evaluate(gamma(p), [BigFloat(1) / 1000]) - gamma(center + BigFloat(1) / 10000)) < BigFloat(1.0e-16)
            end
            @test_throws DomainError polygamma(2, zero(x))
            for center in (BigFloat(1) / 8, BigFloat(2), BigFloat(20), BigFloat(300)), n in (0, 2)
                p = center + x / 100
                i, k = besselix(n, p), besselkx(n, p)
                @test i isa TaylorPolynomial{BigFloat} && k isa TaylorPolynomial{BigFloat}
                @test DifferentialAlgebra.coefficient_norm(i * besselkx(n + 1, p) + besselix(n + 1, p) * k - inv(p)) < 100tol
                @test DifferentialAlgebra.coefficient_norm(besseli(n, p) * exp(-p) - i) < 100tol
                @test DifferentialAlgebra.coefficient_norm(besselk(n, p) * exp(p) - k) < 100tol
            end
            @test DifferentialAlgebra.coefficient_norm(besseli(1, -(2 + x)) + besseli(1, 2 + x)) < tol
            @test DifferentialAlgebra.coefficient_norm(besseli(0, x) - (1 + x^2 / 4 + x^4 / 64)) < tol
            @test_throws DifferentialAlgebra.TaylorError besselix(0, x)
            @test_throws DomainError besselk(0, -one(x))
        end
    end
    DifferentialAlgebra.initialize!(4, 1)
    x = DifferentialAlgebra.variable(1, Float32)
    for f in (gamma, loggamma, digamma, p -> polygamma(2, p), p -> besseli(2, p), p -> besselk(2, p), p -> besselix(2, p), p -> besselkx(2, p))
        @test f(2 + x / 10) isa TaylorPolynomial{Float32}
    end
end
