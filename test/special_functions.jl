using Test, DifferentialAlgebra
using SpecialFunctions
using QuadGK: quadgk

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
            @test_throws DifferentialAlgebra.TaylorError polygamma(2, zero(x))
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
            @test_throws DifferentialAlgebra.TaylorError besselk(0, -one(x))
        end
    end
    DifferentialAlgebra.initialize!(4, 1)
    x = DifferentialAlgebra.variable(1, Float32)
    for f in (gamma, loggamma, digamma, p -> polygamma(2, p), p -> besseli(2, p), p -> besselk(2, p), p -> besselix(2, p), p -> besselkx(2, p))
        @test f(2 + x / 10) isa TaylorPolynomial{Float32}
    end
end

# Differentiate the convergent series at the origin term by term. Exact integer
# binomial/factorial factors and 1024-bit arithmetic provide an oracle independent
# of the neighboring-order recurrence used by the polynomial implementation.
function origin_bessel_coefficient(n, center, k, modified)
    value = zero(center)
    first_term = max(0, cld(k - n, 2))
    for m in first_term:(first_term + 500)
        degree = n + 2m
        numerator = binomial(big(degree), k) * center^(degree - k)
        denominator = big(2)^degree * factorial(big(m)) * factorial(big(n + m))
        term = numerator / denominator
        !modified && isodd(m) && (term = -term)
        value += term
        abs(term) <= eps(center) * abs(value) && return value
    end
    error("Independent Bessel series did not converge")
end

@testset "Regular Bessel functions near the origin" begin
    for budget in (0, 32 * 1024^2)
        x, = variables(1; order = 20, table_bytes = budget)
        for center in (1.0e-3, 0.05, 0.1), n in 0:3
            for (f, modified) in ((besselj, false), (besseli, true))
                p = f(n, center + x)
                expected = setprecision(BigFloat, 1024) do
                    [origin_bessel_coefficient(n, BigFloat(center), k, modified) for k in 0:20]
                end
                for k in 0:20
                    @test coefficient(p, [k]) ≈ Float64(expected[k + 1]) rtol = 2.0e-12
                end
                # Small steps keep truncation below the numerical regression tolerance.
                step = 0.1
                reference = setprecision(BigFloat, 1024) do
                    origin_bessel_coefficient(n, BigFloat(center) + step, 0, modified)
                end
                @test p([step]) ≈ Float64(reference) rtol = 2.0e-12
            end
        end
        tiny = besselj(2, 1.0e-160 + x)
        @test all(isfinite, tiny.coeffs[1:tiny.len])
        @test coefficient(tiny, [2]) ≈ 1 / 8
        @test besselj(3, 1.0e-3 + x)([0.1]) ≈ besselj(3, 0.101) rtol = 1.0e-12
        @test besseli(2, 0.1 + x)([1.0]) ≈ besseli(2, 1.1) rtol = 1.0e-12
    end
end

@testset "Gamma-family poles" begin
    for budget in (0, 32 * 1024^2), T in (Float32, Float64, BigFloat)
        x, = variables(T, 1; order = 3, table_bytes = budget)
        for center in (0, -1, -2), slope in (0, 1)
            p = T(center) + slope * x
            for f in (gamma, loggamma, digamma, p -> polygamma(2, p))
                @test_throws DifferentialAlgebra.TaylorError f(p)
            end
        end
        @test_throws ArgumentError polygamma(-1, x)
    end
end

@testset "Real Bessel singular domains" begin
    for budget in (0, 32 * 1024^2), T in (Float32, Float64, BigFloat)
        x, = variables(T, 1; order = 3, table_bytes = budget)
        for center in (-1, -2), slope in (0, 1), n in (0, 2)
            p = T(center) + slope * x
            for f in (bessely, besselk, besselkx)
                @test_throws DifferentialAlgebra.TaylorError f(n, p)
            end
        end
    end
end

# Independent scalar K0/K1 integrals, positive order recurrence, and the exact
# neighboring-order derivative identity. Neither the polynomial ODE nor its
# vector derivative quadrature enters this 1024-bit coefficient oracle.
# https://dlmf.nist.gov/10.32.E9 and https://dlmf.nist.gov/10.29
function scaled_k_reference(center, max_order)
    values = zeros(BigFloat, max_order + 1)
    for n in 0:1
        integrand(t) = begin
            d = center * 2sinh(t / 2)^2
            isfinite(d) ? exp(-d) * cosh(n * t) : zero(center)
        end
        value, error = quadgk(
            integrand, zero(center), BigFloat(Inf);
            rtol = BigFloat(2)^(-750), order = 64
        )
        @test error <= BigFloat(2)^(-700) * abs(value)
        values[n + 1] = value
    end
    for n in 1:(max_order - 1)
        values[n + 2] = values[n] + 2n * values[n + 1] / center
    end
    return values
end
function scaled_k_coefficient(values, n, k)
    value = zero(first(values))
    for j in 0:k
        derivative = sum(
            binomial(big(j), m) * values[abs(n - j + 2m) + 1]
                for m in 0:j
        )
        value += binomial(big(k), j) * (-1 // big(2))^j * derivative
    end
    return value / factorial(big(k))
end

@testset "Scaled Bessel K across moderate centers and retained orders" begin
    centers = (2, 5, 15, 30, 50, 64, 65)
    references, endpoints = setprecision(1024) do
        coefficients = Dict{Tuple{Int, Int}, Vector{BigFloat}}()
        endpoints = Dict{Tuple{Int, Int}, BigFloat}()
        for center in centers
            values = scaled_k_reference(BigFloat(center), 33)
            endpoint = scaled_k_reference(BigFloat(3center) / 2, 2)
            for n in (0, 2)
                coefficients[(center, n)] = [
                    scaled_k_coefficient(values, n, k) for k in 0:31
                ]
                endpoints[(center, n)] = endpoint[n + 1]
            end
        end
        coefficients, endpoints
    end
    setprecision(128) do
        for T in (Float32, Float64, BigFloat), budget in (0, 32 * 1024^2)
            for order in (12, 30)
                x, = variables(T, 1; order, table_bytes = budget)
                for center in centers, n in (0, 2)
                    result = besselkx(n, T(center) + x)
                    expected = references[(center, n)]
                    tolerance = 500eps(T)
                    # Float32 storage may round the smallest subnormal to zero.
                    floor = T === Float32 ? nextfloat(zero(T)) : zero(T)
                    @test all(0:order) do k
                        isapprox(
                            coefficient(result, [k]), T(expected[k + 1]);
                            rtol = tolerance, atol = floor
                        )
                    end
                    step = BigFloat(center) / 2
                    # Complete monotonicity bounds the remainder at positive
                    # steps by the first omitted term, independently of the kernel.
                    tail = abs(expected[order + 2]) * step^(order + 1)
                    endpoint = endpoints[(center, n)]
                    @test abs(result(T(step)) - endpoint) <=
                        tail + tolerance * abs(endpoint)
                end
            end
        end
    end
end

@testset "Scaled Bessel K at large centers" begin
    references = setprecision(1024) do
        references = Dict{Tuple{Int, Int}, Vector{BigFloat}}()
        for center in (100, 600)
            values = scaled_k_reference(BigFloat(center), 113)
            for n in (0, 2, 100)
                references[(center, n)] = [
                    scaled_k_coefficient(values, n, k) for k in 0:13
                ]
            end
        end
        references
    end
    endpoint = setprecision(1024) do
        values = scaled_k_reference(BigFloat(630), 2)
        values[3]
    end
    for bits in (128, 256)
        setprecision(bits) do
            types = bits == 128 ? (Float32, Float64, BigFloat) : (BigFloat,)
            for T in types, budget in (0, 32 * 1024^2)
                x, = variables(T, 1; order = 12, table_bytes = budget)
                tolerance = 500eps(T)
                for center in (100, 600), n in (0, 2, 100)
                    result = besselkx(n, T(center) + x)
                    @test coefficient_type(result) === T
                    for k in 0:12
                        @test isapprox(
                            coefficient(result, [k]), references[(center, n)][k + 1];
                            rtol = tolerance
                        )
                    end
                end
                result = besselkx(2, T(600) + x)
                tail = abs(references[(600, 2)][14]) * big(30)^13
                @test abs(result(T(30)) - endpoint) <= tail + tolerance * abs(endpoint)
                @test truncation_order() == 12 && iszero(coefficient_tolerance())
                argument = T(600) + T(1 // 100) * x
                set_coefficient_tolerance!(0.01)
                result = besselkx(100, argument)
                expected = references[(600, 100)][2] / 100
                @test coefficient(result, [1]) ≈ expected rtol = tolerance
                @test coefficient_tolerance() == 0.01
            end
        end
    end
end

@testset "Bessel poles reject constant and varying inputs" begin
    for T in (Float32, Float64, BigFloat), budget in (0, 32 * 1024^2)
        x, = variables(T, 1; order = 3, table_bytes = budget)
        for f in (bessely, besselk, besselkx), p in (zero(x), x)
            @test_throws DifferentialAlgebra.TaylorError f(0, p)
        end
    end
end
