using Test, LinearAlgebra, SpecialFunctions, ForwardDiff, GenericLinearAlgebra, DifferentialAlgebra

@testset "Mixed precision reusable buffers" begin
    initialize!(1, 1)
    a, b = TaylorPolynomial(16777217.0), TaylorPolynomial(-16777216.0)
    out = TaylorPolynomial{Float32}(0)
    # Preserve operand precision before the final assignment to the buffer.
    # Narrowing either operand first loses this exactly representable result.
    DifferentialAlgebra.weighted_sum!(out, a, 1.0, b, 1.0)
    @test constant_term(out) == 1.0f0
    DifferentialAlgebra.scale!(out, TaylorPolynomial(1.0e40), 1.0e-40)
    @test constant_term(out) == 1.0f0
end

@testset "Working-order arithmetic" begin
    x, = variables(1; order = 3)
    p = x^3
    with_order(1) do
        @test (degree(p + 1), degree(p + TaylorPolynomial(1.0)), degree(1 * p)) == (0, 0, 0)
        for q in (1 + p, p - 1, 1 - p, +p, -p, p^1, p / 1, muladd(1, p, one(p)))
            @test degree(q) <= 1
        end
        @test degree(copy(p)) == 3
        @test degree(TaylorPolynomial{Float32}(p)) == 3
        @test degree(p) == 3
    end
    @test truncation_order() == 3
    @test degree(p + 1) == 3
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

# Every matrix access checks the public order, including accesses inside the
# eigenpair lifting loop. This deterministically detects hidden set_truncation_order! calls.
struct OrderCheckedMatrix{T} <: AbstractMatrix{T}
    data::Matrix{T}
    order::Int
end
Base.size(A::OrderCheckedMatrix) = size(A.data)
function Base.getindex(A::OrderCheckedMatrix, i::Int, j::Int)
    DifferentialAlgebra.truncation_order() == A.order || error("eigenpairs changed the global truncation order")
    return A.data[i, j]
end

@testset "Substitution, analytic domains and eigenpairs" begin
    for T in (Float32, Float64, BigFloat), budget in (0, 32 * 1024^2)
        DifferentialAlgebra.initialize!(6, 3; table_bytes = budget)
        x, y, z = variable.(1:DifferentialAlgebra.nvariables(), T)
        p = 3 + x + 2y + x^2 * y + 3x * z^2 + x^5
        tolerance = T === Float32 ? T(2.0e-5) : T === Float64 ? T(2.0e-13) : T(2)^(-200)
        @test DifferentialAlgebra.coefficient_norm(DifferentialAlgebra.substitute(p, 1, 2) - (5 + 6y + 6z^2 + 32)) == 0
        @test DifferentialAlgebra.coefficient_norm(DifferentialAlgebra.substitute(p, 1) - (3 + 2y)) == 0
        @test DifferentialAlgebra.coefficient_norm(DifferentialAlgebra.replace_variable(p, 1, 2, 2) - (3 + 4y + 4y^3 + 6y * z^2 + 32y^5)) == 0
        for v in 1:3
            args = [x, y, z]
            args[v] = 2args[v] + 1
            @test DifferentialAlgebra.coefficient_norm(DifferentialAlgebra.translate_variable(p, v, 2, 1) - DifferentialAlgebra.evaluate(p, args)) < tolerance
            args[v] = -[x, y, z][v]
            @test DifferentialAlgebra.coefficient_norm(DifferentialAlgebra.scale_variable(p, v, -1) - DifferentialAlgebra.evaluate(p, args)) < tolerance
        end
        @test_throws ArgumentError DifferentialAlgebra.substitute(p, 0)
        @test_throws ArgumentError DifferentialAlgebra.replace_variable(p, 1, 4)
        @test DifferentialAlgebra.coefficient_dot(p, 2x + 3y + x^2 * y) == T(9)
        @test DifferentialAlgebra.coefficient_norm(DifferentialAlgebra.filter_terms(p, x + y) - (x + 2y)) == 0
        @test DifferentialAlgebra.nterms(p) == 6
        @test DifferentialAlgebra.coefficient_norm(p) == 3
        @test DifferentialAlgebra.coefficient_norm(DifferentialAlgebra.monomial([2, 1, 0], T(3)) - 3x^2 * y) == 0
        @test DifferentialAlgebra.nterms(DifferentialAlgebra.filled(T(2))) == DifferentialAlgebra.nmonomials()
        @test DifferentialAlgebra.filled(T(2)) isa TaylorPolynomial{T}
        q = inv(1 - x / 4)
        estimates, errors = DifferentialAlgebra.estimate_norms(q, 0, 0, 8; errors = true)
        @test maximum(abs.(estimates - [T(1) / 4^i for i in 0:8])) < tolerance
        @test all(>=(0), errors) && length(errors) == 7
        @test DifferentialAlgebra.convergence_radius(q, T(1) / 2^14) ≈ one(T) rtol = tolerance
        @test_throws ArgumentError DifferentialAlgebra.convergence_radius(q, -1)
        @test DifferentialAlgebra.isinitialized()
        DifferentialAlgebra.set_truncation_order!(2); DifferentialAlgebra.push_order!()
        @test DifferentialAlgebra.truncation_order() == 6
        DifferentialAlgebra.pop_order!()
        @test DifferentialAlgebra.truncation_order() == 2
        DifferentialAlgebra.set_truncation_order!()
        @test DifferentialAlgebra.truncation_order() == 6
        @test DifferentialAlgebra.coefficient_norm((p / q) * q - p) < tolerance
        @test DifferentialAlgebra.coefficient_norm(muladd(T(2), p, q) - (2p + q)) < tolerance
        @test DifferentialAlgebra.coefficient_norm(muladd(p, T(2), q) - (2p + q)) < tolerance
        @test DifferentialAlgebra.coefficient_norm(muladd(p, q, p) - (p * q + p)) < tolerance
        @test DifferentialAlgebra.coefficient_norm(sqrt(2 + x + y * z)^2 - (2 + x + y * z)) < tolerance
        @test asin(TaylorPolynomial{T}(1)) isa TaylorPolynomial{T}
        @test DifferentialAlgebra.constant_term(acos(TaylorPolynomial{T}(-1))) == acos(-one(T))
        @test iszero(acosh(TaylorPolynomial{T}(1)))
        @test_throws DifferentialAlgebra.TaylorError asin(1 + x)
        A = OrderCheckedMatrix([2 + x / 10 one(x) / 5;one(x) / 5 4 + y / 10], 6)
        values, vectors = DifferentialAlgebra.eigenpairs(A)
        @test maximum(DifferentialAlgebra.coefficient_norm, A * vectors - vectors * Diagonal(values)) < 20tolerance
        @test maximum(DifferentialAlgebra.coefficient_norm, vectors' * vectors - I) < 20tolerance
        @test isempty(first(DifferentialAlgebra.eigenpairs(Matrix{TaylorPolynomial{T}}(undef, 0, 0))))
        @test size(DifferentialAlgebra.linear_part(TaylorPolynomial{T}[])) == (0, 3)
        @test size(DifferentialAlgebra.hessian_tensor(TaylorPolynomial{T}[])) == (3, 3, 0)
    end

    DifferentialAlgebra.initialize!(5, 2)
    x, y = variable.(1:DifferentialAlgebra.nvariables())
    @test DifferentialAlgebra.constant_term(hypot(TaylorPolynomial(3.0e200), TaylorPolynomial(4.0e200))) ≈ 5.0e200
    @test DifferentialAlgebra.constant_term(hypot(TaylorPolynomial(3.0e-200), TaylorPolynomial(4.0e-200))) ≈ 5.0e-200
    @test DifferentialAlgebra.constant_term(round(TaylorPolynomial(Inf))) == Inf
    @test DifferentialAlgebra.constant_term(trunc(TaylorPolynomial(-Inf))) == -Inf
    @test DifferentialAlgebra.coefficient_norm(mod(-1.5 + x, 2) - (0.5 + x)) == 0
    setprecision(256) do
        DifferentialAlgebra.set_coefficient_tolerance!(BigFloat(1) - BigFloat(2)^(-100))
        @test DifferentialAlgebra.constant_term(TaylorPolynomial(1.0)) == 1
        @test DifferentialAlgebra.constant_term(TaylorPolynomial(1.0f0)) == 1
        DifferentialAlgebra.set_coefficient_tolerance!(0)
    end
    @test DifferentialAlgebra.coefficient((1.0e308x + 1.0e-308y)^2, [1, 1]) ≈ 2
    @test_throws ArgumentError DifferentialAlgebra.multiindices(1000, 1000)
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
