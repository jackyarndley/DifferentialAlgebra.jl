using Test, DifferentialAlgebra
using LinearAlgebra, Random, SpecialFunctions

include("support/allocations.jl")

function multiplication_allocations(out, a, b)
    mul!(out, a, b)
    return @allocated mul!(out, a, b)
end
function evaluation_allocations(out, map, point, work)
    DifferentialAlgebra.evaluate!(out, map, point, work)
    return @allocated DifferentialAlgebra.evaluate!(out, map, point, work)
end

@testset "Stable coefficient norms" begin
    for T in (Int, Float32, Float64, BigFloat)
        x, y = variables(T, 2; order = 2)
        p = 3x + 4y
        R = typeof(float(one(T)))
        for (power, expected) in ((0, 4), (1, 7), (2, 5), (3, cbrt(R(91))))
            value = coefficient_norm(p, power)
            grouped = DifferentialAlgebra.degree_norms(p, 0, power)
            @test value isa R
            @test grouped isa Vector{R}
            @test value ≈ R(expected)
            @test grouped ≈ R[0, expected, 0]
        end
        @test_throws ArgumentError coefficient_norm(p, -1)
        @test_throws ArgumentError coefficient_norm(p, 1.5)
        @test_throws ArgumentError DifferentialAlgebra.degree_norms(p, 0, -1)
        scales = T === Int ? (typemax(Int) ÷ 8,) : (sqrt(floatmax(T)), sqrt(floatmin(T)))
        for scale in scales
            p = scale * x + scale * y
            for power in (1, 2, 3)
                # Compare scaled answers: neither the oracle nor the test forms
                # potentially overflowing or underflowing coefficient powers.
                expected = R(2)^(one(R) / power)
                @test coefficient_norm(p, power) / scale ≈ expected
                @test DifferentialAlgebra.degree_norms(p, 0, power)[2] / scale ≈ expected
                @test DifferentialAlgebra.degree_norms(p, 1, power)[1:2] ≈ R[scale, scale]
            end
        end
    end
end

@testset "Typed Julia coefficients" begin
    for T in (Float32, Float64, BigFloat)
        DifferentialAlgebra.initialize!(5, 3)
        x, y, z = variable.(1:DifferentialAlgebra.nvariables(), T)
        tolerance = T === Float32 ? T(2.0e-5) : T === Float64 ? T(2.0e-13) : T(2)^(-200)
        p = T(2) + x / T(3) + y / T(5) + z * z / T(7)
        @test x isa TaylorPolynomial{T}
        @test (x + 1) * (2 - y) isa TaylorPolynomial{T}
        @test x / 3 isa TaylorPolynomial{T}
        @test DifferentialAlgebra.coefficient(x / 3, [1, 0, 0]) == one(T) / 3
        for S in (Bool, Int, BigInt, Float32, Float64, BigFloat, Rational{Int}, Rational{BigInt})
            @test promote_type(TaylorPolynomial{T}, S) === TaylorPolynomial{promote_type(T, S)}
            @test promote_type(S, TaylorPolynomial{T}) === TaylorPolynomial{promote_type(T, S)}
        end
        @test promote_type(TaylorPolynomial, TaylorPolynomial{T}) === TaylorPolynomial
        @test promote_type(TaylorPolynomial{T}, TaylorPolynomial) === TaylorPolynomial
        @test @inferred(sin(p)) isa TaylorPolynomial{T}
        @test @inferred(exp(p)) isa TaylorPolynomial{T}
        for f in (sin, cos, tan, exp, log, log2, log10, sqrt, cbrt, sinh, cosh, tanh, asinh, acosh, inv)
            value = f(p)
            @test value isa TaylorPolynomial{T}
            @test DifferentialAlgebra.constant_term(value) ≈ f(T(2)) rtol = tolerance
        end
        for f in (asin, acos, atan, atanh)
            @test f(x / T(3)) isa TaylorPolynomial{T}
        end
        @test DifferentialAlgebra.coefficient_norm(exp(log(p)) - p) < tolerance
        @test DifferentialAlgebra.coefficient_norm(sin(p)^2 + cos(p)^2 - 1) < tolerance
        @test DifferentialAlgebra.coefficient_norm(sqrt(p)^2 - p) < tolerance
        @test DifferentialAlgebra.coefficient_norm(p * inv(p) - 1) < tolerance
        @test DifferentialAlgebra.differentiate(p, 1) isa TaylorPolynomial{T}
        @test DifferentialAlgebra.integrate(p, 1) isa TaylorPolynomial{T}
        @test DifferentialAlgebra.trim(p, 1) isa TaylorPolynomial{T}
        @test DifferentialAlgebra.degree_norms(p) isa Vector{T}
        @test DifferentialAlgebra.raw_moments(exp(x), 3)[2] isa Vector{T}
        @test DifferentialAlgebra.evaluate(p, T[0.1, 0.2, 0.3]) isa T
        map = DifferentialAlgebra.compile([p, sin(p)])
        @test DifferentialAlgebra.evaluate(map, T[0.1, 0.2, 0.3]) isa Vector{T}
        @test DifferentialAlgebra.coefficient_norm(DifferentialAlgebra.evaluate(p, [x, y, z]) - p) < tolerance
        f = [x + T(0.02) * y * y, y + T(0.03) * z * z, z + T(0.01) * x * y]
        inverse = DifferentialAlgebra.invert(f)
        @test eltype(inverse) === TaylorPolynomial{T}
        @test maximum(DifferentialAlgebra.coefficient_norm, DifferentialAlgebra.evaluate(f, inverse) - [x, y, z]) < tolerance
        for f in (erf, erfc, a -> besselj(2, a), a -> bessely(2, a))
            @test f(p) isa TaylorPolynomial{T}
            @test DifferentialAlgebra.constant_term(f(p)) ≈ f(T(2)) rtol = tolerance
        end
        for f in (gamma, loggamma, digamma)
            @test f(TaylorPolynomial{T}(2)) isa TaylorPolynomial{T}
            @test DifferentialAlgebra.constant_term(f(TaylorPolynomial{T}(2))) ≈ f(T(2)) rtol = tolerance
        end
        if CHECK_ZERO_ALLOCATIONS && T !== BigFloat
            out = zero(p)
            @test multiplication_allocations(out, p, p) == 0
            result, work = zeros(T, 2), zeros(T, DifferentialAlgebra.degree(map) + 1)
            @test evaluation_allocations(result, map, T[0.1, 0.2, 0.3], work) == 0
        end
    end

    setprecision(256) do
        DifferentialAlgebra.initialize!(6, 2)
        x, y = variable.(1:DifferentialAlgebra.nvariables(), BigFloat)
        c = one(BigFloat) + BigFloat(2)^(-150)
        @test DifferentialAlgebra.constant_term(TaylorPolynomial{BigFloat}(c)) == c != BigFloat(Float64(c))
        @test DifferentialAlgebra.coefficient(exp(x), [6, 0]) ≈ inv(BigFloat(factorial(6))) rtol = big"1e-70"
        tiny = big"1e-500"
        @test DifferentialAlgebra.coefficient(tiny * x, [1, 0]) == tiny
        DifferentialAlgebra.set_coefficient_tolerance!(big"1e-501")
        @test DifferentialAlgebra.coefficient_tolerance() == big"1e-501"
        @test DifferentialAlgebra.coefficient(tiny * x, [1, 0]) == tiny
        @test iszero(big"1e-502" * x)
        DifferentialAlgebra.set_coefficient_tolerance!(0.0)
        @test DifferentialAlgebra.evaluate(x, [big"0.123456789012345678901234567890"]) == big"0.123456789012345678901234567890"
        # Abstract argument vectors must preserve the precision of TaylorPolynomial entries.
        map = DifferentialAlgebra.compile([DifferentialAlgebra.variable(1, Float64)])
        composed = only(DifferentialAlgebra.evaluate(map, Real[c * x, 0]))
        @test composed isa TaylorPolynomial{BigFloat}
        @test DifferentialAlgebra.coefficient(composed, [1, 0]) == c
        small = TaylorPolynomial{Float32}(1)
        @test small^2.0 isa TaylorPolynomial{Float64}
        @test small^big"2.0" isa TaylorPolynomial{BigFloat}
        @test small^(1 // 2) isa TaylorPolynomial{Float32}
        @test atan(small, c) isa TaylorPolynomial{BigFloat}
        @test hypot(c, small) isa TaylorPolynomial{BigFloat}
        @test DifferentialAlgebra.constant_term(atan(small, c)) == atan(one(BigFloat), c)
        @test DifferentialAlgebra.constant_term(hypot(c, small)) ≈ hypot(c, one(BigFloat)) rtol = 4eps(BigFloat)
    end

    # Exhaustive coefficient oracle: multiply sparse integer coefficient maps by
    # adding exponent tuples, independent of either lookup/ranking implementation.
    for budget in (0, 32 * 1024^2)
        DifferentialAlgebra.initialize!(4, 3; table_bytes = budget)
        rng = Xoshiro(901)
        indices = DifferentialAlgebra.multiindices(4, 3)
        for trial in 1:12
            ac, bc = rand(rng, -2:2, length(indices)), rand(rng, -2:2, length(indices))
            a, b = TaylorPolynomial(), TaylorPolynomial()
            for (i, jj) in enumerate(indices)
                DifferentialAlgebra.set_coefficient!(a, jj, ac[i]); DifferentialAlgebra.set_coefficient!(b, jj, bc[i])
            end
            expected = Dict{NTuple{3, Int}, Int}()
            for (i, alpha) in enumerate(indices), (j, beta) in enumerate(indices)
                jj = Tuple(Int.(alpha .+ beta))
                sum(jj) <= 4 || continue
                expected[jj] = get(expected, jj, 0) + ac[i] * bc[j]
            end
            product = a * b
            @test all(DifferentialAlgebra.coefficient(product, jj) == get(expected, Tuple(Int.(jj)), 0) for jj in indices)
            inplace = copy(a); mul!(inplace, inplace, b)
            @test DifferentialAlgebra.coefficient_norm(inplace - product) == 0
            DifferentialAlgebra.add!(inplace, inplace, b)
            @test DifferentialAlgebra.coefficient_norm(inplace - (product + b)) == 0
        end
        x, y, z = variable.(1:DifferentialAlgebra.nvariables())
        @test DifferentialAlgebra.coefficient_norm(exp(log(2 + x + y)) - (2 + x + y)) < 1.0e-13
        @test DifferentialAlgebra.coefficient_norm(sqrt(2 + x)^2 - (2 + x)) < 1.0e-13
        @test maximum(DifferentialAlgebra.coefficient_norm, DifferentialAlgebra.evaluate([x * y, z * z], [x + y, y + z, z - x]) - [(x + y) * (y + z), (z - x)^2]) < 1.0e-13
        # Small Taylor factors can produce significant coefficients after scaling.
        # Applying epsilon to those factors prematurely broke the DAIOD example.
        DifferentialAlgebra.set_coefficient_tolerance!(1.0e-24)
        @test DifferentialAlgebra.coefficient(inv(1.0e8 + 1.0e6 * x), [3, 0, 0]) ≈ -1.0e-14 rtol = 1.0e-13
        DifferentialAlgebra.set_coefficient_tolerance!(0.0)
    end
    DifferentialAlgebra.initialize!(70, 1)
    @test_throws ArgumentError besselj(typemin(Int), TaylorPolynomial(1.0))
end

@testset "Coefficient conversion trims exact underflow only" begin
    x, = variables(1; order = 3)
    p = 1.0e-50x
    q = TaylorPolynomial{Float32}(p)
    @test iszero(q) && degree(q) == 0
    @test q.len == 1
    @test coefficient(p, [1]) == 1.0e-50
    p = 1 + 1.0e-8x
    set_coefficient_tolerance!(1.0e-6)
    q = TaylorPolynomial{Float32}(p)
    @test coefficient(q, [1]) == Float32(1.0e-8)
    @test degree(q) == 1
end

@testset "Scalar conversion cannot discard polynomial terms" begin
    x, y = variables(2; order = 3)
    for T in (Float32, Float64, BigFloat)
        @test T(TaylorPolynomial(2)) == T(2)
        @test convert(T, TaylorPolynomial(2)) == T(2)
        @test_throws InexactError T(1 + x)
        @test_throws InexactError convert(T, 1 + x)
        @test_throws InexactError T[1 + x, y]
        values = zeros(T, 2)
        @test_throws InexactError values[1] = 1 + x
        @test values == zeros(T, 2)
    end
end

@testset "Inferred coefficient operations" begin
    for T in (
                Int, BigInt, UInt8, UInt64, UInt128, Rational{Int}, Rational{BigInt},
                Float32, Float64, BigFloat,
            ), budget in (0, 32 * 1024^2)
        x, y = variables(T, 2; order = 4, table_bytes = budget)
        p, q = one(x) + x + T(2) * y, T(3) * one(x) + y
        @test @inferred(p + q) isa TaylorPolynomial{T}
        @test @inferred(p * q) isa TaylorPolynomial{T}
        R = typeof(inv(one(T)))
        @test @inferred(inv(q)) isa TaylorPolynomial{R}
        @test @inferred(evaluate(p, T[1, 2])) == T(6)
        @test @inferred(coefficient(p, [1, 0])) == one(T)
    end
end

@testset "Scalar working types for roots and hypotenuses" begin
    for T in (
                Int, BigInt, Rational{Int}, Rational{BigInt}, Float32, Float64, BigFloat,
            ), budget in (0, 32 * 1024^2)
        x, = variables(T, 1; order = 4, table_bytes = budget)
        p, q = T(3) + x + x^4, T(4) + 2x
        z = zero(x)
        root_cases = (
            (sqrt, T(4) + x + x^4, 4, 1 // 4, -1 // 64),
            (cbrt, T(8) + x + x^4, 8, 1 // 12, -1 // 288),
        )
        with_order(2) do
            for (f, argument, center, linear, quadratic) in root_cases
                R = typeof(f(T(center)))
                value = @inferred f(argument)
                @test value isa TaylorPolynomial{R}
                @test constant_term(value) == R(2)
                @test coefficient(value, [1]) ≈ R(linear) rtol = 8eps(R)
                @test coefficient(value, [2]) ≈ R(quadratic) rtol = 8eps(R)
                @test degree(value) == 2
                zero_root = @inferred f(z)
                @test zero_root isa TaylorPolynomial{R}
                @test iszero(zero_root)
                @test zero_root.coeffs !== z.coeffs
                set_coefficient!(value, [1], R(3))
                @test coefficient(argument, [1]) == one(T)
                @test coefficient(argument, [4]) == one(T)

                identity_root = DifferentialAlgebra.nthroot(argument, 1)
                @test identity_root isa TaylorPolynomial{T}
                @test identity_root.coeffs !== argument.coeffs
                @test coefficient(identity_root, [4]) == one(T)
            end
            R = typeof(hypot(T(3), T(4)))
            value = @inferred hypot(p, q)
            @test value isa TaylorPolynomial{R}
            # sqrt(25 + 22x + 5x^2) has these exact first coefficients.
            @test constant_term(value) == R(5)
            @test coefficient(value, [1]) ≈ R(11 // 5) rtol = 8eps(R)
            @test coefficient(value, [2]) ≈ R(2 // 125) rtol = 64eps(R)
            @test degree(value) == 2
            @test value.coeffs !== p.coeffs && value.coeffs !== q.coeffs
            @test @inferred(hypot(z, z)) isa TaylorPolynomial{R}
            @test iszero(hypot(z, z))
            for argument in (p, -p)
                for axis_value in (@inferred(hypot(argument, z)), @inferred(hypot(z, argument)))
                    @test axis_value isa TaylorPolynomial{R}
                    @test constant_term(axis_value) == R(3)
                    @test coefficient(axis_value, [1]) == one(R)
                    @test degree(axis_value) == 1
                    @test axis_value.coeffs !== argument.coeffs
                end
            end
        end
        @test coefficient(p, [4]) == one(T)
    end
    x, = variables(Float32, 1; order = 2)
    big_zero = TaylorPolynomial{BigFloat}(zero(x))
    @test @inferred(hypot(3 + x, big_zero)) isa TaylorPolynomial{BigFloat}
    @test @inferred(hypot(big_zero, 3 + x)) isa TaylorPolynomial{BigFloat}
end
