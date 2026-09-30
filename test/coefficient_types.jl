using Test, LinearAlgebra, Random, GenericLinearAlgebra, SpecialFunctions

function multiplication_allocations(out, a, b)
    mul!(out, a, b)
    return @allocated mul!(out, a, b)
end
function evaluation_allocations(out, map, point, work)
    DifferentialAlgebra.evaluate!(out, map, point, work)
    return @allocated DifferentialAlgebra.evaluate!(out, map, point, work)
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
        A = [T(2) + x / T(10) one(x) / T(5); one(x) / T(5) T(4) + y / T(10)]
        eigenvalues, eigenvectors = DifferentialAlgebra.eigenpairs(A)
        @test eltype(eigenvalues) === TaylorPolynomial{T}
        @test maximum(DifferentialAlgebra.coefficient_norm, A * eigenvectors - eigenvectors * Diagonal(eigenvalues)) < 10tolerance
        @test maximum(DifferentialAlgebra.coefficient_norm, eigenvectors' * eigenvectors - I) < 10tolerance
        for f in (erf, erfc, a -> besselj(2, a), a -> bessely(2, a))
            @test f(p) isa TaylorPolynomial{T}
            @test DifferentialAlgebra.constant_term(f(p)) ≈ f(T(2)) rtol = tolerance
        end
        for f in (gamma, loggamma, digamma)
            @test f(TaylorPolynomial{T}(2)) isa TaylorPolynomial{T}
            @test DifferentialAlgebra.constant_term(f(TaylorPolynomial{T}(2))) ≈ f(T(2)) rtol = tolerance
        end
        if T !== BigFloat
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
    @test isfinite(DifferentialAlgebra.random_polynomial())
    @test_throws ArgumentError DifferentialAlgebra.besselj(typemin(Int), TaylorPolynomial(1.0))
end
