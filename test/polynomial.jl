using LinearAlgebra, Random, SpecialFunctions

# TaylorPolynomial equality intentionally compares constants. Check every coefficient instead.
polyerror(a::TaylorPolynomial, b::Real) = DifferentialAlgebra.coefficient_norm(a - b, 0)
polyerror(a::AbstractArray, b::AbstractArray) = maximum(polyerror.(a, b); init = 0.0)
@noinline temporary_buffer(p) = WeakRef(sin(p).coeffs)

# A valid AbstractVector whose elements have no owner until materialized.
struct GeneratedPolynomials <: AbstractVector{TaylorPolynomial}
    n::Int
end
Base.size(v::GeneratedPolynomials) = (v.n,)
function Base.getindex(v::GeneratedPolynomials, i::Int)
    checkbounds(v, i)
    GC.gc(false)
    return TaylorPolynomial(i, 1.0)^2
end

@testset "Polynomial storage and arithmetic" begin
    @test DifferentialAlgebra.eval(:(1 + 2)) == 3
    @testset "Constructors, coefficients and ownership" begin
        DifferentialAlgebra.initialize!(5, 3)
        DifferentialAlgebra.set_coefficient_tolerance!(0.0)
        x, y, z = variable.(1:DifferentialAlgebra.nvariables())
        @test DifferentialAlgebra.constant_term(TaylorPolynomial(2)) == 2.0
        @test DifferentialAlgebra.constant_term(TaylorPolynomial(0, 2.5)) == 2.5
        @test DifferentialAlgebra.constant_term(TaylorPolynomial(true)) == 1.0
        @test DifferentialAlgebra.constant_term(TaylorPolynomial(1 // 2)) == 0.5
        p = 2 + 3x + 4x * y + 5z^3
        @test DifferentialAlgebra.nmonomials() == binomial(8, 3)
        @test DifferentialAlgebra.coefficient(p, [1, 0, 0]) == 3
        @test_throws DimensionMismatch DifferentialAlgebra.coefficient(p, [1, 1, 0, 9])
        @test_throws ArgumentError DifferentialAlgebra.coefficient(p, [6, 0, 0])
        @test DifferentialAlgebra.linear_part(p) == [3, 0, 0]
        q = copy(p)
        DifferentialAlgebra.set_coefficient!(q, [1, 1], 7)
        @test DifferentialAlgebra.coefficient(p, [1, 1, 0]) == 4
        @test DifferentialAlgebra.coefficient(q, [1, 1, 0]) == 7
        ms = DifferentialAlgebra.monomials(p)
        @test DifferentialAlgebra.degree.(ms) == [0, 1, 2, 3]
        @test DifferentialAlgebra.coefficient.(ms) == [2, 3, 4, 5]
        @test DifferentialAlgebra.exponents(last(ms)) == [0, 0, 3]
        @test DifferentialAlgebra.coefficient(DifferentialAlgebra.monomial(p, 1)) == 2
        @test_throws BoundsError DifferentialAlgebra.monomial(p, 0)
        @test_throws ArgumentError DifferentialAlgebra.set_coefficient!(p, [-1], 2)
        @test_throws ArgumentError TaylorPolynomial(4, 1.0)
        @test_throws ArgumentError DifferentialAlgebra.variable(0)
        # The coefficient storage must stay alive across GC.
        for _ in 1:3
            GC.gc(true)
            @test polyerror(p, 2 + 3x + 4x * y + 5z^3) == 0
        end
        weakbuffer = temporary_buffer(p)
        GC.gc(true)
        @test weakbuffer.value === nothing
        @test Base.summarysize(p) >= sizeof(p.coeffs)
        grown = TaylorPolynomial()
        # Exercise coefficient storage growth, including across GC.
        indices = DifferentialAlgebra.multiindices(5, 3)
        for (i, jj) in enumerate(indices)
            DifferentialAlgebra.set_coefficient!(grown, jj, i)
            i % 7 == 0 && GC.gc(false)
        end
        @test DifferentialAlgebra.coefficient.(Ref(grown), indices) == collect(1:length(indices))
        @test polyerror(copy(grown), grown) == 0
        compiled = DifferentialAlgebra.compile(p)
        DifferentialAlgebra.initialize!(2, 1)
        @test_throws ArgumentError p + 1
        @test_throws ArgumentError DifferentialAlgebra.evaluate(compiled, [TaylorPolynomial(1, 1.0)])
        # Numeric evaluation owns its data and does not depend on algebra state.
        @test only(DifferentialAlgebra.evaluate(compiled, [0.1, 0.2, 0.3])) ≈ 2.515
        @test_throws ArgumentError DifferentialAlgebra.initialize!(0, 1)
        @test_throws ArgumentError DifferentialAlgebra.initialize!(typemax(Int), 2)
        @test DifferentialAlgebra.constant_term(TaylorPolynomial(7)) == 7 # Invalid init arguments leave the algebra intact.
    end

    @testset "Scalar API and error recovery" begin
        DifferentialAlgebra.initialize!(5, 2)
        DifferentialAlgebra.set_coefficient_tolerance!(0.0)
        x, y = variable.(1:DifferentialAlgebra.nvariables())
        p = 2 + x + 0.2y
        for a in (2, 2.0, 2 // 1, π, ℯ)
            @test polyerror(a^(1 + x), exp(log(a) * (1 + x))) < 1.0e-14
        end
        @test polyerror(p^(1 + y), exp(log(p) * (1 + y))) < 1.0e-14
        @test polyerror(p^(1 // 2), sqrt(p)) < 1.0e-14
        @test polyerror(log(3.0, p), log(p) / log(3.0)) < 1.0e-14
        @test polyerror(log(ℯ, p), log(p)) == 0
        @test polyerror(atan(x, p), atan(x / p)) < 1.0e-14
        @test polyerror(hypot(x, p), sqrt(x * x + p * p)) < 1.0e-14
        @test polyerror(mod(p + 3, 3), p) == 0
        @test DifferentialAlgebra.constant_term(round(TaylorPolynomial(2.3) + x)) == 2
        @test DifferentialAlgebra.constant_term(trunc(TaylorPolynomial(-2.3) + x)) == -2
        @test iszero(sqrt(TaylorPolynomial())) && iszero(cbrt(TaylorPolynomial()))
        @test iszero(hypot(TaylorPolynomial(), TaylorPolynomial()))
        @test isless(1.0, p) && !isless(p, 1.0)
        @test isequal(2.0, p) && isequal(p, 2.0)
        @test (p == π) == (π == p) == false
        @test !iszero(x) && x == 0 # Legacy comparisons use the constant part.
        @test !iszero(TaylorPolynomial(NaN)) && !iszero(TaylorPolynomial(Inf))
        @test iszero(p - p)
        @test_throws DifferentialAlgebra.TaylorError log(-p)
        @test_throws DifferentialAlgebra.TaylorError inv(x)
        @test_throws DifferentialAlgebra.TaylorError sqrt(x)
        @test_throws ArgumentError p^typemin(Cint)
        @test polyerror(exp(log(p)), p) < 1.0e-13
        @test polyerror(DifferentialAlgebra.differentiate(p^3, [2, 1]), 1.2) < 1.0e-13
        @test polyerror(DifferentialAlgebra.differentiate(DifferentialAlgebra.integrate(p, 1), 1), p) < 1.0e-14
        @test polyerror(DifferentialAlgebra.differentiate(DifferentialAlgebra.integrate(p, [1, 1]), [1, 1]), p) < 1.0e-14
        @test polyerror(DifferentialAlgebra.divide_variable(x^2 * y, 1, 2), y) == 0
        @test polyerror(DifferentialAlgebra.coefficient_product(2 + 3x, 4 + 5x + x * x), 8 + 15x) == 0
        @test DifferentialAlgebra.coefficient_norm(2 + 3x - 4x * y, 0) == 4
        @test DifferentialAlgebra.coefficient_norm(2 + 3x - 4x * y, 1) == 9
        @test DifferentialAlgebra.degree_norms(2 + 3x - 4x * y, 0, 1) == [2, 3, 4, 0, 0, 0]
        @test length(DifferentialAlgebra.estimate_norms(exp(x), 0, 1, 8)) == 9
        @test_logs (:warn, r"estimate") DifferentialAlgebra.estimate_norms(x)
        @test iszero(DifferentialAlgebra.coefficient_product(x, y))
        @test iszero(DifferentialAlgebra.coefficient_product(x, TaylorPolynomial()))
        # Sparse outputs reserve only as many terms as the operation can create.
        @test polyerror((2 + x) * (3 + y), 6 + 3x + 2y + x * y) == 0
        @test polyerror((2 + x) / TaylorPolynomial(2), 1 + 0.5x) == 0
        @test polyerror(TaylorPolynomial() - x, -x) == 0
        bounds = DifferentialAlgebra.bounds(2 + 3x)
        @test bounds.lower <= -1 && bounds.upper >= 5
        @test isfinite(DifferentialAlgebra.random_polynomial(; density = 0.2, scaled = false))
        @test eps(Float64) ≈ eps(Float64)
        @test DifferentialAlgebra.set_coefficient_tolerance!(1.0e-20) == 0
        @test DifferentialAlgebra.coefficient_tolerance() == 1.0e-20
        @test DifferentialAlgebra.set_truncation_order!(3) == 5
        DifferentialAlgebra.push_order!(2); DifferentialAlgebra.push_order!(1)
        @test DifferentialAlgebra.truncation_order() == 1
        DifferentialAlgebra.pop_order!(); DifferentialAlgebra.pop_order!()
        @test DifferentialAlgebra.truncation_order() == 3
        @test_throws ArgumentError DifferentialAlgebra.pop_order!()
        @test_throws ArgumentError DifferentialAlgebra.set_truncation_order!(6)
        DifferentialAlgebra.set_truncation_order!(5)
    end

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
        @test polyerror(DifferentialAlgebra.trim(DifferentialAlgebra.differentiate(loggamma(x), 1) - digamma(x), 0, 3), 0) < 1.0e-12
        @test polyerror(DifferentialAlgebra.trim(DifferentialAlgebra.differentiate(besselj(1, x), 1) - (besselj(0, x) - besselj(2, x)) / 2, 0, 3), 0) < 1.0e-12
    end

    @testset "Compiled maps and inversion" begin
        DifferentialAlgebra.initialize!(5, 3)
        DifferentialAlgebra.set_coefficient_tolerance!(0.0)
        x, y, z = variable.(1:DifferentialAlgebra.nvariables())
        f = [2 + x + 2y + x * z, z^3 + y * y]
        compiled = CompiledMap(f)
        @test (DifferentialAlgebra.noutputs(compiled), DifferentialAlgebra.nvariables(compiled), DifferentialAlgebra.degree(compiled)) == (2, 3, 3)
        @test DifferentialAlgebra.nnodes(compiled) > 1
        @test DifferentialAlgebra.evaluate(compiled, [0.1, 0.2, 0.3]) ≈ [2.53, 0.067]
        @test DifferentialAlgebra.evaluate(compiled, [0.1]) ≈ [2.1, 0.0]
        @test DifferentialAlgebra.evaluate(compiled, [0.1, 0.2, 0.3, 99]) ≈ [2.53, 0.067]
        @test DifferentialAlgebra.evaluate(f, 0.1) ≈ [2.1, 0.0]
        @test DifferentialAlgebra.evaluate(x, 0.3) == 0.3
        @test only(DifferentialAlgebra.evaluate(DifferentialAlgebra.compile(z), [1, 2, 3])) == 3
        @test DifferentialAlgebra.evaluate(DifferentialAlgebra.compile(GeneratedPolynomials(3)), [1, 2, 3]) == [1, 4, 9]
        @test DifferentialAlgebra.evaluate(TaylorPolynomial(7), Float64[]) == 7
        @test polyerror(DifferentialAlgebra.evaluate(compiled, [x, y, z]), f) == 0
        @test polyerror(DifferentialAlgebra.evaluate(compiled, Real[x, 0, 0]), [2 + x, TaylorPolynomial()]) == 0
        out, work = zeros(2), zeros(4)
        @test DifferentialAlgebra.evaluate!(out, compiled, [0.1, 0.2, 0.3], work) ≈ [2.53, 0.067]
        @test DifferentialAlgebra.evaluate!(out, compiled, [0.1, 0.2, 0.3]) === out
        @test_throws DimensionMismatch DifferentialAlgebra.evaluate!(zeros(1), compiled, [1.0, 2.0, 3.0], work)
        @test_throws ArgumentError DifferentialAlgebra.evaluate!(out, compiled, out, work)
        @test_throws ArgumentError DifferentialAlgebra.evaluate!(out, compiled, [x, y, z], work)
        @test_throws ArgumentError DifferentialAlgebra.compile(TaylorPolynomial[])
        # A coupled map with a nontrivial Jacobian and a parameter coordinate.
        f = [2x + y + 0.1x * y, x + 3y + 0.2x * x, z + 0.1x * z]
        inverse = DifferentialAlgebra.invert(f)
        @test polyerror(DifferentialAlgebra.evaluate(f, inverse), [x, y, z]) < 1.0e-12
        @test polyerror(DifferentialAlgebra.evaluate(inverse, f), [x, y, z]) < 1.0e-12
        partial = [x + 0.2x * y + z * z, y + 0.1x * z]
        inverse = DifferentialAlgebra.invert(partial)
        @test polyerror(DifferentialAlgebra.evaluate(partial, vcat(inverse, z)), [x, y]) < 1.0e-12
        @test_throws SingularException DifferentialAlgebra.invert([x, x, z])
        @test DifferentialAlgebra.truncation_order() == 5
        DifferentialAlgebra.initialize!(1, 1)
        x = TaylorPolynomial(1, 1.0)
        @test polyerror(only(DifferentialAlgebra.invert([2 + 3x])), (x - 2) / 3) < 1.0e-15
    end

    @testset "Arrays and higher-order symmetric eigenpairs" begin
        DifferentialAlgebra.initialize!(4, 2)
        DifferentialAlgebra.set_coefficient_tolerance!(1.0e-24)
        x, y = variable.(1:DifferentialAlgebra.nvariables())
        M = [2 + x y; y 4 - x]
        @test polyerror(det(M), (2 + x) * (4 - x) - y * y) < 1.0e-14
        @test polyerror(M * inv(M), TaylorPolynomial.(Matrix{Float64}(I, 2, 2))) < 1.0e-13
        @test polyerror(transpose(M), M) == 0
        @test polyerror(norm(M), sqrt(sum(abs2, M))) < 1.0e-13
        A = [2 + x 0.3 + y 0.2x; 0.3 + y 4 - y 0.1 + x * y; 0.2x 0.1 + x * y 7 + x + y]
        values, vectors = DifferentialAlgebra.eigenpairs(A)
        @test polyerror(A * vectors, vectors * Diagonal(values)) < 2.0e-12
        @test polyerror(transpose(vectors) * vectors, TaylorPolynomial.(Matrix{Float64}(I, 3, 3))) < 2.0e-12
        @test DifferentialAlgebra.truncation_order() == 4
        # Evaluate the Taylor eigenvalues against independent numeric LAPACK.
        for point in ([1.0e-3, 2.0e-3], [-2.0e-3, 1.0e-3])
            numeric = [DifferentialAlgebra.evaluate(a, point) for a in A]
            @test DifferentialAlgebra.evaluate(values, point) ≈ eigvals(Symmetric(numeric)) atol = 2.0e-13
        end
        repeated, basis = DifferentialAlgebra.eigenpairs([1 + x TaylorPolynomial(); TaylorPolynomial() 1 + y])
        @test polyerror(repeated, [1 + x, 1 + y]) == 0
        @test polyerror(basis, TaylorPolynomial.(Matrix{Float64}(I, 2, 2))) == 0
        @test_throws ArgumentError DifferentialAlgebra.eigenpairs([1 + x y; y 1 - x])
        @test_throws ArgumentError DifferentialAlgebra.eigenpairs([x y; x y])
        @test size(DifferentialAlgebra.hessian_tensor([x * x + y, x * y])) == (2, 2, 2)
        @test polyerror(DifferentialAlgebra.hessian(x * x + y)[1, 1], 2) == 0
        @test polyerror(DifferentialAlgebra.jacobian([x * x + y, x * y])[1, 2], 1) == 0
    end

    @testset "Concurrent callers" begin
        DifferentialAlgebra.initialize!(4, 2)
        x, y = variable.(1:DifferentialAlgebra.nvariables())
        expected = DifferentialAlgebra.coefficient(sin(1 + x) * exp(y), [2, 1])
        jobs = [
            Threads.@spawn begin
                for _ in 1:100
                    value = DifferentialAlgebra.coefficient(sin(1 + x) * exp(y), [2, 1])
                    value == expected || error("Cross-task algebra state corruption")
                end
                true
            end for _ in 1:8
        ]
        @test all(fetch, jobs)
    end
    # Check our methods against Julia's numeric API. Combinations of distinct AD
    # scalar types (e.g. TaylorPolynomial with ForwardDiff.Dual) need their own promotion policy.
    standard(m) = m.module in (DifferentialAlgebra, Base, Core, SpecialFunctions) || parentmodule(m.module) === Base
    @test isempty(filter(pair -> all(standard, pair), Test.detect_ambiguities(DifferentialAlgebra; recursive = true)))
end
