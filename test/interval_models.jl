using Test, DifferentialAlgebra, IntervalArithmetic, LinearAlgebra

const IA = IntervalArithmetic
const DA = DifferentialAlgebra
const IF = Interval{Float64}
const IB = Interval{BigFloat}
subset(a, b) = IA.issubset_interval(a, b)
interval_contains(a, x) = IA.in_interval(x, a)
sameinterval(a, b) = IA.isequal_interval(a, b)
guaranteed_coefficients(p) = all(IA.isguaranteed, p.coeffs[1:p.len])

@testset "Optional extension and loading" begin
    @test Base.get_extension(DA, :DifferentialAlgebraIntervalArithmeticExt) !== nothing
    for script in (
            "using DifferentialAlgebra; @assert !DifferentialAlgebra.isinitialized(); @assert !any(m -> nameof(m) == :IntervalArithmetic, values(Base.loaded_modules)); using IntervalArithmetic; @assert !DifferentialAlgebra.isinitialized()",
            "using IntervalArithmetic; using DifferentialAlgebra; @assert !DifferentialAlgebra.isinitialized()",
            "using DifferentialAlgebra; x,=variables(1;order=2); using IntervalArithmetic; @assert degree(x*x)==2; @assert truncation_order()==2",
        )
        @test success(`$(Base.julia_cmd()) --startup-file=no --project=$(dirname(Base.active_project())) -e $script`)
    end
end

@testset "Interval coefficients and stored-polynomial enclosure" begin
    for I in (IF, IB), budget in (0, 32 * 1024^2)
        T = IA.numtype(I)
        x, y = variables(I, 2; order = 4, table_bytes = budget)
        z = IA.interval(T, -1, 1)
        p = z * x
        @test degree(p) == 1
        @test DA.nterms(p) == 1
        @test sameinterval(coefficient(p, [1, 0]), z)
        @test !iszero(p)
        @test iszero(zero(x))
        @test length(monomials(p)) == 1
        @test guaranteed_coefficients(p * p)
        @test subset(IA.interval(T, -1, 1), coefficient(p * p, [2, 0]))
        @test sameinterval(coefficient(IA.interval(T, 0) * x, [1, 0]), zero(I))
        q = 2 + x / 3 + y / 5
        @test guaranteed_coefficients(q)
        @test promote_type(typeof(q), I) === typeof(q)
        @test promote_type(I, typeof(q)) === typeof(q)
        @test promote_type(typeof(q), typeof(IA.exact(1))) === typeof(q)
        @test guaranteed_coefficients(q + IA.exact(0.1))
        @test IA.isguaranteed(constant_term(convert(TaylorPolynomial{I}, IA.exact(0.1))))
        @test IA.isguaranteed(constant_term(TaylorPolynomial{I}(IA.exact(0.1))))
        @test subset(IA.interval(T, 1 // 3), coefficient(q, [1, 0]))
        @test guaranteed_coefficients(muladd(2, p, q))
        @test guaranteed_coefficients(differentiate(q, 1))
        @test guaranteed_coefficients(integrate(q, 1))
        @test guaranteed_coefficients(DA.translate_variable(q, 1, 2, 1))
        @test guaranteed_coefficients(evaluate(q, [x + y, y]))
        # Exact constants stay guaranteed in both recurrence/series division,
        # mixed compiled maps, numeric evaluation, substitution and buffer paths.
        @test guaranteed_coefficients(3 / q)
        @test guaranteed_coefficients(TaylorPolynomial{Int}(3) / q)
        @test guaranteed_coefficients(3 / TaylorPolynomial{I}(2))
        @test !guaranteed_coefficients(0.1 / q)
        @test interval_contains(constant_term(3 / q), 3 // 2)
        @test IA.isguaranteed(evaluate(q, [1, 0]))
        @test !IA.isguaranteed(evaluate(q, [1.0, 0.0]))
        integer_polynomial = TaylorPolynomial{Int}(3)
        @test IA.isguaranteed(evaluate(integer_polynomial, [IA.interval(T, 1), IA.interval(T, 0)]))
        @test all(IA.isguaranteed, compile([integer_polynomial, q]).coefficients)
        integer_polynomial += variable(1, Int) + variable(2, Int)
        @test guaranteed_coefficients(DA.translate_variable(integer_polynomial, 1, IA.interval(T, 1), IA.interval(T, 0)))
        interval_buffer = zero(q)
        mul!(interval_buffer, integer_polynomial, integer_polynomial)
        @test guaranteed_coefficients(interval_buffer)
        @test interval_contains(coefficient(interval_buffer, [1, 1]), 2)
        # No midpoint extraction or silent guarantee repair during conversion.
        ng = convert(I, 1.0)
        @test !IA.isguaranteed(constant_term(TaylorPolynomial{I}(ng)))
        @test !IA.isguaranteed(constant_term(TaylorPolynomial{I}(1.0)))
        @test !IA.isguaranteed(coefficient(x * 0.1, [1, 0]))
        @test !IA.isguaranteed(enclose(TaylorPolynomial{I}(ng), [z, z]))
        ngzero = convert(I, 0)
        nz = ngzero * x
        @test !IA.isguaranteed(coefficient(nz * y, [1, 1]))
        @test !IA.isguaranteed(enclose(nz, [z, z]))
        degraded_zero = IA.setdecoration(zero(I), IA.def)
        @test IA.decoration(enclose(degraded_zero * x, [z, z])) == IA.def
        for f in (exp, log, sin, cos, sqrt, inv)
            @test guaranteed_coefficients(f(2 + x / 4 + y / 8))
        end
        for f in (log, inv, sqrt)
            @test_throws Exception f(TaylorPolynomial{I}(z) + x)
        end
        @test_throws Exception Float64(p)
        set_coefficient_tolerance!(1)
        @test degree(IA.interval(T, 1 // 1000) * x) == 1
        set_coefficient_tolerance!(0)
        @test sameinterval(enclose(x^2, [z, z]), IA.interval(T, 0, 1))
        @test subset(IA.interval(T, -1, 1), enclose(x * y, [z, z]))
        @test interval_contains(evaluate(q, [IA.interval(T, 1 // 2), IA.interval(T, 0)]), 13 // 6)
        @test_throws DimensionMismatch enclose(q, [z])
        @test_throws DimensionMismatch enclose(q, [z, z, z])
        @test_throws ArgumentError enclose(q, [IA.emptyinterval(I), z])
        @test_throws ArgumentError enclose(q, [IA.entireinterval(I), z])
        @test_throws ArgumentError enclose(q, [IA.nai(I), z])
        @test_throws ArgumentError enclose(q, [IA.bareinterval(-1, 1), z])
        @test_throws ArgumentError enclose(TaylorPolynomial{I}(sqrt(IA.interval(T, -1, 1))), [z, z])
        invalid = TaylorPolynomial{I}(IA.nai(I))
        @test IA.isnai(constant_term(invalid + 1))
        @test IA.isnai(constant_term(invalid * zero(x)))
        @test_throws ArgumentError enclose(invalid, [z, z])

        # Independent exact rational coefficient oracle, including mixed terms.
        terms_a = [((0, 0), 2 // 1), ((1, 0), 1 // 3), ((0, 1), -2 // 5), ((1, 1), 3 // 7)]
        terms_b = [((0, 0), -1 // 1), ((1, 0), 2 // 3), ((0, 2), 1 // 5)]
        a, b = zero(x), zero(x)
        for (ex, c) in terms_a
            set_coefficient!(a, collect(ex), IA.interval(T, c))
        end
        for (ex, c) in terms_b
            set_coefficient!(b, collect(ex), IA.interval(T, c))
        end
        expected = Dict{Tuple{Int, Int}, Rational{Int}}()
        for (alpha, ac) in terms_a, (beta, bc) in terms_b
            ex = (alpha[1] + beta[1], alpha[2] + beta[2])
            expected[ex] = get(expected, ex, 0 // 1) + ac * bc
        end
        product = a * b
        for ex in DA.multiindices(4, 2)
            @test interval_contains(coefficient(product, ex), get(expected, Tuple(Int.(ex)), 0 // 1))
        end
        point = (1 // 2, -1 // 4)
        exact_value = sum(c * point[1]^ex[1] * point[2]^ex[2] for (ex, c) in expected)
        @test interval_contains(enclose(product, IA.interval.(T, point)), exact_value)
        out = copy(a)
        mul!(out, out, b)
        @test interval_contains(enclose(out, IA.interval.(T, point)), exact_value)
    end
    # Stored binary value versus the desired decimal real.
    x, = variables(1; order = 2)
    p = TaylorPolynomial(0.1)
    @test interval_contains(enclose(p, [IA.interval(0)]), Rational{BigInt}(0.1))
    @test !interval_contains(enclose(p, [IA.interval(0)]), 1 // 10)
    @test interval_contains(IA.interval(1 // 10), 1 // 10)
end

@testset "Native Taylor-model arithmetic" begin
    for T in (Float64, BigFloat), budget in (0, 32 * 1024^2)
        z = IA.interval(T, -1, 1)
        x, y = taylor_models([z, z]; order = 1, table_bytes = budget)
        xx = x * x
        @test iszero(polynomial(xx))
        @test sameinterval(remainder(xx), IA.interval(T, 0, 1))
        @test sameinterval(enclose(xx), IA.interval(T, 0, 1))
        @test sameinterval(remainder(x * y), z)
        @test sameinterval(remainder((x + y)^2), IA.interval(T, -2, 4))
        @test sameinterval(enclose(x - x), IA.interval(T, 0))
        @test sameinterval(remainder(2 + x), IA.interval(T, 0))
        @test sameinterval(enclose(TaylorModel(1 // 3, x)), IA.interval(T, 1 // 3))
        a = TaylorModel(polynomial(x), IA.interval(T, 1, 2), x)
        b = TaylorModel(polynomial(y), IA.interval(T, -3, -2), x)
        # Exact formula: [-1,1] + [-3,3] + [-2,2] + [-6,-2].
        @test sameinterval(remainder(a * b), IA.interval(T, -12, 4))
        @test interval_contains(evaluate(a * b, [0, 0]), -6)
        @test interval_contains(evaluate(a * b, [0, 0]), -2)
        @test sameinterval(enclose(+a), enclose(a))
        @test sameinterval(enclose(-a), -enclose(a))
        set_coefficient_tolerance!(100)
        @test degree(x / 1000) == 1
        @test IA.isguaranteed(enclose(x / 1000))
        set_coefficient_tolerance!(0)
        with_order(0) do # zero selects order one
            @test sameinterval(remainder(x * y), z)
        end
        @test_throws ArgumentError TaylorModel(polynomial(x), x)
        @test_throws ArgumentError TaylorModel(convert(Interval{T}, 1), x)
        @test_throws ArgumentError TaylorModel(polynomial(x), convert(Interval{T}, 0), x)
        @test_throws ArgumentError TaylorModel(polynomial(x), sqrt(IA.interval(T, -1, 1)), x)
        @test_throws ArgumentError x + polynomial(x)
        @test_throws ArgumentError TaylorModel(IA.exact(polynomial(x)), x)
        @test_throws ArgumentError x < 1
        @test_throws ArgumentError x == x
        @test_throws ArgumentError Float64(x)
        @test sameinterval(enclose(compile(x)), enclose(x))
        @test_throws ArgumentError compile([x, y])
        @test_throws ArgumentError differentiate(x, 1)
        @test_throws ArgumentError invert([x, y])

        x, y = taylor_models([z, z]; order = 4, table_bytes = budget)
        exact = (1 + x + 2y) * (3 - x + y)
        @test sameinterval(remainder(exact), IA.interval(T, 0))
        for point in ((-1, -1), (0, 0), (1, 1), (1 // 2, -1 // 4))
            u, v = point
            @test interval_contains(evaluate(exact, collect(point)), (1 + u + 2v) * (3 - u + v))
        end
        with_order(1) do
            @test_throws ArgumentError x * x
            @test_throws ArgumentError x + 1
            @test_throws ArgumentError exp(x)
            @test_throws ArgumentError +x
            @test sameinterval(enclose(exact), enclose(copy(exact)))
        end
        original = enclose(exact)
        extracted = polynomial(exact)
        set_coefficient!(extracted, [0, 0], 123)
        box = domain(exact)
        box[1] = IA.interval(T, 0)
        @test sameinterval(enclose(exact), original)
        copied = copy(exact)
        @test copied._polynomial !== exact._polynomial
        @test copied._polynomial.coeffs !== exact._polynomial.coeffs
        @test sameinterval(enclose(copied + exact), enclose(exact + exact))
        deep = deepcopy(exact)
        @test deep._polynomial !== exact._polynomial
        @test sameinterval(enclose(deep + exact), enclose(exact + exact))
        asserted_polynomial = polynomial(x)
        asserted = TaylorModel(asserted_polynomial, IA.interval(T, 0), x)
        set_coefficient!(asserted_polynomial, [0, 0], 123)
        @test sameinterval(enclose(asserted), enclose(x))
        ordinary = polynomial(x)^4
        with_order(1) do
            @test_throws ArgumentError TaylorModel(ordinary, IA.interval(T, 0), x)
        end
        variables(2; order = 4)
        @test_throws ArgumentError enclose(exact)
        @test_throws ArgumentError polynomial(exact)
        @test_throws ArgumentError copy(exact)
        @test_throws ArgumentError enclose(deep)
    end
end

@testset "Physical coordinates, fixed axes and boundaries" begin
    for T in (Float64, BigFloat)
        box = [IA.interval(T, 2, 4), IA.interval(T, 7)]
        x, fixed = taylor_models(box; order = 3, names = (:distance, :fixed))
        @test sameinterval(enclose(x), box[1])
        @test sameinterval(enclose(fixed), box[2])
        @test interval_contains(evaluate(x + fixed, [2, 7]), 9)
        @test IA.decoration(evaluate(x + fixed, [2, 7])) == IA.com
        @test interval_contains(evaluate(x + fixed, [4, 7]), 11)
        @test sameinterval(enclose(x, [IA.interval(T, 3, 4), box[2]]), IA.interval(T, 3, 4))
        r = TaylorModel(polynomial(x), IA.interval(T, -1, 1), x)
        @test subset(IA.interval(T, 2, 4), enclose(r, [IA.interval(T, 3), box[2]]))
        @test_throws DomainError evaluate(x, [1, 7])
        @test_throws DomainError evaluate(fixed, [3, 8])
        @test_throws DomainError enclose(x, [IA.interval(T, 1, 3), box[2]])
        @test_throws DimensionMismatch evaluate(x, [3])
        @test_throws DimensionMismatch evaluate(x, [3, 7, 0])
        @test_throws ArgumentError evaluate(x, [IA.interval(T, 2, 4), box[2]])
        @test_throws ArgumentError enclose(x, [IA.emptyinterval(Interval{T}), box[2]])
        @test_throws ArgumentError taylor_models([convert(Interval{T}, 1)]; order = 2)
        @test_throws ArgumentError taylor_models([IA.interval(Float32, -1, 1)]; order = 2)
        @test_throws ArgumentError taylor_models([IA.interval(T, 1, 2, IA.trv)]; order = 2)
        @test_throws ArgumentError evaluate(x, [NaN, 7])
        # Midpoint can round to an endpoint. Outward radius still covers both.
        lo = one(T)
        hi = nextfloat(lo)
        tiny, = taylor_models([IA.interval(T, lo, hi)]; order = 2)
        @test interval_contains(evaluate(tiny, [lo]), lo)
        @test interval_contains(evaluate(tiny, [hi]), hi)
        @test subset(IA.interval(T, lo, hi), enclose(tiny))
        huge, = taylor_models([IA.interval(T, 2^50, 2^50 + 1)]; order = 2)
        @test interval_contains(evaluate(huge + 1 // 1024, [2^50]), (big(2)^50 + 1 // 1024))
        # Distinct coordinate identities cannot be combined, even if domains
        # appear identical. Ordinary reinitialization also invalidates old models.
        previous = tiny
        other, = taylor_models([IA.interval(T, 0, 1)]; order = 2)
        @test_throws ArgumentError previous + other
    end
end

@testset "Taylor-theorem remainders and full-enclosure domain checks" begin
    for T in (Float64, BigFloat), budget in (0, 32 * 1024^2)
        h = IA.interval(T, -1 // 4, 1 // 4)
        x, = taylor_models([h]; order = 1, table_bytes = budget)
        for f in (exp, sin, cos)
            m = f(x)
            @test IA.isguaranteed(enclose(m))
            @test guaranteed_coefficients(polynomial(m))
            @test interval_contains(evaluate(m, [0]), f === cos || f === exp ? 1 : 0)
            @test IA.decoration(evaluate(m, [0])) == IA.com
        end
        # Independent analytical uniform remainder bounds. For |u|<=1/4:
        # |sin(u)-u|<=|u|^3/6; -u^2/2<=cos(u)-1<=0.
        @test subset(IA.interval(T, -1 // 384, 1 // 384), remainder(sin(x)))
        @test subset(IA.interval(T, -1 // 32, 0), remainder(cos(x)))
        # log(1+u)-u in [-1/24,0] by integrating -u/(1+u).
        @test subset(IA.interval(T, -1 // 24, 0), remainder(log(1 + x)))
        # sqrt(1+u)-1-u/2 = -u^2/(2*(1+sqrt(1+u))^2),
        # and sqrt(1+u)>=3/4, giving [-1/98,0].
        @test subset(IA.interval(T, -1 // 98, 0), remainder(sqrt(1 + x)))
        for f in (inv, log, sqrt)
            @test IA.isguaranteed(enclose(f(1 + x)))
        end
        x, = taylor_models([h]; order = 2, table_bytes = budget)
        # Exact geometric-series identity: inv(1-u)-(1+u+u^2)
        # = u^3/(1-u), uniformly in [-1/48,1/48].
        reciprocal = inv(1 - x)
        @test subset(IA.interval(T, -1 // 48, 1 // 48), remainder(reciprocal))
        @test interval_contains(coefficient(polynomial(reciprocal), [2]), 1 // 16) # x = ξ/4
        @test interval_contains(evaluate((1 - x)^(-2), [0]), 1)
        @test interval_contains(evaluate(1 / (1 - x), [0]), 1)
        @test interval_contains(evaluate(sqrt(1 + x), [0]), 1)
        @test sameinterval(enclose(sqrt(zero(x))), zero(Interval{T}))
        @test interval_contains(evaluate(inv(-2 + x), [0]), -1 // 2)
        incoming = TaylorModel(polynomial(x), IA.interval(T, -1 // 8, 1 // 8), x)
        @test subset(enclose(exp(x)), enclose(exp(incoming)))
        wide, = taylor_models([IA.interval(T, -2, 2)]; order = 3)
        for f in (inv, log, sqrt)
            @test_throws DomainError f(1 + wide)
        end
        @test_throws DomainError sqrt(wide * wide)
        @test_throws DomainError inv(TaylorModel(polynomial(wide) + 3, IA.interval(T, -2, 0), wide))
    end
end

@testset "Independent retained and discarded rational product oracle" begin
    for T in (Float64, BigFloat), budget in (0, 32 * 1024^2)
        z = IA.interval(T, -1, 1)
        x, y = taylor_models([z, z]; order = 2, table_bytes = budget)
        terms_a = [((0, 0), 1 // 1), ((1, 0), 1 // 3), ((0, 1), -2 // 5), ((1, 1), 3 // 7), ((2, 0), 1 // 11)]
        terms_b = [((0, 0), -1 // 1), ((1, 0), 2 // 3), ((0, 2), 1 // 5)]
        a, b = zero(polynomial(x)), zero(polynomial(x))
        for (ex, c) in terms_a
            set_coefficient!(a, collect(ex), IA.interval(T, c))
        end
        for (ex, c) in terms_b
            set_coefficient!(b, collect(ex), IA.interval(T, c))
        end
        expected = Dict{Tuple{Int, Int}, Rational{Int}}()
        for (alpha, ac) in terms_a, (beta, bc) in terms_b
            ex = (alpha[1] + beta[1], alpha[2] + beta[2])
            expected[ex] = get(expected, ex, 0 // 1) + ac * bc
        end
        product = TaylorModel(a, IA.interval(T, 0), x) * TaylorModel(b, IA.interval(T, 0), x)
        retained = polynomial(product)
        for ex in DA.multiindices(2, 2)
            @test interval_contains(coefficient(retained, ex), get(expected, Tuple(Int.(ex)), 0 // 1))
        end
        lower, upper = 0 // 1, 0 // 1
        for (ex, c) in expected
            sum(ex) > 2 || continue
            if all(iseven, ex)
                lower += min(0, c)
                upper += max(0, c)
            else
                lower -= abs(c)
                upper += abs(c)
            end
        end
        # These exact rational monomial bounds enclose the entire discarded
        # polynomial uniformly. Ordered-product bounding can be wider because
        # it does not first collect and cancel its discarded coefficients.
        @test interval_contains(remainder(product), lower)
        @test interval_contains(remainder(product), upper)
    end
end

@testset "Validated rounding policy" begin
    x, = taylor_models([IA.interval(-1, 1)]; order = 2)
    p = polynomial(x)
    for mode in (:none, :ulp)
        IA.configure(; rounding = mode)
        try
            @test_throws ArgumentError enclose(p, [IA.interval(-1, 1)])
            @test_throws ArgumentError enclose(x)
            @test_throws ArgumentError enclose(x, [IA.interval(0)])
            @test_throws ArgumentError evaluate(x, [0])
            @test_throws ArgumentError x + x
            @test_throws ArgumentError x * x
            @test_throws ArgumentError inv(2 + x)
            @test_throws ArgumentError taylor_models([IA.interval(-1, 1)]; order = 2)
        finally
            IA.configure(; rounding = :correct)
        end
    end
    @test sameinterval(enclose(x * x), IA.interval(0, 1))
end
