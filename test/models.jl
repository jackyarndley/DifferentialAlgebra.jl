using Test, DifferentialAlgebra, IntervalArithmetic
using ..Fixtures
include("support/interval_products.jl")
const IA = IntervalArithmetic
const DA = DifferentialAlgebra
contains(a, x) = IA.in_interval(x, a)
sameinterval(a, b) = IA.isequal_interval(a, b)

@testset "Interval norms, outward rounding and metadata" begin
    for T in (Float64, BigFloat)
        x, y = variables(Interval{T}, 2; order = 2)
        for (c, lower, upper) in (
                (interval(T, -3, 3), 4, 5),
                (interval(T, 3, 5), 5, sqrt(big(41))),
                (interval(T, 3), 5, 5),
            )
            p = c * x + 4y
            norm = coefficient_norm(p, 2)
            @test contains(norm, lower) && contains(norm, upper)
            @test IA.inf(norm) >= lower * (1 - 10eps(T))
            @test IA.sup(norm) <= upper * (1 + 10eps(T))
            @test sameinterval(norm, DA.degree_norms(p, 0, 2)[2])
        end
        for c in (IA.setdecoration(interval(T, 3), IA.def), convert(Interval{T}, T(3)))
            result = coefficient_norm(c * x + 4y, 2)
            @test IA.decoration(result) == IA.decoration(c)
            @test IA.isguaranteed(result) == IA.isguaranteed(c)
        end
        for scale in (ldexp(one(T), -1000), ldexp(one(T), 1000))
            @test sameinterval(
                coefficient_norm(3scale * x + 4scale * y, 2), interval(T, 5scale)
            )
        end
        @test !isinteger(TaylorPolynomial{Interval{T}}(convert(Interval{T}, T(2))))
    end
end

@testset "Absolute remainders independently of certificates" begin
    for budget in (0, 32 * 1024^2)
        x, y = taylor_models(
            [interval(1, 2), interval(0, 1)]; order = 2, table_bytes = budget
        )
        results = ((x + y)^3, inv(x + y + 1), log(x + y + 1), sqrt(x + y + 1))
        functions = (t -> t^3, t -> inv(t + 1), t -> log(t + 1), t -> sqrt(t + 1))
        # Exact dyadic coordinates; high-precision scalar values are independent
        # of both the model polynomial and its separate range certificate.
        for (m, f) in zip(results, functions)
            P, R = polynomial(m), remainder(m)
            for (a, b) in ((1 // 1, 0 // 1), (3 // 2, 1 // 2), (2 // 1, 1 // 1))
                value = setprecision(512) do
                    f(BigFloat(a) + BigFloat(b))
                end
                raw = P([interval(2a - 3), interval(2b - 1)]) + R
                @test contains(raw, value)
                @test contains(evaluate(m, [a, b]), value)
            end
            @test IA.isguaranteed(R) && IA.isguaranteed(enclose(m))
        end
        # Accessors, scalar shifts and coordinate domains cannot alias the model.
        before = enclose(x)
        p, box = polynomial(x), domain(x)
        set_coefficient!(p, [0, 0], 99)
        box[1] = interval(99)
        @test sameinterval(enclose(x), before)
        @test contains(enclose(x + 2), 3) && contains(enclose(x + 2), 4)
        @test_throws ArgumentError x + 0.1
        @test_throws DomainError evaluate(x, [3, 0])
        with_order(1) do
            @test_throws ArgumentError x * x
            @test sameinterval(enclose(x), before)
        end
    end
end

@testset "Wide compositions, fixed axes and model snapshots" begin
    x, = taylor_models([interval(-10, 10)]; order = 12)
    wide = cos(cos(cos(x)))
    @test all(isfinite, (IA.inf(remainder(wide)), IA.sup(remainder(wide))))
    for point in (-10, 0, 10)
        reference = cos(cos(cos(BigFloat(point))))
        @test contains(
            polynomial(wide)([interval(point // 10)]) + remainder(wide), reference
        )
    end
    x, fixed = taylor_models([interval(2, 4), interval(7)]; order = 2)
    snapshot = compile(x + fixed)
    @test contains(snapshot([2, 7]), 9) && contains(snapshot([4, 7]), 11)
    initialize!(2, 1)
    @test contains(snapshot([3, 7]), 10)
    @test_throws ArgumentError enclose(x)
    x, = taylor_models([interval(16, 1.0e10)]; order = 7)
    wide = log(log(log(x)))
    @test isfinite(IA.sup(abs(remainder(wide))))
    @test contains(evaluate(wide, [16]), log(log(log(BigFloat(16)))))
end

@testset "Coordinate certificates include the rounded normalized cover" begin
    x, = taylor_models([interval(0.1, 0.3)]; order = 2)
    # This labelled invariant test bypasses physical-domain queries and checks
    # the full declared normalized box, which may round beyond physical endpoints.
    coordinates = x._coordinates
    cover = coordinates.center[1] + coordinates.radius[1] * coordinates.normalized[1]
    @test IA.issubset_interval(cover, enclose(x))
end

@testset "Normalized powers and discarded product oracle" begin
    ext = Base.get_extension(DA, :DifferentialAlgebraIntervalArithmeticExt)
    for (a, power, endpoints) in (
            (interval(-1, 1), 4, (0, 1)), (interval(-1, 1), 3, (-1, 1)),
            (interval(0), 2, (0, 0)), (interval(1, 2), -1, (1 // 2, 1)),
        )
        @test sameinterval(ext.coordinate_power(a, power), interval(endpoints...))
    end
    # x²*y² is discarded at order two; its unit-box range is exactly [0,1].
    for budget in (0, 32 * 1024^2)
        x, y = taylor_models(
            [interval(-1, 1), interval(-1, 1)]; order = 2, table_bytes = budget
        )
        m = (x * x) * (y * y)
        @test sameinterval(remainder(m), interval(0, 1))
        @test all(iszero, coefficients(polynomial(m)))
        @test contains(enclose(m), 1)
    end
    x, y = variables(Interval{Float64}, 2; order = 2)
    aterms = [((2, 0), (-1 // 2, 1 // 4)), ((0, 1), (1 // 2, 1 // 2))]
    bterms = [((1, 0), (1 // 4, 3 // 4)), ((0, 2), (-1 // 4, -1 // 4))]
    p = sum(interval(c...) * x^i * y^j for ((i, j), c) in aterms)
    q = sum(interval(c...) * x^i * y^j for ((i, j), c) in bterms)
    powers = ext.power_cache([interval(-1, 1), interval(-1, 1)], 4, Interval{Float64})
    bounds = ext.unit_product_bounds(powers, 2, Interval{Float64})
    expected = rational_discarded_product(aterms, bterms, ((-1, 1), (-1, 1)), 2)
    actual = ext.discarded_unit_bound(p, q, 2, bounds, Interval{Float64})
    @test sameinterval(actual, interval(Float64, expected...))
    @test sameinterval(
        actual, ext.decorated_discarded_unit_bound(p, q, 2, bounds, Interval{Float64})
    )
end

@testset "Unsafe interval rounding fails before evaluation" begin
    x, = taylor_models([interval(-1, 1)]; order = 2)
    IA.configure(; rounding = :none)
    try
        @test_throws ArgumentError enclose(x)
        @test_throws ArgumentError coefficient_norm(polynomial(x), 2)
    finally
        IA.configure(; rounding = :correct)
    end
end
