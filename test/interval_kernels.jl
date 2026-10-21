using Test, DifferentialAlgebra, IntervalArithmetic, Random

include("support/interval_products.jl")

const DA = DifferentialAlgebra
const IA = IntervalArithmetic
const EXT = Base.get_extension(DA, :DifferentialAlgebraIntervalArithmeticExt)

function same_result(a, b)
    return IA.isequal_interval(a, b) && IA.decoration(a) == IA.decoration(b) &&
        IA.isguaranteed(a) == IA.isguaranteed(b)
end

@testset "Exact normalized coordinate powers" begin
    for T in (Float64, BigFloat), endpoints in ((-1, 1), (0, 0), (-1, 0), (1, 2)),
            decoration in (IA.com, IA.dac, IA.def), guaranteed in (false, true)
        x = IA.setdecoration(interval(T, endpoints...), decoration)
        guaranteed || (x *= T(1))
        for k in -2:12
            @test same_result(@inferred(EXT.coordinate_power(x, k)), IA.pown(x, k))
        end
    end
end

@testset "Unit-box kernel against exact rational bounds" begin
    rng = MersenneTwister(12431)
    I = Interval{Float64}
    for budget in (0, 32 * 1024^2), fixed in (false, true)
        x = variables(I, 3; order = 3, table_bytes = budget)
        basis = first(x).algebra.basis
        count = basis.ends[end]
        p, q = DA.allocate(first(x).algebra, I, count),
            DA.allocate(first(x).algebra, I, count)
        aterms, bterms = [], []
        # Dyadic inputs make this independent exact oracle representable in
        # Float64. Include all coefficient signs, widths and degrees.
        for (poly, terms) in ((p, aterms), (q, bterms))
            for (k, exponents) in enumerate(DA.multiindices(3, 3))
                endpoints = extrema(rand(rng, -8:8, 2) .// 8)
                poly.coeffs[k] = interval(Float64, endpoints...)
                push!(terms, (Tuple(Int.(exponents)), endpoints))
            end
            poly.len = count
        end
        box = [interval(-1, 1), interval(-1, 1), interval(fixed ? 0 : -1, fixed ? 0 : 1)]
        powers = EXT.power_cache(box, 6, I)
        bounds = EXT.unit_product_bounds(powers, 3, I)
        for order in 0:6
            exact_box = ((-1, 1), (-1, 1), fixed ? (0, 0) : (-1, 1))
            lo, hi = rational_discarded_product(aterms, bterms, exact_box, order)
            actual = @inferred EXT.discarded_unit_bound(p, q, order, bounds, I)
            @test same_result(actual, interval(Float64, lo, hi))
            @test same_result(
                actual,
                EXT.decorated_discarded_unit_bound(p, q, order, bounds, I)
            )
        end
    end
end

@testset "Kernel rounding extremes and metadata" begin
    I = Interval{Float64}
    extremes = (
        (interval(floatmax(Float64)), interval(2)),
        (interval(-floatmax(Float64)), interval(2)),
        (interval(-floatmax(Float64), floatmax(Float64)), interval(2)),
        (interval(floatmin(Float64)), interval(floatmin(Float64))),
        (interval(nextfloat(0.0)), interval(0.5)),
        (interval(-nextfloat(0.0)), interval(0.5)),
        (interval(-0.5, 0.25), interval(-0.25, 0.5)),
        (IA.setdecoration(interval(0), IA.def), interval(2)),
    )
    for budget in (0, 32 * 1024^2), fixed in (false, true),
            decoration in (IA.com, IA.def), guaranteed in (false, true)
        x, = variables(I, 1; order = 3, table_bytes = budget)
        axis = IA.setdecoration(interval(fixed ? 0 : -1, fixed ? 0 : 1), decoration)
        powers = EXT.power_cache([axis], 6, I)
        bounds = EXT.unit_product_bounds(powers, 1, I)
        for (ac, bc) in extremes, (degree_a, degree_b) in ((2, 2), (2, 3), (3, 3))
            guaranteed || (ac *= 1.0)
            p, q = zero(x), zero(x)
            set_coefficient!(p, [degree_a], ac)
            set_coefficient!(q, [degree_b], bc)
            @test EXT.finite_guaranteed_coefficients(p) == guaranteed
            actual = @inferred EXT.discarded_unit_bound(p, q, 3, bounds, I)
            expected = EXT.decorated_discarded_unit_bound(p, q, 3, bounds, I)
            @test same_result(actual, expected)
            exact_product = rational_interval_product(
                Rational{BigInt}.(IA.bounds(ac)), Rational{BigInt}.(IA.bounds(bc))
            )
            factor = fixed ? (0, 0) : iseven(degree_a + degree_b) ? (0, 1) : (-1, 1)
            lo, hi = rational_interval_product(exact_product, factor)
            @test IA.in_interval(lo, actual)
            @test IA.in_interval(hi, actual)
            if fixed && IA.isunbounded(ac * bc)
                @test IA.isthinzero(actual)
                @test IA.decoration(actual) == min(decoration, IA.dac)
            end
        end
    end
end
