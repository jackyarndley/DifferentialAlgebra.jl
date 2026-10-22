using Test, DifferentialAlgebra
using ForwardDiff, IntervalArithmetic, LinearAlgebra

include("support/processes.jl")

include("support/intervals.jl")
using .IntervalTestSupport: IA, DA, subset, interval_contains, sameinterval

exact_derivative(f, x, n) = n == 0 ? f(x) : ForwardDiff.derivative(t -> exact_derivative(f, t, n - 1), x)
const RContinuity = Rational{BigInt}

@testset "Irrational boundary routing agrees across ordinary geometries" begin
    f(v) = v[1]^2
    for upper in (Float64(π), nextfloat(Float64(π)))
        box = adaptive_map(f, [0, -1], [upper, 1]; order = 2)
        polygon = ConvexPolygon([(0, -1), (upper, -1), (upper, 1), (0, 1)])
        oriented = adaptive_map(f, polygon; directions = [1 1; -1 1], order = 2)
        for fit in (box, oriented), query in (fit, continuous_map(f, fit))
            if upper < BigFloat(π)
                @test_throws DomainError query(Real[π, 0])
            else
                @test query(Real[π, 0]) ≈ Float64(π)^2 rtol = 16eps()
            end
            @test query([BigFloat(upper), 0]) ≈ upper^2 rtol = 16eps()
        end
    end
end

@testset "Irrational physical coordinates" begin
    f(v) = v[1]^2 + v[2] / 2
    box = adaptive_map(f, [-4, -4], [4, 4]; order = 2)
    polygon = ConvexPolygon([(-4, -4), (4, -4), (4, 4), (-4, 4)])
    oriented = adaptive_map(f, polygon; directions = [1 1; -1 1], order = 2)
    for fit in (box, oriented)
        smooth = continuous_map(f, fit)
        for point in ((π, 0), [π, 0])
            @test smooth(point) ≈ f(point)
        end
    end
end

@testset "Compact tapers and exact one-sided continuity oracles" begin
    # Independent power-basis coefficients; endpoint jets prove the regularity,
    # rather than relying on finite differences of plotted samples.
    coefficients = ([1, -1], [1, 0, -3, 2], [1, 0, 0, -10, 15, -6])
    for k in 0:2
        c = coefficients[k + 1]
        oracle(t) = sum(c[i] * t^(i - 1) for i in eachindex(c))
        for t in RContinuity[0, 1 // 4, 1 // 2, 3 // 4, 1]
            @test DifferentialAlgebra.continuity_taper(t, k) == oracle(t)
            @test 0 <= oracle(t) <= 1
        end
        for d in 1:k, x in RContinuity[0, 1]
            jet = sum(c[i] * prod((i - d):(i - 1)) * x^(i - 1 - d) for i in (d + 1):length(c))
            @test jet == 0
        end
        # f=x³, degree-two local fits on [-1,1/4] and [-1/4,1].
        # Taylor's formula gives P_c(x)=3cx²-3c²x+c³ exactly.
        polynomial(c, x) = 3c * x^2 - 3c^2 * x + c^3
        left(x) = polynomial(-3 // 8, x)
        right(x) = polynomial(3 // 8, x)
        inner(x) = (oracle(4x) * left(x) + right(x)) / (1 + oracle(4x))
        face = RContinuity(1 // 4)
        for d in 0:k
            @test exact_derivative(inner, face, d) == exact_derivative(right, face, d)
        end
        @test exact_derivative(inner, face, k + 1) != exact_derivative(right, face, k + 1)
        core = adaptive_map(v -> v[1]^3, [-1.0], [1.0]; order = 2, atol = 1 // 2)
        @test length(core.patches) == 2
        smooth = continuous_map(v -> v[1]^3, core; continuity = (:c0, :c1, :c2)[k + 1])
        for x in (0.0, 0.125, 0.25), d in 0:k
            expected = exact_derivative(x < 0.25 ? inner : right, RContinuity(x), d)
            @test exact_derivative(t -> smooth([t]), x, d) ≈ Float64(expected) atol = 2.0e-13
        end
    end
end

@testset "Owned continuous ADS, geometry and AD" begin
    for T in (Float64, BigFloat), budget in (0, 32 * 1024^2)
        f(v) = [v[1]^3 + 2v[2]^3, 2 + v[1] / 3 - v[2] / 5]
        caller, = variables(1; order = 4)
        ctx = DifferentialAlgebra.CURRENT_ALGEBRA[]
        for oriented in (false, true)
            fit = adaptive_map(f, T[-1, -1], T[1, 1]; order = 2, atol = 1 // 2, splitter = oriented ? :oriented : :tail, table_bytes = budget)
            s = continuous_map(f, fit; table_bytes = budget)
            @test DifferentialAlgebra.CURRENT_ALGEBRA[] === ctx && degree(caller) == 1
            @test s.continuity == :c2 && nvariables(s) == 2 && noutputs(s) == 2 && max_order(s) == 2
            @test all(isfinite, s.error_estimate)
            @test_throws ArgumentError s.error_bounds
            @test_throws ArgumentError degree(s)
            @test_throws ArgumentError enclose(s)
            for x in (T[0, 0], T[1, 1], T[-1, -1], T[1 // 3, -1 // 7])
                w = blend_weights(s, x)
                @test all(w .>= 0) && sum(w) ≈ 1
                # Affine reproduction is exact in real arithmetic. This checks
                # normalization/mapping at faces, intersections and corners.
                @test s(x)[2] ≈ f(x)[2] atol = 20eps(T)
                @test all(isfinite, ForwardDiff.jacobian(s, x))
            end
            x = T[0, 0]
            J = ForwardDiff.jacobian(s, x)
            @test J[2, :] ≈ T[1 // 3, -1 // 5] atol = 30eps(T)
            hessian = ForwardDiff.hessian(z -> s(z)[2], x)
            @test maximum(abs, hessian) <= 200eps(T)
            @test_throws DimensionMismatch s([0])
            @test_throws DimensionMismatch blend_weights(s, [0, 0, 0])
            @test_throws DomainError s([2, 0])
            @test_throws ArgumentError s([interval(0), interval(0)])
            before = s(T[1 // 4, 1 // 4])
            b = copy(s)
            @test b._patches[1].map.coefficients !== s._patches[1].map.coefficients
            if !oriented
                d = domain(s); d.lower[1] = 100
                fit.patches[1].map.coefficients .= 100
            end
            initialize!(2, 3)
            @test s(T[1 // 4, 1 // 4]) == before == b(T[1 // 4, 1 // 4])
            # Restore the caller before the next independent construction.
            caller, = variables(1; order = 4); ctx = DifferentialAlgebra.CURRENT_ALGEBRA[]
        end
        fixed = adaptive_map(v -> v[1] + v[2], T[-1, 2], T[1, 2]; order = 1)
        fixed_s = continuous_map(v -> v[1] + v[2], fixed; order = 1)
        @test fixed_s(T[1 // 4, 2]) ≈ 9 // 4
        @test_throws DomainError fixed_s([0, 3])
    end
    core = adaptive_map(v -> v[1]^2, [-1.0], [1.0]; order = 1, atol = 1 // 4)
    @test_throws ArgumentError continuous_map(identity, core; continuity = :c3)
    @test_throws ArgumentError continuous_map(identity, core; overlap = 0)
    @test_throws ArgumentError continuous_map(identity, core; overlap = Inf)
    @test_throws ArgumentError continuous_map(v -> v[1]^2, core; overlap = 1 // big(10)^1000)
    @test_throws ArgumentError continuous_map(identity, core; estimator = IntervalBound())
    @test_throws ArgumentError continuous_map(v -> (set_coefficient_tolerance!(1); v[1]), core)
    @test_throws ArgumentError continuous_map(v -> (set_truncation_order!(1); v[1]), core)
    @test_throws DimensionMismatch continuous_map(v -> [v[1], v[1]], core)
    @test_throws ArgumentError continuous_map(v -> v[1]^2, core; atol = 1 // 1000)
    script = raw"""
    using DifferentialAlgebra
    @assert !DifferentialAlgebra.isinitialized()
    @assert !any(m -> nameof(m) in (:IntervalArithmetic, :ForwardDiff),
        values(Base.loaded_modules))
    using ForwardDiff
    @assert !DifferentialAlgebra.isinitialized()
    a = adaptive_map(v -> v[1]^2, [-1.0], [1.0]; order = 2)
    s = continuous_map(v -> v[1]^2, a)
    @assert ForwardDiff.derivative(x -> s([x]), 0.25) ≈ 0.5
    @assert !DifferentialAlgebra.isinitialized()
    """
    @test fresh_julia(script)
end

@testset "Certified function bounds for continuous surrogates" begin
    for T in (Float64, BigFloat), budget in (0, 32 * 1024^2)
        f(v) = v[1]^2
        fit = adaptive_map(f, [interval(T, -1, 1)]; estimator = IntervalBound(), order = 1, atol = 1 // 16, table_bytes = budget)
        s = continuous_map(f, fit; table_bytes = budget)
        # Independent uniform oracle f-P_c=(x-c)² on the overlap interval.
        # The largest radius is 3/8; convex blending preserves [0,9/64].
        @test isequal_interval(only(s.error_bounds), interval(T, 0, 9 // 64))
        @test isguaranteed(only(s.error_bounds)) && IntervalArithmetic.decoration(only(s.error_bounds)) == IntervalArithmetic.com
        @test_throws ArgumentError continuous_map(f, fit; atol = 1 // 16)
        @test continuous_map(f, fit; atol = 9 // 64) isa ContinuousTaylorMap
        for x in RContinuity[-1, -1 // 2, -1 // 4, 0, 1 // 4, 1 // 2, 1]
            @test in_interval(x^2, enclose(s, [x]))
            @test isguaranteed(enclose(s, [x])) && IntervalArithmetic.decoration(enclose(s, [x])) == IntervalArithmetic.com
        end
        @test issubset_interval(interval(T, 0, 1), enclose(s))
        @test issubset_interval(interval(T, 0, 1 // 4), enclose(s, [interval(T, -1 // 2, 1 // 2)]))
        @test_throws DimensionMismatch enclose(s, [])
        @test_throws DimensionMismatch enclose(s, ConvexPolygon([(0, 0), (1, 0), (0, 1)]))
        @test_throws DomainError enclose(s, [interval(T, 0, 2)])
        @test_throws DomainError enclose(s, [2])
        @test_throws ArgumentError enclose(s, [convert(Interval{T}, zero(T))])
        @test_throws ArgumentError enclose(s, [bareinterval(-1, 1)])
        @test_throws DomainError continuous_map(v -> log(v[1]), fit)
        @test_throws ArgumentError continuous_map(v -> convert(Interval{T}, zero(T)), fit)
        before = enclose(s)
        b = copy(s)
        @test b._patches[1].certificates[1]._coefficients !== s._patches[1].certificates[1]._coefficients
        initialize!(2, 3)
        @test isequal_interval(enclose(s), before) && isequal_interval(enclose(b), before)
        triangle = ConvexPolygon([(zero(T), zero(T)), (one(T), zero(T)), (zero(T), one(T))])
        p = adaptive_map(v -> (v[1] + v[2])^2, triangle; estimator = IntervalBound(), directions = [1 1; -1 1], order = 1, atol = 1 // 64, table_bytes = budget)
        t = continuous_map(v -> (v[1] + v[2])^2, p; table_bytes = budget)
        @test in_interval(1 // 4, enclose(t, [1 // 4, 1 // 4]))
        @test issubset_interval(interval(T, 0, 1), enclose(t, triangle))
        @test_throws DomainError t([3 // 4, 3 // 4])
        @test_throws DomainError enclose(t, [interval(T, 0, 1), interval(T, 0, 1)])
        @test_throws DomainError enclose(t, ConvexPolygon([(0, 0), (2, 0), (0, 2)]))
        @test all(isfinite, ForwardDiff.hessian(t, T[1 // 4, 1 // 4]))
        six = adaptive_map(v -> sum(v), fill(interval(T, -1, 1), 6); estimator = IntervalBound(), order = 1)
        six_s = continuous_map(v -> sum(v), six)
        @test in_interval(2, enclose(six_s, fill(1 // 3, 6)))
        @test ForwardDiff.gradient(six_s, zeros(T, 6)) == ones(T, 6)
    end
end

# Full scans use independent power-basis tapers and exact geometry. They do not
# use the support index or its floating filters.
function scan_continuity_weights(a, point)
    T = promote_type(eltype(first(a._patches).center), eltype(point))
    z = T.(a._directions) * point
    x = RContinuity.(DA.continuity_primal.(point))
    exact = [sum(a._directions[i, j] * x[j] for j in eachindex(x)) for i in eachindex(x)]
    c = a.continuity == :c0 ? [1, -1] :
        a.continuity == :c1 ? [1, 0, -3, 2] : [1, 0, 0, -10, 15, -6]
    weights = map(a._patches) do p
        weight = one(first(z))
        for i in eachindex(z)
            p.core_lower[i] <= exact[i] <= p.core_upper[i] && continue
            p.lower[i] < exact[i] < p.upper[i] || return zero(weight)
            t = if exact[i] < p.core_lower[i]
                (T(p.core_lower[i]) - z[i]) / T(p.core_lower[i] - p.lower[i])
            else
                (z[i] - T(p.core_upper[i])) / T(p.upper[i] - p.core_upper[i])
            end
            t = clamp(t, zero(t), one(t))
            weight *= sum(c[k] * t^(k - 1) for k in eachindex(c))
        end
        weight
    end
    return z, weights ./ sum(weights)
end
function scan_continuity(a, point)
    z, weights = scan_continuity_weights(a, point)
    R = promote_type(eltype(first(a._patches).map.coefficients), eltype(z))
    value = zeros(R, noutputs(a))
    for (p, weight) in zip(a._patches, weights)
        iszero(DA.continuity_primal(weight)) && continue
        normalized = [
            iszero(p.radius[i]) ? zero(z[i]) : (z[i] - p.center[i]) / p.radius[i]
                for i in eachindex(z)
        ]
        value .+= weight .* evaluate(p.map, normalized)
    end
    return a._scalar ? only(value) : value
end
function scan_continuity_enclosure(a, query)
    result = nothing
    for p in a._patches
        lower = max.(RContinuity.(IA.inf.(query)), p.lower)
        upper = min.(RContinuity.(IA.sup.(query)), p.upper)
        any(lower .> upper) && continue
        T = IA.numtype(typeof(first(p.certificates)._remainder))
        box = [
            IA.intersect_interval(interval(T, lo, hi), b; dec = :auto)
                for (lo, hi, b) in zip(lower, upper, first(p.certificates)._coordinates.box)
        ]
        values = [enclose(m, box) for m in p.certificates]
        result = result === nothing ? values : IA.hull.(result, values; dec = :auto)
    end
    return a._scalar ? only(result) : result
end

@testset "Conservative floating projection filters" begin
    bigtwo = big(2)^1100
    matrices = (
        RContinuity[1 0; 0 1], RContinuity[1 // 3 -7 // 9; 5 // 11 2 // 17],
        RContinuity[bigtwo 1; -bigtwo -1],
        RContinuity[1 // bigtwo -1 // bigtwo; 1 -1],
    )
    points = (
        [0.0, -0.0], [prevfloat(1.0), nextfloat(1.0)],
        [nextfloat(0.0), -nextfloat(0.0)], [floatmax(Float64), -floatmax(Float64)],
        RContinuity[1 // 3, 1 // 7], BigFloat[1 // 3, 1 // 7],
        RContinuity[bigtwo, bigtwo + 1],
    )
    for B in matrices, point in points
        lower, upper = DA.ads_project_bounds(
            DA.ads_float_bounds(B),
            DA.ads_float_bounds(point)
        )
        exact = B * RContinuity.(point)
        for i in eachindex(exact)
            @test lower[i] == -Inf || RContinuity(lower[i]) <= exact[i]
            @test upper[i] == Inf || exact[i] <= RContinuity(upper[i])
        end
    end
    polygons = (
        ConvexPolygon([(-1, -1), (1, -1), (1, 1), (-1, 1)]),
        ConvexPolygon([(0, 0), (1, 1 // 3), (0, 1)]),
        ConvexPolygon([(1, 1), (2, 1), (2, 1 + 1 // big(10)^30)]),
    )
    for polygon in polygons
        planes = DA.polygon_filter(polygon)
        vertices = polygon_vertices(polygon)
        for x in vcat(vertices, [(1 // 3, 1 // 3), (0, 0), (2, 2)])
            exact = DA.polygon_contains(polygon, RContinuity.(x))
            @test DA.polygon_contains_filtered(
                polygon, x,
                DA.ads_float_bounds(x), planes
            ) == exact
        end
        for vertex in vertices, direction in (-1, 0, 1)
            x = Float64.(vertex)
            shifted = (
                direction < 0 ? prevfloat(x[1]) :
                    direction > 0 ? nextfloat(x[1]) : x[1], x[2],
            )
            exact = DA.polygon_contains(polygon, RContinuity.(shifted))
            @test DA.polygon_contains_filtered(
                polygon, shifted,
                DA.ads_float_bounds(shifted), planes
            ) == exact
        end
    end
end

@testset "Indexed supports agree with independent full scans" begin
    f(v) = v[1]^2 + 2v[2]^2
    polygon = ConvexPolygon([(-1, -1), (1, -1), (1, 1), (-1, 1)])
    for budget in (0, 32 * 1024^2), oriented in (false, true)
        fit = if oriented
            adaptive_map(
                f, polygon; directions = [1 1 // 3; -1 // 3 1],
                order = 1, atol = 1 // 8, table_bytes = budget
            )
        else
            adaptive_map(
                f, [-1.0, -1.0], [1.0, 1.0];
                order = 1, atol = 1 // 8, table_bytes = budget
            )
        end
        for continuity in (:c0, :c1, :c2)
            smooth = continuous_map(f, fit; continuity, table_bytes = budget)
            @test length(smooth._nodes) == 2length(smooth._patches) - 1
            for point in (
                    [0.123, 0.234], [0.0, 0.0], [-1.0, 1.0],
                    [prevfloat(0.25), nextfloat(0.25)], [0.25, 0.25],
                )
                _, expected = scan_continuity_weights(smooth, point)
                @test blend_weights(smooth, point) ≈ expected atol = 2.0e-14
                @test smooth(point) ≈ scan_continuity(smooth, point) atol = 2.0e-14
                exact = smooth._directions * RContinuity.(point)
                matching = findall(smooth._patches) do p
                    all(p.lower .<= exact .<= p.upper)
                end
                _, bounds = DA.continuity_query(smooth, point)
                indices = DA.continuity_indices(smooth, bounds)
                @test issubset(matching, indices)
            end
            if continuity == :c2
                point = [0.123, 0.234]
                @test ForwardDiff.gradient(smooth, point) ≈
                    ForwardDiff.gradient(x -> scan_continuity(smooth, x), point)
                @test ForwardDiff.hessian(smooth, point) ≈
                    ForwardDiff.hessian(x -> scan_continuity(smooth, x), point)
            end
            owned = copy(smooth)
            @test owned._nodes !== smooth._nodes
            @test owned._float_directions !== smooth._float_directions
            @test owned._patches[1].bounds[1] !== smooth._patches[1].bounds[1]
            @test owned._patches[1].tapers !== smooth._patches[1].tapers
            @test_throws DomainError smooth([nextfloat(1.0), 0.0])
        end
    end
    # A local query must not visit the whole partition. This verifies the index
    # independently of timing, and blend_weights still returns one entry per fit.
    fit = adaptive_map(f, [-1.0, -1.0], [1.0, 1.0]; order = 1, atol = 1 // 64)
    smooth = continuous_map(f, fit)
    _, bounds = DA.continuity_query(smooth, [0.123, 0.234])
    @test length(DA.continuity_indices(smooth, bounds)) < length(smooth._patches) ÷ 8
    @test length(blend_weights(smooth, [0.123, 0.234])) == length(smooth._patches)
end

@testset "Exact fallback at faces and beyond Float64 precision" begin
    for budget in (0, 32 * 1024^2)
        f(v) = v[1]^3
        fit = adaptive_map(
            f, [-1.0], [1.0]; order = 2, atol = 1 // 2,
            table_bytes = budget
        )
        smooth = continuous_map(f, fit; table_bytes = budget)
        for face in unique(
                vcat(
                    [
                        vcat(p.lower, p.core_lower, p.core_upper, p.upper)
                            for p in smooth._patches
                    ]...
                )
            )
            for x in (prevfloat(Float64(face)), Float64(face), nextfloat(Float64(face)))
                -1 <= x <= 1 || continue
                @test blend_weights(smooth, [x]) ≈
                    last(scan_continuity_weights(smooth, [x]))
                @test smooth([x]) ≈ scan_continuity(smooth, [x]) atol = 2.0e-15
            end
        end
        fixed = adaptive_map(
            v -> v[1]^2 + v[2], BigFloat[-1, 1 // 3],
            BigFloat[1, 1 // 3]; order = 1, atol = 1 // 8, table_bytes = budget
        )
        fixed_smooth = continuous_map(v -> v[1]^2 + v[2], fixed; table_bytes = budget)
        point = BigFloat[1 // 7, 1 // 3]
        @test fixed_smooth(point) ≈ scan_continuity(fixed_smooth, point)
        @test_throws DomainError fixed_smooth([point[1], nextfloat(point[2])])
        # Build in one precision, query in another: cached Float64 directions
        # and tapers must never determine BigFloat or nested-dual arithmetic.
        polygon = ConvexPolygon([(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)])
        f2(v) = v[1]^2 + v[2]^2
        source = adaptive_map(
            f2, polygon; directions = [1 1 // 3; -1 // 3 1],
            order = 1, atol = 1 // 8, table_bytes = budget
        )
        high = continuous_map(f2, source; table_bytes = budget)
        setprecision(512) do
            point = BigFloat[1 // 7, 2 // 9]
            @test high(point) isa BigFloat
            @test high(point) ≈ scan_continuity(high, point) atol = 200eps(BigFloat)
            @test ForwardDiff.gradient(high, point) ≈
                ForwardDiff.gradient(x -> scan_continuity(high, x), point)
            @test ForwardDiff.hessian(high, point) ≈
                ForwardDiff.hessian(x -> scan_continuity(high, x), point)
        end
    end
end

@testset "Indexed interval enclosures retain closed support intersections" begin
    for T in (Float64, BigFloat), budget in (0, 32 * 1024^2)
        f(v) = v[1]^2
        box = [interval(T, -1, 1)]
        fit = adaptive_map(
            f, box; estimator = IntervalBound(), order = 1,
            atol = 1 // 64, table_bytes = budget
        )
        smooth = continuous_map(f, fit; table_bytes = budget)
        queries = [[interval(T, x)] for x in (-1, -1 // 4, 0, 1 // 4, 1)]
        append!(queries, ([interval(T, -1 // 2, 1 // 2)], box))
        for query in queries
            @test sameinterval(
                enclose(smooth, query),
                scan_continuity_enclosure(smooth, query)
            )
            @test isguaranteed(enclose(smooth, query))
        end
    end
end
