using Test, DifferentialAlgebra, IntervalArithmetic, ForwardDiff
using Main.Fixtures: ShiftedVector
const DA = DifferentialAlgebra
const IA = IntervalArithmetic
quadratic(v) = v[1]^2 + v[2]
exponential(v) = exp(v[1])
cubic(v) = v[1]^3
tail_polynomial(v) = v[1]^2 + 3v[1]^3

@testset "Shared geometry contracts and normalized estimators" begin
    root, = variables(1; order = 3)
    square = ConvexPolygon([(-1, -1), (1, -1), (1, 1), (-1, 1)])
    @test domain_area(square) == 4
    @test_throws ArgumentError ConvexPolygon([(0, 0), (1, 0), (2, 0)])
    @test_throws ArgumentError ConvexPolygon(ShiftedVector([(0, 0), (1, 0), (0, 1)]))
    @test_throws ArgumentError ConvexPolygon([ShiftedVector([0, 0]), [1, 0], [0, 1]])
    @test_throws ArgumentError adaptive_map(v -> ShiftedVector([v[1]]), [0.0], [1.0])
    calls = Ref(0)
    callback(v) = (calls[] += 1; quadratic(v))
    @test_throws ArgumentError adaptive_map(callback, [-1, -1], [1, 1]; atol = -1)
    @test iszero(calls[])
    fits = (
        adaptive_map(quadratic, [-1, -1], [1, 1]; order = 2),
        adaptive_map(quadratic, square; order = 2, directions = :axes),
        adaptive_map(quadratic, square; order = 2),
    )
    for fit in fits
        @test fit.converged && length(fit.patches) == 1
        @test all(fit(p) ≈ p[1]^2 + p[2] for p in ((-1, -1), (0, 0), (1, 1)))
        @test_throws DomainError fit([2, 0])
        @test nvariables() == 1 && variable(1).algebra === root.algebra
    end
    @test ForwardDiff.derivative(t -> first(fits)((t, 0)), 1.0) == 2
    # Dual seeds must not move a physical point to the upper child. At zero,
    # the lower patch's quadratic expansion of x⁴ about -1/2 has derivative 1.
    quartic = adaptive_map(v -> v[1]^4, [-1.0], [1.0]; order = 2, atol = 0.6)
    @test ForwardDiff.derivative(t -> quartic((t,)), 0.0) == 1
    # Degree norms 1 and 3 give an extrapolated next norm of 9.
    for (estimator, order, guard_order, expected) in (
            (GuardedTail(), 2, 1, 12), (ExtrapolatedTail(), 3, 0, 9),
            (LastTerms(; degrees = 1, safety = 2), 3, 0, 6),
        )
        fit = adaptive_map(
            tail_polynomial, [-1.0], [1.0];
            estimator, order, guard_order, atol = 100, check_points = false
        )
        @test only(only(fit.patches).error_estimate) ≈ expected
    end
    source = adaptive_map(
        exponential, [-1.0], [1.0]; order = 2,
        atol = 1.0e-4, max_depth = 0, strict = false
    )
    @test !source.converged && only(source.patches).status == :max_depth
    refined = adaptive_map(exponential, source; max_depth = 4, atol = 0.01)
    @test refined.converged && length(source.patches) == 1
    @test refined([0.3]) ≈ exp(0.3) atol = 0.01
    # Exercise the polygon search tree and exact closed faces with multiple leaves.
    refined = adaptive_map(quadratic, square; order = 1, atol = 1 // 4, directions = :axes)
    @test length(refined.patches) > 1
    @test all(
        abs(refined([x, y]) - (x^2 + y)) <= 1 // 4 for
            (x, y) in ((-1, -1), (-1 // 2, 0), (0, 1), (1 // 2, -1), (1, 1))
    )
    @test sum(domain_area(domain(p)) for p in refined.patches) == 4
end

@testset "Frame scaling and representable physical covers" begin
    lo, hi = 2.0 - 1.0e-6, 2.0 + 1.0e-6
    polygon = ConvexPolygon([(lo, lo), (hi, lo), (hi, hi), (lo, hi)])
    for scale in (1.0e308, 1 // big(2)^1100)
        fit = adaptive_map(
            first, polygon; directions = [scale -scale; 0 1],
            order = 1, check_points = false
        )
        @test split_directions(fit) == [1 -1; 0 1]
        @test domain_area(domain(only(fit.patches))) == domain_area(polygon)
        @test fit((2.0, 2.0)) ≈ 2.0
        @test continuous_map(first, fit; check_points = false)((2.0, 2.0)) ≈ 2.0
        @test ForwardDiff.derivative(t -> fit((t, 2.0)), 2.0) ≈ 1
    end
    # Row normalization cannot rescue a nonrepresentable physical cover from
    # nearly parallel directions. Fail before evaluating the callback.
    calls = Ref(0)
    @test_throws ArgumentError adaptive_map(
        v -> (calls[] += 1; first(v)), polygon;
        directions = [1 0; 1 1 // big(2)^1100], order = 1
    )
    @test iszero(calls[])
    @test_throws ArgumentError adaptive_map(first, polygon; directions = [1 1; 2 2])
end

@testset "Certified subdivision and closed-boundary routing" begin
    fit = adaptive_map(
        v -> v[1]^2, [interval(-1, 1)]; order = 1,
        atol = 1 // 16, estimator = IntervalBound(), splitter = :width
    )
    @test fit.converged && length(fit.patches) == 4
    @test all(
        IA.isequal_interval(only(p.error_bounds), interval(0, 1 // 16)) for p in fit.patches
    )
    @test all(IA.in_interval(p^2, fit([p])) for p in (-1, -1 // 2, 0, 1 // 2, 1))
    @test IA.issubset_interval(
        interval(0, 1 // 4), enclose(fit, [interval(-1 // 2, 1 // 2)])
    )
    saved = enclose(fit)
    box = domain(fit)
    box[1] = interval(100)
    initialize!(1, 1)
    @test IA.isequal_interval(enclose(fit), saved)
    square = ConvexPolygon([(-1, -1), (1, -1), (1, 1), (-1, 1)])
    # Every output must include extrema from different leaves of the partition.
    vector_map(v) = [v[1]^2, v[1] + v[2]]
    box_fit = adaptive_map(
        vector_map, [interval(-1, 1), interval(-1, 1)]; order = 1, atol = 1 // 4
    )
    polygon_fit = adaptive_map(
        vector_map, square; order = 1, atol = 1 // 4,
        estimator = IntervalBound(), directions = :axes
    )
    for fit in (box_fit, polygon_fit, continuous_map(vector_map, polygon_fit))
        bounds = enclose(fit)
        @test length(bounds) == 2 && all(IA.isguaranteed, bounds)
        @test IA.issubset_interval(interval(0, 1), bounds[1])
        @test IA.issubset_interval(interval(-2, 2), bounds[2])
        point = [interval(1 // 2), IA.setdecoration(interval(0), IA.def)]
        values = enclose(fit, point)
        @test all(IA.in_interval.([1 // 4, 1 // 2], values))
        @test all(x -> IA.decoration(x) == IA.def && IA.isguaranteed(x), values)
        @test_throws ArgumentError enclose(fit, [convert(Interval{Float64}, 0.5), interval(0)])
    end
    # Exact projection must retain metadata even on an unused coordinate, and
    # subsequent refinement/smoothing must carry it into their fresh models.
    decorated = [interval(-1, 1), IA.setdecoration(interval(-1, 1), IA.def)]
    with_constant(v) = [v[1]^2, 1]
    polygon = adaptive_map(
        with_constant, decorated; order = 2, atol = 1 // 4,
        splitter = :oriented, directions = :axes
    )
    refined = adaptive_map(with_constant, polygon; order = 1)
    source = adaptive_map(with_constant, decorated; order = 2, atol = 1 // 4)
    for fit in (polygon, refined, continuous_map(with_constant, refined), continuous_map(with_constant, source))
        bounds = enclose(fit)
        @test all(x -> IA.decoration(x) == IA.def && IA.isguaranteed(x), bounds)
        @test IA.issubset_interval(interval(0, 1), bounds[1]) && IA.in_interval(1, bounds[2])
    end
    @test_throws ArgumentError adaptive_map(
        with_constant, [interval(-1, 1), convert(Interval{Float64}, 0.0)];
        splitter = :oriented, directions = :axes
    )
    # strict=false retains valid unresolved models, never invalid function domains.
    @test_throws ErrorException adaptive_map(
        v -> log(v[1]), [interval(-1, 1)];
        order = 2, max_depth = 0, strict = false
    )
end

@testset "Certified irrational queries preserve exact domain membership" begin
    reference = setprecision(BigFloat, 512) do
        BigFloat(pi)
    end
    setprecision(BigFloat, 192) do
        broad = adaptive_map(first, [interval(BigFloat, -4, 4)]; order = 1)
        smooth = continuous_map(first, broad)
        @test IA.in_interval(reference, enclose(smooth, (pi,)))
        narrow = adaptive_map(first, [interval(BigFloat, pi)]; order = 1)
        smooth = continuous_map(first, narrow)
        setprecision(BigFloat, 64) do
            @test IA.in_interval(reference, narrow((pi,)))
            @test IA.in_interval(reference, enclose(smooth, Real[pi]))
            @test_throws DomainError enclose(smooth, (-pi,))
        end
    end
    # The narrowest Float64 bracket of π is admitted as a whole, including
    # when transformed into polygon coordinates; the point -π is outside.
    lower = Float64(pi)
    polygon = adaptive_map(
        first, [interval(lower, nextfloat(lower)), interval(-1, 1)];
        order = 1, splitter = :oriented, directions = [1 1; -1 1]
    )
    @test IA.in_interval(reference, polygon((pi, 0)))
    @test IA.in_interval(reference, polygon(Real[pi, 0]))
    @test_throws DomainError polygon((-pi, 0))
    @test IA.in_interval(reference, enclose(continuous_map(first, polygon), (pi, 0)))
    # A rounded-down Float64 endpoint is strictly below the exact scalar π.
    outside = adaptive_map(first, [interval(0, lower)]; order = 1)
    @test_throws DomainError outside((pi,))
    @test_throws DomainError enclose(continuous_map(first, outside), (pi,))
    outside = adaptive_map(
        first, [interval(0, lower), interval(-1, 1)];
        order = 1, splitter = :oriented, directions = :axes
    )
    @test_throws DomainError outside((pi, 0))
end

@testset "Continuity, partition weights and irrational endpoints" begin
    core = adaptive_map(cubic, [-1.0], [1.0]; order = 2, atol = 1 // 2)
    smooth = continuous_map(cubic, core; continuity = :c2)
    @test length(core.patches) == 2
    for (point, expected) in ((0.0, [1 // 2, 1 // 2]), (0.125, [1 // 3, 2 // 3]), (0.25, [0, 1]))
        @test blend_weights(smooth, [point]) ≈ expected
    end
    # At the support endpoint the C2 blend has the right patch's exact jets.
    @test ForwardDiff.derivative(t -> smooth([t]), 0.25) ≈ 3 * (3 / 8) / 2 - 3 * (3 / 8)^2
    @test ForwardDiff.derivative(
        t -> ForwardDiff.derivative(s -> smooth([s]), t), 0.25
    ) ≈ 6 * 3 / 8
    endpoint = Float64(pi)
    fit = adaptive_map(quadratic, [endpoint, -1], [nextfloat(endpoint), 1]; order = 2)
    @test fit(Real[pi, 0]) ≈ Float64(BigFloat(pi)^2)
    @test_throws DomainError fit([prevfloat(endpoint), 0])
end

@testset "Adaptive flow refines and always returns the final epoch" begin
    # Exact flow u'=u² makes spatial truncation independent of ODE integration.
    advance(u, times) = u / (1 - (last(times) - first(times)) * u)
    flow = adaptive_flow(
        advance, first, [0.0], [0.6], [0.0, 0.5, 1.0]; order = 3, atol = 1.0e-3
    )
    @test flow.converged && length(flow.patches) > 1
    @test all(abs(flow([x]) - x / (1 - x)) <= 1.0e-3 for x in (0.1, 0.3, 0.55))
    unresolved = adaptive_flow(
        advance, first, [0.0], [0.6], [0.0, 0.5, 1.0];
        order = 3, atol = 1.0e-3, max_depth = 0, strict = false
    )
    @test !unresolved.converged && only(unresolved.patches).status == :max_depth
    @test unresolved([0.3]) ≈ 0.3 / 0.7
end
