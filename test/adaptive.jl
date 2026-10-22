using Test, DifferentialAlgebra, IntervalArithmetic, ForwardDiff
const DA = DifferentialAlgebra
const IA = IntervalArithmetic
quadratic(v) = v[1]^2 + v[2]
exponential(v) = exp(v[1])
cubic(v) = v[1]^3

@testset "Shared geometry contracts and normalized estimators" begin
    root, = variables(1; order = 3)
    square = ConvexPolygon([(-1, -1), (1, -1), (1, 1), (-1, 1)])
    @test domain_area(square) == 4
    @test_throws ArgumentError ConvexPolygon([(0, 0), (1, 0), (2, 0)])
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
        @test all(fit(collect(p)) ≈ p[1]^2 + p[2] for p in ((-1, -1), (0, 0), (1, 1)))
        @test_throws DomainError fit([2, 0])
        @test nvariables() == 1 && variable(1).algebra === root.algebra
    end
    options = (; order = 1, guard_order = 1, atol = 1, check_points = false)
    guarded = (
        adaptive_map(exponential, [-1, -1], [1, 1]; options...),
        adaptive_map(exponential, square; options...),
    )
    @test [only(f.patches).error_estimate[1] for f in guarded] ≈ [3 / 4, 3 / 4]
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
            x in (-1, -1 // 2, 0, 1 // 2, 1), y in (-1, 0, 1)
    )
    @test sum(domain_area(domain(p)) for p in refined.patches) == 4
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
    polygon = adaptive_map(
        quadratic, square; order = 1, atol = 1 // 4,
        estimator = IntervalBound(), directions = :axes
    )
    smooth = continuous_map(quadratic, polygon)
    @test IA.in_interval(1 // 4, polygon([1 // 2, 0]))
    @test IA.in_interval(1 // 4, enclose(smooth, [interval(1 // 2), interval(0)]))
end

@testset "Continuity, partition weights and irrational endpoints" begin
    core = adaptive_map(cubic, [-1.0], [1.0]; order = 2, atol = 1 // 2)
    smooth = continuous_map(cubic, core; continuity = :c2)
    @test length(core.patches) == 2
    for point in (0.0, 0.125, 0.25)
        @test sum(blend_weights(smooth, [point])) ≈ 1
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

@testset "Adaptive flow composes epochs and the initial map" begin
    advance(u, times) = u .* exp(last(times) - first(times))
    flow = adaptive_flow(
        advance, v -> [1 + v[1]], [-0.1], [0.1], [0.0, 0.2, 0.4];
        order = 2, atol = 1.0e-8
    )
    @test only(flow([0.05])) ≈ 1.05exp(0.4)
end
