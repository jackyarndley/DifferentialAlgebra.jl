@testset "Automatic domain splitting" begin
    old_x, old_y = variables(2; order = 3, names = (:u, :v))
    set_truncation_order!(2)
    set_coefficient_tolerance!(1.0e-14)
    old = 2 + old_x * old_y

    exact = adaptive_map(x -> [1 + x[1]^2, x[1] * x[2]], [-1, -2], [1, 2]; order = 2)
    @test exact.converged
    @test length(exact.patches) == 1
    @test exact([0.3, -0.4]) ≈ [1.09, -0.12]
    @test nvariables(exact) == 2
    @test noutputs(exact) == 2
    @test degree(exact) == 2
    @test truncation_order() == 2
    @test coefficient_tolerance() == 1.0e-14
    @test old([1, 2]) == 4

    f(x) = inv(1 - x[1]) + x[2]
    split = adaptive_map(f, [-0.5, -10.0], [0.5, 10.0]; order = 4, atol = 1.0e-6)
    @test split.converged
    @test length(split.patches) > 1
    @test all(p -> p.lower[2] == -10 && p.upper[2] == 10, split.patches)
    @test maximum(abs(split([x, 0.1]) - f([x, 0.1])) for x in range(-0.5, 0.5; length = 151)) < 1.0e-6
    @test sum(prod(p.upper - p.lower) for p in split.patches) ≈ 20
    @test all(p -> p.status == :converged && p.error_estimate[1] <= 1.0e-6, split.patches)
    @test occursin("converged", sprint(show, split))
    for p in split.patches, point in (p.lower, p.upper, p.center)
        @test abs(split(point) - f(point)) < 1.0e-6
    end

    # Coefficients above the guard order are invisible at the parent center.
    high = adaptive_map(x -> x[1]^10, [-1.0], [1.0]; order = 4, atol = 1.0e-5)
    @test length(high.patches) > 1
    @test maximum(abs(high([x]) - x^10) for x in range(-1, 1; length = 101)) < 1.0e-5
    mixed = adaptive_map(x -> (x[1] * x[2])^5, [-1, -1], [1, 1]; order = 4, atol = 1.0e-3)
    @test mixed.converged
    @test maximum(abs(mixed([x, y]) - (x * y)^5) for x in (-0.91, 0.12, 0.83), y in (-0.87, 0.3, 0.97)) < 1.0e-3
    unchecked = adaptive_map(x -> x[1]^10, [-1.0], [1.0]; order = 4, check_points = false)
    @test length(unchecked.patches) == 1 # Documented limitation of coefficient-only estimates.

    calls = Ref(0)
    counted(x) = (calls[] += 1; exp(x[1]))
    limited = adaptive_map(counted, [-1.0], [1.0]; order = 2, atol = 1.0e-12, max_patches = 3, strict = false)
    @test !limited.converged
    @test length(limited.patches) == 3
    @test all(p -> p.status == :max_patches, limited.patches)
    @test calls[] > 3 # Children reevaluate the original callback.
    @test sum(only(p.upper - p.lower) for p in limited.patches) == 2
    @test_throws ErrorException adaptive_map(exp ∘ first, [-1.0], [1.0]; order = 2, max_depth = 0)
    depth = adaptive_map(exp ∘ first, [-1.0], [1.0]; order = 2, max_depth = 1, strict = false)
    @test all(p -> p.status == :max_depth && p.depth == 1, depth.patches)

    fixed = adaptive_map(x -> [x[1] + x[2], 7], [2.0, -1.0], [2.0, 1.0])
    @test fixed([2.0, 0.5]) == [2.5, 7]
    constant = adaptive_map(x -> 3, [1.0], [1.0])
    @test constant([1.0]) == 3
    @test_throws DomainError fixed([2.01, 0.0])
    @test_throws DomainError split([NaN, 0.0])
    @test_throws DimensionMismatch split([0.0])
    adjacent = adaptive_map(x -> sin((x[1] - 1) / eps()), [1.0], [nextfloat(1.0)]; order = 1, strict = false)
    @test !adjacent.converged
    @test only(adjacent.patches).status == :roundoff
    @test adjacent([nextfloat(1.0)]) == 1.0
    rounded = adaptive_map(first, [1.0], [nextfloat(1.0)])
    @test rounded([nextfloat(1.0)]) == nextfloat(1.0)

    for T in (Float32, BigFloat)
        m = adaptive_map(x -> [exp(x[1]), x[1]^2], T[-0.5], T[0.5]; order = 4, atol = T(1.0e-4))
        @test m(T[0.1]) isa Vector{T}
        @test m(T[0.1]) ≈ [exp(T(0.1)), T(0.1)^2] atol = T(1.0e-4)
    end
    relative = adaptive_map(x -> [1.0e6 * exp(x[1]), exp(x[1])], [-0.5], [0.5]; order = 3, atol = [1.0e-2, 1.0e-8], rtol = 1.0e-5)
    @test relative.converged
    @test relative([0.17]) ≈ [1.0e6 * exp(0.17), exp(0.17)] rtol = 1.0e-5

    out, normalized, work = zeros(2), zeros(2), zeros(degree(exact) + 1)
    @test evaluate!(out, exact, [0.3, -0.4], normalized, work) === out
    @test out ≈ [1.09, -0.12]
    point = [0.3, -0.4]
    evaluate!(out, exact, point, normalized, work)
    @test @allocated(evaluate!(out, exact, point, normalized, work)) == 0
    @test_throws ArgumentError evaluate!(out, exact, out, normalized, work)
    @test_throws ArgumentError evaluate!(out, exact, [0.3, -0.4], out, work)
    @test_throws DimensionMismatch evaluate!(out, exact, [0.3, -0.4], zeros(1), work)
    @test_throws DimensionMismatch evaluate!(out, exact, [0.3, -0.4], normalized, zeros(1))

    for kwargs in (
            (; order = 0), (; guard_order = 0), (; max_depth = -1),
            (; max_patches = 0), (; atol = -1), (; atol = 0), (; rtol = Inf), (; atol = NaN),
        )
        @test_throws ArgumentError adaptive_map(first, [0.0], [1.0]; kwargs...)
    end
    @test_throws ArgumentError adaptive_map(first, Float64[], Float64[])
    @test_throws ArgumentError adaptive_map(first, [1.0], [0.0])
    @test_throws ArgumentError adaptive_map(first, [0.0], [Inf])
    @test_throws DimensionMismatch adaptive_map(first, [0.0], [1.0, 2.0])
    @test_throws DimensionMismatch adaptive_map(x -> [x[1], x[1]], [0.0], [1.0]; atol = [1.0e-3])
    @test_throws ArgumentError adaptive_map(x -> [], [0.0], [1.0])
    @test_throws ArgumentError adaptive_map(x -> Inf, [0.0], [1.0])
    @test_throws DimensionMismatch adaptive_map(x -> x[1] isa TaylorPolynomial ? [x[1]] : x[1], [0.0], [1.0])
    @test_throws ErrorException adaptive_map(x -> error("callback failure"), [0.0], [1.0])
    @test_throws ArgumentError adaptive_map(x -> (variables(1; order = 2); x), [0.0], [1.0])
    @test old([1, 2]) == 4 # Restore the caller's context even after callback failures.
    @test truncation_order() == 2
    @test coefficient_tolerance() == 1.0e-14
    variables(1; order = 2)
    @test exact([0.3, -0.4]) ≈ [1.09, -0.12]
end

@testset "ADS estimators and split directions" begin
    for estimator in (GuardedTail(), ExtrapolatedTail(), LastTerms(), LastTerms(safety = 2.08)),
            splitter in (:tail, :width)
        f(x) = inv(1 - x[1]) + x[2]
        map = adaptive_map(f, [-0.5, -2.0], [0.5, 2.0]; order = 4, atol = 1.0e-4, estimator, splitter)
        @test map.converged
        @test maximum(abs(map([x, y]) - f([x, y])) for x in range(-0.5, 0.5; length = 59), y in (-1.9, 0.3, 1.8)) < 1.0e-4
        @test sum(prod(p.upper - p.lower) for p in map.patches) ≈ 4
        if splitter == :tail
            @test all(p -> p.lower[2] == -2 && p.upper[2] == 2, map.patches)
        end
    end

    # A geometric series has an exactly known exponential fit at the next degree.
    for T in (Float32, Float64, BigFloat)
        map = adaptive_map(
            x -> inv(1 - x[1]), T[-0.5], T[0.5]; order = 4,
            estimator = ExtrapolatedTail(), check_points = false, atol = 1, max_depth = 0
        )
        @test only(only(map.patches).error_estimate) ≈ T(1) / 32
        @test map(T[0.1]) isa T
    end
    for (estimator, expected) in ((LastTerms(), 1 / 8), (LastTerms(degrees = 1, safety = 2), 1 / 8))
        map = adaptive_map(
            x -> inv(1 - x[1]), [-0.5], [0.5]; order = 4,
            estimator, check_points = false, atol = 1
        )
        @test only(only(map.patches).error_estimate) == expected
    end
    for estimator in (ExtrapolatedTail(), LastTerms())
        even = adaptive_map(x -> cos(x[1]), [-1.0], [1.0]; order = 5, atol = 1.0e-4, estimator)
        @test even.converged
        @test maximum(abs(even([x]) - cos(x)) for x in range(-1, 1; length = 71)) < 1.0e-4
        fixed = adaptive_map(x -> x[1]^2, [2.0], [2.0]; order = 1, estimator)
        @test length(fixed.patches) == 1 && fixed([2.0]) == 4
    end
    @test_throws ArgumentError LastTerms(degrees = 0)
    @test_throws ArgumentError LastTerms(safety = 0.5)
    @test_throws ArgumentError LastTerms(safety = Inf)
    @test_throws ArgumentError adaptive_map(first, Float32[0], Float32[1]; estimator = LastTerms(safety = 1.0e100))
    @test_throws ArgumentError adaptive_map(first, [0.0], [1.0]; splitter = :unknown)
    @test_throws ArgumentError adaptive_map(first, [0.0], [1.0]; estimator = LastTerms(), guard_order = 2)
end

@testset "Refine an existing partition" begin
    f(x) = exp(x[1])
    coarse = adaptive_map(f, [-1.0], [1.0]; order = 3, atol = 0.01)
    old_bounds = [(copy(p.lower), copy(p.upper)) for p in coarse.patches]
    old_value = coarse([0.35])
    fine = adaptive_map(f, coarse; atol = 1.0e-7)
    @test fine.converged && length(fine.patches) > length(coarse.patches)
    @test maximum(abs(fine([x]) - exp(x)) for x in range(-1, 1; length = 113)) < 1.0e-7
    @test all(p -> any(q -> all(q.lower .<= p.lower) && all(p.upper .<= q.upper), coarse.patches), fine.patches)
    @test [(p.lower, p.upper) for p in coarse.patches] == old_bounds
    @test coarse([0.35]) == old_value
    @test_throws ArgumentError adaptive_map(f, coarse; max_patches = length(coarse.patches) - 1)
    @test_throws ArgumentError adaptive_map(f, fine; max_depth = 0)
    # Refinement also works after reinitializing the global algebra.
    variables(3; order = 2)
    same = adaptive_map(f, fine; atol = 1.0e-7)
    @test length(same.patches) == length(fine.patches)
    @test same([0.35]) ≈ fine([0.35])
end

@testset "Checkpointed flow splitting" begin
    advance(u, span) = u / (1 - (span[2] - span[1]) * u)
    times = [0.0, 0.25, 0.5, 1.0]
    for estimator in (GuardedTail(), ExtrapolatedTail(), LastTerms())
        flow = adaptive_flow(advance, first, [0.0], [0.6], times; order = 4, atol = 1.0e-5, estimator)
        @test flow.converged && length(flow.patches) > 1
        @test maximum(abs(flow([x]) - x / (1 - x)) for x in range(0, 0.6; length = 101)) < 1.0e-5
        @test sum(only(p.upper - p.lower) for p in flow.patches) ≈ 0.6
    end

    starts = Float64[]
    initial(x) = (x[1] isa TaylorPolynomial && push!(starts, constant_term(x[1])); x[1])
    flow = adaptive_flow(advance, initial, [0.0], [0.6], times; order = 3, atol = 1.0e-6)
    @test length(starts) > 1 # Each child recomputes initial conditions.
    @test all(c -> 0 <= c <= 0.6, starts)
    backward = adaptive_flow(advance, first, [0.0], [0.5], reverse(times); order = 4, atol = 1.0e-6)
    @test backward.converged
    @test backward([0.37]) ≈ 0.37 / 1.37 atol = 1.0e-6

    # The final map can be linear while an intermediate map needs refinement.
    excursion(u, span) = span[1] == 0 ? exp.(u) : log.(u)
    initial_vector(x) = [x[1], 2 + x[1]]
    roundtrip = adaptive_flow(excursion, initial_vector, [-0.5], [0.5], [0, 1, 2]; order = 3, atol = 1.0e-5)
    @test roundtrip.converged && length(roundtrip.patches) > 1
    @test roundtrip([0.27]) ≈ [0.27, 2.27] atol = 1.0e-5
    unresolved = adaptive_flow(
        excursion, initial_vector, [-0.5], [0.5], [0, 1, 2];
        order = 3, atol = 1.0e-10, max_depth = 0, strict = false
    )
    @test !unresolved.converged
    @test only(unresolved.patches).status == :max_depth
    @test unresolved([0.0]) ≈ [0.0, 2.0] atol = 1.0e-14 # Final, not stopped-at-first-checkpoint state.
    @test maximum(only(unresolved.patches).error_estimate) > 1.0e-10
    @test_throws ErrorException adaptive_flow(advance, first, [0.0], [0.6], times; order = 2, max_depth = 0)
    budget = adaptive_flow(
        advance, first, [0.0], [0.6], times; order = 2, atol = 1.0e-12,
        max_patches = 3, strict = false
    )
    @test !budget.converged && length(budget.patches) == 3
    @test sum(only(p.upper - p.lower) for p in budget.patches) ≈ 0.6

    for T in (Float32, BigFloat)
        flow = adaptive_flow(advance, first, T[0], T[0.4], T[0, 0.5, 1]; order = 4, atol = T(1.0e-4))
        @test flow(T[0.3]) isa T
        @test flow(T[0.3]) ≈ T(0.3) / (1 - T(0.3)) atol = T(1.0e-4)
    end
    old, = variables(1; order = 3)
    for bad_times in ([0.0], [0.0, NaN], [0.0, 0.0], [0.0, 1.0, 0.5])
        @test_throws ArgumentError adaptive_flow(advance, first, [0.0], [0.5], bad_times)
    end
    @test_throws DimensionMismatch adaptive_flow((u, _) -> [u], first, [0.0], [0.5], [0, 1])
    @test_throws ErrorException adaptive_flow((u, _) -> error("failed integration"), first, [0.0], [0.5], [0, 1])
    @test old(0.5) == 0.5
end

@testset "Online flow splitting" begin
    function advance(u, span, monitor)
        previous = first(span)
        for time in range(span...; length = 5)[2:end]
            u = u / (1 - (time - previous) * u)
            monitor(u, time) && return (; state = u, time)
            previous = time
        end
        return (; state = u, time = last(span))
    end
    for estimator in (GuardedTail(), ExtrapolatedTail(), LastTerms()), checked in (false, true)
        map = adaptive_flow(
            advance, first, [0.0], [0.6], (0, 1);
            order = 4, atol = 1.0e-5, estimator, check_points = checked
        )
        @test map.converged && length(map.patches) > 1
        @test maximum(abs(map([x]) - x / (1 - x)) for x in range(0, 0.6; length = 87)) < 1.0e-5
    end
    back = adaptive_flow(advance, first, [0.0], [0.5], (1, 0); order = 4, atol = 1.0e-6)
    @test back.converged
    @test back([0.31]) ≈ 0.31 / 1.31 atol = 1.0e-6
    limited = adaptive_flow(
        advance, first, [0.0], [0.6], (0, 1);
        order = 2, atol = 1.0e-12, max_patches = 2, strict = false
    )
    @test !limited.converged && length(limited.patches) == 2
    for p in limited.patches
        c = only(p.center)
        @test limited([c]) ≈ c / (1 - c) atol = 1.0e-13 # Each leaf reaches final time.
    end
    # Check the endpoint even if the propagator omits its last callback.
    endpoint(u, span, monitor) = (; state = u / (1 - (span[2] - span[1]) * u), time = last(span))
    final = adaptive_flow(endpoint, first, [0.0], [0.6], (0, 1); order = 4, atol = 1.0e-5)
    @test final.converged && length(final.patches) > 1
    @test final([0.43]) ≈ 0.43 / 0.57 atol = 1.0e-5
    @test_throws ArgumentError adaptive_flow(advance, first, [0.0], [0.6], (1, 1))
    @test_throws ArgumentError adaptive_flow(advance, first, [0.0], [0.6], (0, Inf))
    @test_throws ArgumentError adaptive_flow((u, span, check) -> u, first, [0.0], [0.5], (0, 1))
    @test_throws ArgumentError adaptive_flow((u, span, check) -> (; state = u, time = 0.5), first, [0.0], [0.5], (0, 1))
    @test_throws ArgumentError adaptive_flow((u, span, check) -> (check(u, 2); (; state = u, time = 2)), first, [0.0], [0.5], (0, 1))
end
