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
