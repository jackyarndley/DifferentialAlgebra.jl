using Test, DifferentialAlgebra, IntervalArithmetic
isdefined(@__MODULE__, :IntervalTestSupport) || include("support/intervals.jl")
using .IntervalTestSupport: IA, DA, subset, interval_contains, sameinterval

@testset "Certified box ADS and portable model snapshots" begin
    for T in (Float64, BigFloat), budget in (0, 32 * 1024^2)
        z = IA.interval(T, -1, 1)
        caller, = variables(1; order = 4)
        previous = caller^3
        ctx = DA.CURRENT_ALGEBRA[]
        calls = Ref(0)
        f(v) = (calls[] += 1; v[1]^2)
        a = validated_adaptive_map(f, [z]; order = 1, atol = 1 // 16, table_bytes = budget)
        @test DA.CURRENT_ALGEBRA[] === ctx && degree(previous) == 3
        @test a.converged && length(a.patches) == 4 && calls[] == 7
        @test noutputs(a) == 1 && nvariables(a) == 1 && max_order(a) == 1
        @test all(p -> p.status == :converged && p.depth == 2, a.patches)
        # Independent uniform oracle: on c+r*xi, x² equals c²+2cr*xi+r²*xi².
        # These dyadic boxes have exact retained coefficients and r=1/4.
        for patch in a.patches
            @test sameinterval(only(patch.error_bounds), IA.interval(T, 0, 1 // 16))
            @test sameinterval(remainder(only(patch.models)), IA.interval(T, 0, 1 // 16))
            @test IA.isguaranteed(enclose(only(patch.models)))
        end
        @test interval_contains(enclose(a), 0) && interval_contains(enclose(a), 1)
        @test IA.decoration(enclose(a)) == IA.com && IA.isguaranteed(enclose(a))
        @test interval_contains(evaluate(a, [0]), 0)
        @test interval_contains(a([1]), 1) && interval_contains(a([-1]), 1)
        @test subset(IA.interval(T, 0, 1 // 4), enclose(a, [IA.interval(T, -1 // 2, 1 // 2)]))
        @test_throws DimensionMismatch evaluate(a, [])
        @test_throws DimensionMismatch enclose(a, [z, z])
        @test_throws DomainError evaluate(a, [2])
        @test_throws DomainError enclose(a, [IA.interval(T, 0, 2)])
        @test_throws ArgumentError evaluate(a, [z])
        @test_throws ArgumentError enclose(a, [IA.bareinterval(-1, 1)])
        @test_throws ArgumentError evaluate(a, [convert(Interval{T}, T(0))])
        saved = enclose(a)
        owned = domain(a); owned[1] = IA.interval(T, 100)
        @test sameinterval(enclose(a), saved)
        b = copy(a)
        @test b.patches[1].models[1]._coefficients !== a.patches[1].models[1]._coefficients
        @test sameinterval(enclose(deepcopy(a)), saved)
        initialize!(2, 2)
        @test sameinterval(enclose(a), saved) && sameinterval(enclose(b), saved)
        # Mixed discarded terms are bounded over the full square, not samples.
        mixed = validated_adaptive_map(
            v -> [v[1] * v[2], 3, v[1] - v[1]], [z, z];
            order = 1, atol = [1 // 4, 1 // 100, 1 // 100], table_bytes = budget
        )
        @test mixed.converged && length(mixed.patches) == 4
        @test noutputs(mixed) == 3
        @test all(p -> sameinterval(p.error_bounds[1], IA.interval(T, -1 // 4, 1 // 4)), mixed.patches)
        @test all(p -> IA.isthinzero(p.error_bounds[2]) && IA.isthinzero(p.error_bounds[3]), mixed.patches)
        values = evaluate(mixed, [1 // 2, -1 // 2])
        @test interval_contains(values[1], -1 // 4) && interval_contains(values[2], 3) && IA.isthinzero(values[3])
        constant = validated_adaptive_map(v -> 2 // 3, [z]; atol = 1 // 100)
        @test constant.converged && length(constant.patches) == 1
        @test interval_contains(enclose(constant), 2 // 3)
        fixed = validated_adaptive_map(v -> v[1] * v[2], [z, IA.interval(T, 2)]; order = 1)
        @test fixed.converged && length(fixed.patches) == 1
        @test interval_contains(evaluate(fixed, [1, 2]), 2)
        @test_throws DomainError evaluate(fixed, [0, 3])
        # Coefficient uncertainty must enter acceptance even with zero remainder.
        uncertain = validated_adaptive_map(
            v -> IA.interval(T, -1 // 8, 1 // 8), [z];
            atol = 1 // 16, max_depth = 0, strict = false
        )
        @test !uncertain.converged && only(uncertain.patches).status == :max_depth
        @test IA.isthinzero(remainder(only(only(uncertain.patches).models)))
        @test sameinterval(only(only(uncertain.patches).error_bounds), IA.interval(T, -1 // 8, 1 // 8))
        incoming = validated_adaptive_map(
            v -> TaylorModel(polynomial(v[1]), IA.interval(T, -1 // 8, 1 // 8), v[1]), [z];
            atol = 1 // 16, max_patches = 2, strict = false
        )
        @test !incoming.converged && length(incoming.patches) == 2
        @test all(p -> p.status == :max_patches && sameinterval(only(p.error_bounds), IA.interval(T, -1 // 8, 1 // 8)), incoming.patches)
        @test_throws ErrorException validated_adaptive_map(v -> v[1]^2, [z]; order = 1, atol = 1 // 16, max_depth = 0)
        @test_throws ErrorException validated_adaptive_map(v -> v[1]^2, [z]; order = 1, atol = 1 // 16, max_patches = 1)
        adjacent = IA.interval(T, one(T), nextfloat(one(T)))
        unsplittable = validated_adaptive_map(v -> IA.interval(T, -1, 1), [adjacent]; atol = 1 // 2, strict = false)
        @test only(unsplittable.patches).status == :roundoff
        for fn in (inv, log, sqrt)
            @test_throws DomainError validated_adaptive_map(v -> fn(1 + v[1]), [IA.interval(T, -2, 2)])
        end
        caller, = variables(1; order = 3)
        ctx = DA.CURRENT_ALGEBRA[]
        @test_throws ArgumentError validated_adaptive_map(v -> (set_truncation_order!(1); v[1]), [z]; order = 3)
        @test_throws ArgumentError validated_adaptive_map(v -> (set_coefficient_tolerance!(1); v[1]), [z])
        @test DA.CURRENT_ALGEBRA[] === ctx && degree(caller) == 1
        @test_throws ArgumentError validated_adaptive_map(v -> polynomial(v[1]), [z])
        @test_throws ArgumentError validated_adaptive_map(v -> convert(Interval{T}, T(1)), [z])
        @test_throws ArgumentError validated_adaptive_map(v -> v[1], [z]; atol = 0)
        @test_throws DimensionMismatch validated_adaptive_map(v -> [v[1], v[1]], [z]; atol = [1])
        @test_throws ArgumentError validated_adaptive_map(v -> v[1], [z]; max_patches = 0)
        @test_throws ArgumentError validated_adaptive_map(v -> v[1], [z]; order = 0)
        @test_throws ArgumentError validated_adaptive_map(v -> v[1], [IA.emptyinterval(T)])
        x, = taylor_models([z]; order = 2)
        source = x^2
        snapshot = compile(source)
        before = enclose(snapshot)
        @test sameinterval(before, enclose(source))
        @test snapshot._coefficients !== source._polynomial.coeffs
        @test sameinterval(enclose(copy(snapshot)), before)
        @test_throws DomainError snapshot([2])
        @test_throws DimensionMismatch snapshot([])
        @test_throws ArgumentError snapshot([z])
        initialize!(1, 1)
        @test sameinterval(enclose(snapshot), before)
        @test interval_contains(snapshot([1 // 2]), 1 // 4)
    end
    # Analytical Taylor-theorem oracle: the degree-one expansion of exp on a
    # child radius r has remainder at most exp(upper)*r²/2, plus enclosed retained
    # coefficients. Acceptance uses these validated quantities uniformly.
    a = validated_adaptive_map(v -> exp(v[1]), [IA.interval(Float64, -1 // 4, 1 // 4)]; order = 1, atol = 1 // 100)
    @test a.converged && length(a.patches) > 1
    @test all(p -> IA.sup(abs(only(p.error_bounds))) <= 1 // 100, a.patches)
    for patch in a.patches
        box = only(domain(patch))
        theoretical = exp(box) * IA.pown(box - IA.interval(IA.mid(box)), 2) / IA.interval(2)
        r = remainder(only(patch.models))
        @test IA.inf(r) >= 0 # Convexity of exp, with exact dyadic normalization.
        @test IA.sup(r) <= IA.sup(theoretical) # Tightness permits a narrower valid remainder.
    end
    a = validated_adaptive_map(v -> v[1]^2, [IA.interval(-1, 1)]; order = 1, atol = 1 // 16)
    IA.configure(; rounding = :none)
    try
        @test_throws ArgumentError enclose(a)
        @test_throws ArgumentError evaluate(a, [0])
        @test_throws ArgumentError validated_adaptive_map(v -> v[1], [IA.interval(-1, 1)])
    finally
        IA.configure(; rounding = :correct)
    end
end
