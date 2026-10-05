using Test, DifferentialAlgebra, ForwardDiff, IntervalArithmetic, LinearAlgebra

exact_derivative(f, x, n) = n == 0 ? f(x) : ForwardDiff.derivative(t -> exact_derivative(f, t, n - 1), x)
const RContinuity = Rational{BigInt}

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
    script = "using DifferentialAlgebra; @assert !DifferentialAlgebra.isinitialized(); @assert !any(m->nameof(m) in (:IntervalArithmetic,:ForwardDiff),values(Base.loaded_modules)); using ForwardDiff; @assert !DifferentialAlgebra.isinitialized(); a=adaptive_map(v->v[1]^2,[-1.],[1.];order=2); s=continuous_map(v->v[1]^2,a); @assert ForwardDiff.derivative(x->s([x]),.25)≈.5; @assert !DifferentialAlgebra.isinitialized()"
    @test success(`$(Base.julia_cmd()) --startup-file=no --project=$(dirname(Base.active_project())) -e $script`)
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
