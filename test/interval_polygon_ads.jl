using Test, DifferentialAlgebra, IntervalArithmetic, LinearAlgebra

@testset "IntervalBound selector and certified polygon ADS" begin
    for T in (Float64, BigFloat), budget in (0, 32 * 1024^2)
        z = IA.interval(T, -1, 1)
        square = ConvexPolygon([(-one(T), -one(T)), (one(T), -one(T)), (one(T), one(T)), (-one(T), one(T))])
        f(v) = (v[1] + v[2])^2
        caller, = variables(1; order = 4)
        ctx = DA.CURRENT_ALGEBRA[]
        a = adaptive_map(f, [z, z]; estimator = IntervalBound(), splitter = :oriented, order = 1, atol = 1 // 16, table_bytes = budget)
        @test DA.CURRENT_ALGEBRA[] === ctx && degree(caller) == 1
        @test a.converged && a isa PiecewisePolygonMap && length(a.patches) == 8
        @test split_directions(a) == [1 1; -1 1]
        @test sum(p -> domain_area(domain(p)), a.patches) == 4
        @test all(p -> IA.isguaranteed(only(p.error_bounds)) && IA.decoration(only(p.error_bounds)) == IA.com, a.patches)
        # Independent exact identity f(x)=z1² and dyadic normalization. Each
        # patch has r1=1/4, so the full discarded polynomial is (1/16)*xi1².
        for patch in a.patches
            @test sameinterval(only(patch.error_bounds), IA.interval(T, 0, 1 // 16))
            @test sameinterval(remainder(only(patch._patch.models)), IA.interval(T, 0, 1 // 16))
        end
        @test subset(IA.interval(T, 0, 4), enclose(a))
        # Exact rational point oracles supplement the uniform remainder proof.
        for point in ([0, 0], [1, 1], [-1, -1], [1, -1], [1 // 2, 1 // 4])
            value = a(point)
            @test interval_contains(value, sum(point)^2)
            @test IA.isguaranteed(value) && IA.decoration(value) == IA.com
        end
        @test interval_contains(a([IA.interval(T, 0), IA.interval(T, 0)]), 0)
        @test subset(IA.interval(T, 0, 1 // 4), enclose(a, [IA.interval(T, -1 // 2, 1 // 2), IA.interval(T, 0)]))
        @test interval_contains(enclose(a, [IA.interval(T, 1), IA.interval(T, 1)]), 4)
        triangle = ConvexPolygon([(zero(T), zero(T)), (one(T), zero(T)), (zero(T), one(T))])
        @test subset(IA.interval(T, 0, 1), enclose(a, triangle))
        @test_throws DimensionMismatch a([0])
        @test_throws DimensionMismatch enclose(a, [z])
        @test_throws DomainError a([0, 2])
        @test_throws DomainError enclose(a, [z, IA.interval(T, 0, 2)])
        @test_throws ArgumentError a([z, z])
        @test_throws ArgumentError a([convert(Interval{T}, zero(T)), 0])
        @test_throws ArgumentError enclose(a, [IA.bareinterval(-1, 1), IA.bareinterval(-1, 1)])
        @test_throws ArgumentError enclose(a, [IA.emptyinterval(T), z])
        before = enclose(a)
        dirs = split_directions(a); dirs[1, 1] = 100
        vertices = polygon_vertices(domain(a)); vertices[1] = (100 // 1, 100 // 1)
        errors = first(a.patches).error_bounds
        @test sameinterval(only(errors), only(first(a.patches)._patch.error_bounds))
        if T === BigFloat
            @test IA.inf(only(errors)) !== IA.inf(only(first(a.patches)._patch.error_bounds))
        end
        b = copy(a)
        @test first(b.patches)._patch.models[1]._coefficients !== first(a.patches)._patch.models[1]._coefficients
        initialize!(2, 3)
        @test sameinterval(enclose(a), before) && sameinterval(enclose(b), before) && sameinterval(enclose(deepcopy(a)), before)
        refined = adaptive_map(f, a; order = 1, atol = 1 // 64, table_bytes = budget)
        @test refined.converged && length(refined.patches) == 16
        @test sum(p -> domain_area(domain(p)), refined.patches) == 4
        @test all(p -> IA.sup(abs(only(p.error_bounds))) <= 1 // 64, refined.patches)
        same = validated_adaptive_map(f, [z, z]; order = 1, atol = 1 // 16, splitter = :oriented, table_bytes = budget)
        @test sameinterval(enclose(same), before)
        endpoints = adaptive_map(f, T[-1, -1], T[1, 1]; estimator = IntervalBound(), splitter = :oriented, order = 1, atol = 1 // 16, table_bytes = budget)
        @test sameinterval(enclose(endpoints), before)
        axis = adaptive_map(f, [z, z]; estimator = IntervalBound(), order = 1, atol = 1 // 16, table_bytes = budget)
        @test axis isa PiecewiseTaylorModel && axis.converged && length(axis.patches) > length(a.patches)
        # A general nonorthogonal frame is inverted exactly, not transposed.
        affine = adaptive_map(v -> [v[1] / 3 - v[2] / 5, 2 // 3], triangle; estimator = IntervalBound(), order = 1, directions = [1 2; -3 1], table_bytes = budget)
        @test length(affine.patches) == 1 && affine.converged
        @test interval_contains(affine([1 // 4, 1 // 4])[1], 1 // 30)
        @test interval_contains(affine([0, 0])[2], 2 // 3)
        @test_throws DomainError affine([3 // 4, 3 // 4])
        @test_throws DomainError enclose(affine, [IA.interval(T, 0, 1), IA.interval(T, 0, 1)])
        # The function is positive on the physical triangle, but its axis
        # cover includes x+y=2. A valid center or polygon is not a cover proof.
        @test_throws DomainError adaptive_map(v -> log(3 // 2 - v[1] - v[2]), triangle; estimator = IntervalBound(), directions = :axes)
        zero_fit = adaptive_map(v -> v[1] - v[1], square; estimator = IntervalBound(), order = 1)
        @test length(zero_fit.patches) == 1 && IA.isthinzero(enclose(zero_fit))
        constant = adaptive_map(v -> 2 // 3, square; estimator = IntervalBound(), order = 1)
        @test length(constant.patches) == 1 && interval_contains(enclose(constant), 2 // 3)
        # Includes a genuinely mixed omitted monomial; oracle |xi1*xi2|<=1.
        mixed = adaptive_map(v -> v[1] * v[2], square; estimator = IntervalBound(), order = 1, atol = 1, directions = [1 0; 1 1], table_bytes = budget)
        @test mixed.converged
        @test interval_contains(mixed([1, -1]), -1) && interval_contains(mixed([1, 1]), 1)
        incoming = adaptive_map(
            v -> TaylorModel(polynomial(v[1]), IA.interval(T, -1 // 8, 1 // 8), v[1]),
            square; estimator = IntervalBound(), order = 1, atol = 1 // 16,
            directions = :axes, max_depth = 0, strict = false, table_bytes = budget
        )
        @test !incoming.converged && sameinterval(only(only(incoming.patches).error_bounds), IA.interval(T, -1 // 8, 1 // 8))
        for fn in (inv, log, sqrt)
            @test_throws DomainError adaptive_map(v -> fn(1 + v[1]), square; estimator = IntervalBound(), directions = :axes)
        end
        @test_throws ArgumentError adaptive_map(v -> polynomial(v[1]), square; estimator = IntervalBound(), directions = :axes)
        @test_throws ArgumentError adaptive_map(v -> convert(Interval{T}, one(T)), square; estimator = IntervalBound(), directions = :axes)
        @test_throws ArgumentError adaptive_map(v -> IA.interval(T, 1, 2, IA.trv), square; estimator = IntervalBound(), directions = :axes)
        @test_throws ArgumentError adaptive_map(identity, square; estimator = IntervalBound(), rtol = 1 // 10)
        @test_throws ArgumentError adaptive_map(identity, square; estimator = IntervalBound(), check_points = true)
        @test_throws ArgumentError adaptive_map(v -> (set_coefficient_tolerance!(1); v[1]), square; estimator = IntervalBound(), directions = :axes)
        @test_throws ArgumentError adaptive_map(v -> (set_truncation_order!(1); v[1]), square; estimator = IntervalBound(), directions = :axes)
        @test_throws DimensionMismatch adaptive_map(v -> [v[1], v[2]], square; estimator = IntervalBound(), atol = [1], directions = :axes)
        @test_throws ArgumentError adaptive_map(v -> IA.interval(T, 1) * v[1], square; estimator = GuardedTail(), directions = :axes)
        @test_throws ArgumentError adaptive_map(v -> v[1], [z, z]; estimator = LastTerms())
        @test_throws ArgumentError adaptive_map(v -> v[1], [z, z]; estimator = IntervalBound(), guard_order = 1)
        @test_throws ArgumentError adaptive_map(v -> v[1], [z, z]; estimator = IntervalBound(), splitter = :unknown)
        # Six-variable box ADS still uses independent coordinates and the
        # native engine; the new geometry does not redirect this path.
        six = adaptive_map(v -> sum(v), fill(z, 6); estimator = IntervalBound(), order = 1, table_bytes = budget)
        @test nvariables(six) == 6 && length(six.patches) == 1
        @test interval_contains(six(fill(1 // 3, 6)), 2)
    end
    setprecision(BigFloat, 160) do
        z = IA.interval(BigFloat, -1, 1)
        a = adaptive_map(v -> (v[1] + v[2])^2, [z, z]; estimator = IntervalBound(), splitter = :oriented, order = 1, atol = 1 // 16)
        setprecision(BigFloat, 64) do
            @test interval_contains(a([1 // 3, 1 // 3]), 4 // 9)
            @test subset(IA.interval(BigFloat, 0, 4), enclose(a))
        end
    end
    @test_throws ArgumentError adaptive_map(identity, [IA.interval(-1, 1)]; estimator = IntervalBound(), splitter = :oriented)
    IA.configure(; rounding = :none)
    try
        @test_throws ArgumentError adaptive_map(identity, [IA.interval(-1, 1), IA.interval(-1, 1)]; estimator = IntervalBound(), splitter = :oriented)
    finally
        IA.configure(; rounding = :correct)
    end
end
