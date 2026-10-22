using Test, DifferentialAlgebra
using LinearAlgebra

include("support/processes.jl")

@testset "Exact polygon geometry and ordinary oriented ADS" begin
    DA = DifferentialAlgebra
    vertices = [[-1, -1], [1, -1], [1, 1], [-1, 1]]
    square = ConvexPolygon(vertices)
    vertices[1][1] = 100
    @test domain_area(square) == 4
    @test polygon_vertices(square) == polygon_vertices(ConvexPolygon(reverse(polygon_vertices(square))))
    owned = polygon_vertices(square); owned[1] = (100 // 1, 100 // 1)
    @test first(polygon_vertices(square)) == (-1, -1)
    @test domain_area(copy(square)) == domain_area(deepcopy(square)) == 4
    @test polygon_vertices(ConvexPolygon([(0.1, 0), (1, 0), (1, 1)]))[1][1] == Rational{BigInt}(0.1)
    for v in (
            [(0, 0), (1, 0)], [(0, 0), (1, 0), (2, 0)],
            [(0, 0), (1, 0), (0, 0)], [(0, 0), (1, 1), (0, 1), (1, 0)],
            [(0, 0), (2, 0), (1, 1 // 2), (2, 2), (0, 2)],
            [(0, 0), (Inf, 0), (0, 1)],
        )
        @test_throws ArgumentError ConvexPolygon(v)
    end
    @test_throws DimensionMismatch ConvexPolygon([(0, 0, 0), (1, 0), (0, 1)])
    B, A = DA.polygon_frame([1 1; -1 1])
    @test B * A == Matrix{Rational{BigInt}}(I, 2, 2)
    left = DA.polygon_clip(square, (1, 1), 0)
    right = DA.polygon_clip(square, (-1, -1), 0)
    @test domain_area(left) == domain_area(right) == 2
    @test DA.polygon_intersection(left, right) === nothing
    @test DA.polygon_contains(left, (0, 0)) && DA.polygon_contains(right, (0, 0))
    @test length(DA.polygon_clip_points(polygon_vertices(left), (-1, -1), 0)) == 2
    @test_throws ArgumentError DA.polygon_frame([1 1; 2 2])
    @test_throws DimensionMismatch DA.polygon_frame([1, 1])

    for estimator in (GuardedTail(), ExtrapolatedTail(), LastTerms())
        caller, = variables(1; order = 4)
        ctx = DA.CURRENT_ALGEBRA[]
        f(v) = (v[1] + v[2])^2
        a = adaptive_map(
            f, [-1, -1], [1, 1]; order = 1, atol = 1 // 16, estimator, splitter = :oriented
        )
        patches = a.patches
        @test DA.CURRENT_ALGEBRA[] === ctx && degree(caller) == 1
        @test a isa PiecewisePolygonMap && a.converged
        @test split_directions(a) == [1 1; -1 1]
        @test max_order(a) == 1 && degree(a) <= 1 && nvariables(a) == 2 && noutputs(a) == 1
        @test sum(p -> domain_area(domain(p)), a.patches) == domain_area(domain(a))
        # Coverage and disjoint interiors are checked exactly, independently of
        # function samples. Every pair has at most a shared edge or vertex.
        @test all(p -> all(v -> DA.polygon_contains(square, v), polygon_vertices(domain(p))), a.patches)
        @test all(
            i -> all(
                j -> DA.polygon_intersection(
                    domain(patches[i]), domain(patches[j])
                ) === nothing,
                (i + 1):length(patches)
            ), eachindex(patches)
        )
        @test all(x -> any(p -> DA.polygon_contains(domain(p), x), a.patches), [(x // 4, y // 4) for x in -4:4 for y in -4:4])
        # Independent analytical error on every parallelogram: the omitted
        # term of z1² is r1²*xi1². This is a uniform bound, not a sampling proof.
        for patch in a.patches
            @test patch._patch.radius[1]^2 <= 1 // 16
            @test patch.error_estimate[1] <= 1 // 16
        end
        @test abs(a([1 // 2, 1 // 4]) - 9 // 16) <= 1 // 16
        @test abs(a([1, 1]) - 4) <= 1 // 16
        @test_throws DimensionMismatch a([0])
        @test_throws DomainError a([2, 0])
        @test_throws ArgumentError enclose(a)
        @test_throws ArgumentError a.patches[1].error_bounds
        before = a([1 // 2, 1 // 4])
        owned = split_directions(a); owned[1, 1] = 100
        b = copy(a)
        @test b.patches[1]._patch.map !== a.patches[1]._patch.map
        initialize!(2, 3)
        @test a([1 // 2, 1 // 4]) == b([1 // 2, 1 // 4]) == before
        refined = adaptive_map(f, a; order = 1, atol = 1 // 64)
        @test length(refined.patches) >= length(a.patches)
        @test sum(p -> domain_area(domain(p)), refined.patches) == 4
        @test all(p -> any(q -> all(v -> DA.polygon_contains(domain(q), v), polygon_vertices(domain(p))), a.patches), refined.patches)
    end
    triangle = ConvexPolygon([(0, 0), (1, 0), (0, 1)])
    affine = adaptive_map(v -> [2 + v[1] / 3 - v[2] / 5, 3], triangle; order = 1, directions = [1 2; -3 1])
    @test length(affine.patches) == 1
    @test affine([1 // 4, 1 // 4]) ≈ [2 + 1 / 12 - 1 / 20, 3]
    @test_throws DomainError affine([3 // 4, 3 // 4])
    unresolved = adaptive_map(v -> v[1]^2, square; order = 1, atol = 1 // 100, max_depth = 0, strict = false, directions = :axes)
    @test !unresolved.converged && only(unresolved.patches).status == :max_depth
    @test_throws ErrorException adaptive_map(v -> v[1]^2, square; order = 1, atol = 1 // 100, max_patches = 1, directions = :axes)
    @test_throws ArgumentError adaptive_map(identity, [0, 0, 0], [1, 1, 1]; splitter = :oriented)
    @test_throws ArgumentError adaptive_map(identity, [0, 0], [0, 1]; splitter = :oriented)
    @test_throws ArgumentError adaptive_map(identity, [0, 0], [1, 1]; directions = :axes)
    @test_throws ArgumentError adaptive_map(identity, square; directions = [1 1; 2 2])
    @test_throws ArgumentError adaptive_map(v -> (set_coefficient_tolerance!(1); v[1]), square; directions = :axes)
    @test_throws ArgumentError adaptive_map(v -> (set_truncation_order!(1); v[1]), square; directions = :axes)
    # Ordinary polygon geometry works without loading the optional extension.
    script = raw"""
    using DifferentialAlgebra
    @assert !DifferentialAlgebra.isinitialized()
    a = adaptive_map(v -> v[1] + v[2], [-1, -1], [1, 1];
        splitter = :oriented, order = 1)
    @assert a([0.5, 0.25]) == 0.75
    @assert !DifferentialAlgebra.isinitialized()
    @assert !any(m -> nameof(m) == :IntervalArithmetic, values(Base.loaded_modules))
    try
        adaptive_map(identity, [-1, -1], [1, 1]; estimator = IntervalBound())
        error("missing extension failure")
    catch err
        @assert err isa ArgumentError
    end
    """
    @test fresh_julia(script)
end

@testset "Thin polygon covers round outward" begin
    tri = ConvexPolygon([(1.0, 1.0), (2.0, 1.0), (2.0, 1.0 + 1.0e-10)])
    fit = adaptive_map(
        v -> inv(v[1] - 0.999999), tri;
        order = 3, atol = 1.0e-12, directions = :axes, strict = false, max_depth = 20
    )
    @test !fit.converged
    @test all(p -> p.status in (:converged, :roundoff, :max_depth, :max_patches), fit.patches)
    @test sum(p -> domain_area(domain(p)), fit.patches) == domain_area(tri)
    B = split_directions(fit)
    @test all(fit.patches) do patch
        box = patch._patch
        all(polygon_vertices(domain(patch))) do vertex
            all(1:2) do i
                projected = sum(B[i, j] * vertex[j] for j in 1:2)
                Rational{BigInt}(box.lower[i]) <= projected <= Rational{BigInt}(box.upper[i])
            end
        end
    end
    @test all(v -> isfinite(fit(v)), polygon_vertices(tri))
end


include("support/polygon_queries.jl")
using .PolygonQueryTestSupport: exact_patch, exact_evaluate, face_points

@testset "Filtered polygon routing retains exact decisions across frames" begin
    DA = DifferentialAlgebra
    f(v) = v[1]^2 + v[2]^2
    for T in (Float32, Float64, BigFloat)
        polygon = ConvexPolygon(
            [
                (-one(T), -one(T)), (one(T), -one(T)),
                (one(T), one(T)), (-one(T), one(T)),
            ]
        )
        original = adaptive_map(
            f, polygon; directions = :axes,
            order = 1, atol = 1 // 8
        )
        fit = adaptive_map(
            f, original; directions = [1 1 // 3; -1 // 3 1],
            atol = 1 // 16
        )
        refined = adaptive_map(
            f, fit; directions = [1 -1 // 5; 1 // 5 1],
            atol = 1 // 32
        )
        for a in (original, fit, refined)
            query = DA.polygon_query(a, [0.123, 0.234])
            @test DA.ads_point_patch(a.nodes, query) == exact_patch(a, [0.123, 0.234])
            @test query.exact === nothing
            inside, outside = [], []
            for face in face_points(a)
                points = (collect(face), Float64[face...])
                append!(inside, points)
                point = Float64[face...]
                for axis in 1:2, neighbor in (prevfloat, nextfloat)
                    x = copy(point)
                    x[axis] = neighbor(x[axis])
                    if DA.polygon_contains(a._domain, Rational{BigInt}.(x))
                        push!(inside, x)
                    else
                        push!(outside, x)
                    end
                end
            end
            @test all(inside) do point
                DA.ads_point_patch(a.nodes, DA.polygon_query(a, point)) ==
                    exact_patch(a, point) && isapprox(a(point), exact_evaluate(a, point))
            end
            @test all(outside) do point
                try
                    a(point)
                    false
                catch err
                    err isa DomainError
                end
            end
            b = copy(a)
            @test b._lookup.directions !== a._lookup.directions
            @test b._lookup.bounds[1] !== a._lookup.bounds[1]
            @test b._lookup.domain[1] !== a._lookup.domain[1]
            @test a(T[1 // 7, 2 // 9]) ≈ exact_evaluate(a, T[1 // 7, 2 // 9])
        end
        @test size(refined._lookup.bounds[1], 1) >= 4
        # The root oblique cut crosses the origin. A nonzero tangent point
        # needs the exact fallback to decide the side of its rounded projection.
        oblique = adaptive_map(
            f, polygon; directions = [1 1 // 3; -1 // 3 1],
            order = 1, atol = 1 // 8
        )
        n = first(oblique.nodes).normal
        point = [n[2] / 7, -n[1] / 7]
        query = DA.polygon_query(oblique, point)
        @test DA.ads_point_patch(oblique.nodes, query) == exact_patch(oblique, point)
        @test query.exact !== nothing
        for invalid in ([Inf, 0], [NaN, 0])
            @test_throws ArgumentError refined(invalid)
        end
    end
end

@testset "Automatic frames use nonlinear sensitivity" begin
    for T in (Float64, BigFloat)
        polygon = ConvexPolygon(
            [
                (-one(T), -one(T)), (one(T), -one(T)),
                (one(T), one(T)), (-one(T), one(T)),
            ]
        )
        for phi in (exp, sin, x -> x^3)
            nonlinear(v) = phi(v[1] + 2v[2])
            shifted(v) = nonlinear(v) + 1000v[1] + 250v[2] + 7
            # Every nonconstant coefficient of the nonlinear gradient is
            # parallel to (1,2), independently of the affine part or estimator.
            expected = T[1 // 2 1; -1 1 // 2]
            for estimator in (GuardedTail(), LastTerms(), ExtrapolatedTail())
                a = adaptive_map(
                    nonlinear, polygon; order = 3, atol = 1 // 4,
                    estimator, max_patches = 512
                )
                b = adaptive_map(
                    shifted, polygon; order = 3, atol = 1 // 4,
                    estimator, max_patches = 512
                )
                @test T.(split_directions(a)) ≈ expected atol = 128eps(T)
                @test T.(split_directions(b)) ≈ expected atol = 128eps(T)
                @test a.converged && b.converged
            end
        end
        affine = adaptive_map(
            v -> 2v[1] + v[2], polygon; order = 1,
        )
        @test split_directions(affine) == Matrix{Int}(I, 2, 2)
        # Retaining a larger working basis must not hide the requested order.
        orders = Int[]
        callback(v) = (push!(orders, truncation_order()); v[1] + v[2])
        adaptive_map(
            callback, polygon; order = 1, estimator = LastTerms(),
            check_points = false
        )
        @test orders == [3, 1]
        empty!(orders)
        cubic(v) = (push!(orders, truncation_order()); (v[1] + 2v[2])^3)
        fit = adaptive_map(
            cubic, polygon; order = 2, atol = 1 // 2,
            check_points = false
        )
        @test length(orders) == 2length(fit.patches)
        @test T.(split_directions(fit)) ≈ T[1 // 2 1; -1 1 // 2] atol = 128eps(T)
    end
end
