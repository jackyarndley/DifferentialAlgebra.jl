using Test, DifferentialAlgebra, LinearAlgebra

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

    for budget in (0, 32 * 1024^2), estimator in (GuardedTail(), ExtrapolatedTail(), LastTerms())
        caller, = variables(1; order = 4)
        ctx = DA.CURRENT_ALGEBRA[]
        f(v) = (v[1] + v[2])^2
        a = adaptive_map(f, [-1, -1], [1, 1]; order = 1, atol = 1 // 16, estimator, splitter = :oriented, table_bytes = budget)
        @test DA.CURRENT_ALGEBRA[] === ctx && degree(caller) == 1
        @test a isa PiecewisePolygonMap && a.converged
        @test split_directions(a) == [1 1; -1 1]
        @test max_order(a) == 1 && degree(a) <= 1 && nvariables(a) == 2 && noutputs(a) == 1
        @test sum(p -> domain_area(domain(p)), a.patches) == domain_area(domain(a))
        # Coverage and disjoint interiors are checked exactly, independently of
        # function samples. Every pair has at most a shared edge or vertex.
        @test all(p -> all(v -> DA.polygon_contains(square, v), polygon_vertices(domain(p))), a.patches)
        @test all(i -> all(j -> DA.polygon_intersection(domain(a.patches[i]), domain(a.patches[j])) === nothing, (i + 1):length(a.patches)), eachindex(a.patches))
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
        refined = adaptive_map(f, a; order = 1, atol = 1 // 64, table_bytes = budget)
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
    script = "using DifferentialAlgebra; @assert !DifferentialAlgebra.isinitialized(); a=adaptive_map(v->v[1]+v[2],[-1,-1],[1,1];splitter=:oriented,order=1); @assert a([0.5,0.25])==0.75; @assert !DifferentialAlgebra.isinitialized(); @assert !any(m->nameof(m)==:IntervalArithmetic,values(Base.loaded_modules)); try adaptive_map(identity,[-1,-1],[1,1];estimator=IntervalBound()); error(\"missing failure\"); catch e; @assert e isa ArgumentError; end"
    @test success(`$(Base.julia_cmd()) --startup-file=no --project=$(dirname(Base.active_project())) -e $script`)
end
