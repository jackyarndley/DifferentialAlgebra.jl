using Test, DifferentialAlgebra, IntervalArithmetic
isdefined(@__MODULE__, :IntervalTestSupport) || include("support/intervals.jl")
using .IntervalTestSupport: DA, IA, IF, IB, subset, interval_contains, sameinterval

function contract_contains(parent::NamedTuple, child::NamedTuple)
    return all(parent.lower .<= child.lower) && all(child.upper .<= parent.upper)
end
contract_contains(parent::ConvexPolygon, child::ConvexPolygon) = all(v -> DA.polygon_contains(parent, v), polygon_vertices(child))
contract_contains(parent::AbstractVector{<:Interval}, child::AbstractVector{<:Interval}) = all(subset.(child, parent))
contract_value(actual::Interval, exact) = interval_contains(actual, exact)
contract_value(actual::Real, exact) = isapprox(actual, exact; atol = 1.0e-12)
contract_same(a::Interval, b::Interval) = sameinterval(a, b)
contract_same(a::Real, b::Real) = a == b
contract_same(a::AbstractVector, b::AbstractVector) = all(contract_same.(a, b))

@testset "Shared ADS construction contracts" begin
    square = ConvexPolygon([(-1, -1), (1, -1), (1, 1), (-1, 1)])
    for estimator in (GuardedTail(), ExtrapolatedTail(), LastTerms(), IntervalBound()), polygon in (false, true)
        @testset "$(typeof(estimator)) / $(polygon ? :polygon : :box)" begin
            build(f; kwargs...) = polygon ?
                adaptive_map(f, square; estimator, directions = :axes, check_points = false, kwargs...) :
                adaptive_map(f, [-1.0, -1.0], [1.0, 1.0]; estimator, check_points = false, kwargs...)
            caller, = variables(1; order = 4)
            ctx = DA.CURRENT_ALGEBRA[]
            affine(v) = [2 + v[1] / 2 - v[2] / 4, 3]
            exact = build(affine; order = 4, atol = 1.0e-5)
            @test exact.converged && length(exact.patches) == 1
            @test degree(exact) == 1 && max_order(exact) == 4
            @test nvariables(exact) == 2 && noutputs(exact) == 2
            @test exact._estimator === estimator
            @test DA.CURRENT_ALGEBRA[] === ctx && degree(caller) == 1
            @test all(p -> nvariables(p) == 2 && noutputs(p) == 2 && max_order(p) == 4 && p.status == :converged, exact.patches)
            for point in ([0, 0], [-1, -1], [1, 1], [1 // 4, -1 // 2])
                @test all(contract_value.(exact(point), affine(point)))
            end
            saved_domain = domain(exact)
            unchanged = exact([1 // 4, -1 // 2])
            low_degree = adaptive_map(affine, exact; atol = 1.0e-5, check_points = false)
            @test max_order(low_degree) == 4 && degree(low_degree) == 1
            @test low_degree._estimator === estimator
            @test contract_contains(domain(exact), domain(low_degree))
            @test contract_contains(domain(low_degree), saved_domain)
            @test contract_same(exact([1 // 4, -1 // 2]), unchanged)
            if !(estimator isa IntervalBound)
                @test_throws ArgumentError enclose(exact)
            else
                @test all(p -> all(IA.isguaranteed, p.error_bounds), exact.patches)
            end

            f(v) = exp(v[2])
            for (limits, status) in ((; max_depth = 0) => :max_depth, (; max_patches = 1) => :max_patches)
                unresolved = build(f; order = 2, atol = 1.0e-12, strict = false, limits...)
                @test !unresolved.converged && only(unresolved.patches).status == status
                @test_throws ErrorException build(f; order = 2, atol = 1.0e-12, limits...)
            end
            # Every coordinate's rounded midpoint equals an endpoint. Exact
            # polygon area is positive, but no representable fitted cut remains.
            lo, hi = [1.0, 1.0], fill(nextfloat(1.0), 2)
            narrow(v) = exp((v[1] - 1) / eps())
            roundoff = polygon ?
                adaptive_map(narrow, DA.polygon_box(lo, hi); estimator, directions = :axes, order = 1, atol = 1.0e-12, check_points = false, strict = false) :
                adaptive_map(narrow, lo, hi; estimator, order = 1, atol = 1.0e-12, check_points = false, strict = false)
            @test !roundoff.converged && only(roundoff.patches).status == :roundoff
            if !polygon
                fixed = adaptive_map(affine, [2.0, -1.0], [2.0, 1.0]; estimator, order = 4, atol = 1.0e-5, check_points = false)
                @test all(contract_value.(fixed([2, 0]), affine([2, 0])))
                @test_throws DomainError fixed([3, 0])
            end
            calls = Ref(0)
            original(v) = (calls[] += 1; exp(v[2]))
            coarse = build(original; order = 3, atol = 0.01)
            previous_calls = calls[]
            previous_value = coarse([1 // 4, 1 // 4])
            previous_domains = domain.(coarse.patches)
            fine = adaptive_map(original, coarse; atol = 1.0e-4, check_points = false)
            @test fine.converged && calls[] - previous_calls >= length(coarse.patches)
            @test max_order(fine) == 3 && fine._estimator === estimator
            @test all(p -> any(q -> contract_contains(q, domain(p)), previous_domains), fine.patches)
            @test all(contract_contains.(previous_domains, domain.(coarse.patches)))
            @test contract_same(coarse([1 // 4, 1 // 4]), previous_value)
            @test_throws ArgumentError adaptive_map(original, coarse; max_patches = length(coarse.patches) - 1)
            @test_throws ArgumentError adaptive_map(original, fine; max_depth = 0)
            before = fine([1 // 4, 1 // 4])
            initialize!(1, 3)
            @test contract_same(fine([1 // 4, 1 // 4]), before)
            @test_throws DomainError fine([2, 0])
            @test_throws DimensionMismatch fine([0])

            # No sampled calls or frame probe: force a shape change on the first child.
            calls[] = 0
            changing(v) = (calls[] += 1; calls[] == 1 ? exp(v[2]) : [exp(v[2])])
            @test_throws DimensionMismatch build(changing; order = 2, atol = 1.0e-12)
            @test_throws DomainError build(v -> throw(DomainError(v, "callback failure")))
            @test nvariables() == 3 && max_order() == 1
            @test_throws ArgumentError build(v -> (initialize!(2, 2); v[1]))
            @test nvariables() == 3 && max_order() == 1
        end
    end
end

@testset "Fixed subnormal box coordinates" begin
    for T in (Float32, Float64)
        endpoint = nextfloat(zero(T))
        fit = adaptive_map(first, T[endpoint], T[endpoint]; order = 3)
        patch = only(fit.patches)
        @test only(patch.center) == endpoint && iszero(only(patch.radius))
        @test degree(fit) == 0 && max_order(fit) == 3
        @test fit(T[endpoint]) == endpoint
        @test_throws DomainError fit(T[0])
        refined = adaptive_map(first, fit)
        @test only(only(refined.patches).center) == endpoint
        @test iszero(only(only(refined.patches).radius))
        @test refined(T[endpoint]) == endpoint
    end
end

@testset "Ordinary BigFloat domain ownership" begin
    for estimator in (GuardedTail(), ExtrapolatedTail(), LastTerms())
        lower, upper = BigFloat[-1, 3], BigFloat[1, 3]
        f(v) = exp(v[1]) + v[2]
        fit = adaptive_map(f, lower, upper; estimator, order = 3, atol = 0.01, check_points = false)
        refined = adaptive_map(f, fit; atol = 1.0e-4, check_points = false)
        @test length(refined.patches) > length(fit.patches)
        for owned in (fit, refined, copy(fit))
            for (bounds, input) in ((owned.lower, lower), (owned.upper, upper))
                @test bounds !== input && bounds == input
                @test all(i -> bounds[i] !== input[i], eachindex(input))
            end
            returned = domain(owned)
            @test returned.lower == lower && returned.upper == upper
            @test all(i -> returned.lower[i] !== owned.lower[i], eachindex(lower))
            @test all(i -> returned.upper[i] !== owned.upper[i], eachindex(upper))
        end
        @test all(i -> refined.lower[i] !== fit.lower[i], eachindex(lower))
        @test all(i -> refined.upper[i] !== fit.upper[i], eachindex(upper))
        source_endpoints = [x for p in fit.patches for x in (p.lower..., p.upper...)]
        @test all(p -> all(x -> all(y -> x !== y, source_endpoints), (p.lower..., p.upper...)), refined.patches)
        owned_copy = copy(fit)
        @test all(i -> owned_copy.lower[i] !== fit.lower[i], eachindex(lower))
        @test all(i -> owned_copy.upper[i] !== fit.upper[i], eachindex(upper))
        saved, saved_refined = fit(BigFloat[0, 3]), refined(BigFloat[0, 3])
        lower[1], upper[1] = -100, 100
        @test domain(fit) == (lower = BigFloat[-1, 3], upper = BigFloat[1, 3])
        @test fit(BigFloat[0, 3]) == saved && refined(BigFloat[0, 3]) == saved_refined
        @test_throws DomainError fit(BigFloat[2, 3])
    end
end

@testset "Interval coefficients discard; models enclose discarded terms" begin
    for I in (IF, IB), budget in (0, 32 * 1024^2)
        T = IA.numtype(I)
        box = [interval(T, -1, 1)]
        p, = variables(I, 1; order = 1, table_bytes = budget)
        @test iszero(p * p) && degree(p * p) == 0
        @test IA.isthinzero(enclose(p * p, box))
        x, = taylor_models(box; order = 1, table_bytes = budget)
        m = x * x
        @test iszero(polynomial(m)) && max_order(m) == 1
        @test sameinterval(remainder(m), interval(T, 0, 1))
        @test sameinterval(enclose(m), interval(T, 0, 1))
        @test sameinterval(enclose(compile(m)), enclose(m))
    end
end

@testset "Validated error-driven split directions" begin
    for T in (Float64, BigFloat), budget in (0, 32 * 1024^2)
        for f in (v -> [2^40 * v[1], exp(v[2])], v -> 2^20 * v[1] + exp(v[2]))
            fit = adaptive_map(
                f, T[-1, -1], T[1, 1]; estimator = IntervalBound(), splitter = :tail,
                order = 3, atol = 1.0e-4, max_depth = 1, strict = false, table_bytes = budget
            )
            @test first(fit.nodes).axis == 2 && first(fit.nodes).midpoint == 0
            @test all(p -> sameinterval(domain(p)[1], interval(T, -1, 1)), fit.patches)
            @test !fit.converged && all(p -> p.status == :max_depth, fit.patches)
            # Tightness check, separate from the certified inclusion proof:
            # halving the exponential coordinate reduces the parent's ~0.113 bound.
            @test maximum(p -> maximum(IA.sup.(abs.(p.error_bounds))), fit.patches) < 0.01
            @test all(interval_contains.(fit([0, 0]) isa AbstractVector ? fit([0, 0]) : [fit([0, 0])], f([0, 0]) isa AbstractVector ? f([0, 0]) : [f([0, 0])]))
        end
        fixed = adaptive_map(v -> exp(v[2]), [interval(T, 2), interval(T, -1, 1)]; estimator = IntervalBound(), order = 3, atol = 1.0e-4, max_depth = 1, strict = false)
        @test first(fixed.nodes).axis == 2
        adjacent = interval(T, one(T), nextfloat(one(T)))
        unresolved = adaptive_map(v -> interval(T, -1, 1), [interval(T, 2), adjacent]; estimator = IntervalBound(), atol = 0.1, strict = false)
        @test only(unresolved.patches).status == :roundoff
        @test_throws DomainError fixed([3, 0])
    end
    # No nonlinear tail or directional uncertainty: relative widths decide.
    fallback = adaptive_map(v -> TaylorModel(polynomial(v[1]), interval(-1, 1), v[1]), [-1.0, -2.0], [1.0, 2.0]; estimator = IntervalBound(), order = 3, atol = 0.1, max_depth = 1, strict = false)
    @test first(fallback.nodes).axis == 1
end

@testset "Closed tree queries and shared multi-output snapshots" begin
    legacy = validated_adaptive_map(v -> 2, [interval(-1, 1)])
    canonical = adaptive_map(v -> 2, [interval(-1, 1)]; estimator = IntervalBound())
    @test max_order(legacy) == 3 && max_order(canonical) == 5
    @test degree(legacy) == degree(canonical) == 0
    @test sameinterval(enclose(legacy), enclose(canonical))
    for T in (Float64, BigFloat)
        f(v) = [v[1]^2, v[1]^3, 7]
        box = [interval(T, -1, 1)]
        fit = adaptive_map(f, box; estimator = IntervalBound(), order = 1, atol = 1 // 16, splitter = :width)
        @test length(fit.nodes) == 2length(fit.patches) - 1
        visited = Int[]
        DA.ads_visit_intersections(i -> push!(visited, i), fit.nodes, ([T(1 // 7)], [T(1 // 7)]))
        @test length(visited) == 1
        empty!(visited)
        DA.ads_visit_intersections(i -> push!(visited, i), fit.nodes, ([zero(T)], [zero(T)]))
        @test length(visited) == 2
        for query in ([interval(T, 0)], [interval(T, -1 // 2, 1 // 2)], box)
            # Independent scan oracle for branch coverage, including closed faces.
            scan = nothing
            for patch in fit.patches
                intersection = [IA.intersect_interval(only(query), only(domain(patch)); dec = :auto)]
                IA.isempty_interval(only(intersection)) && continue
                bounds = [enclose(m, intersection) for m in patch.models]
                scan = scan === nothing ? bounds : [IA.hull(x, y; dec = :auto) for (x, y) in zip(scan, bounds)]
            end
            @test all(sameinterval.(enclose(fit, query), scan))
        end
        for patch in fit.patches
            a, b, constant = patch.models
            @test a._coordinates === b._coordinates === constant._coordinates
            @test a._exponents === b._exponents === constant._exponents
            @test degree(constant) == 0 && max_order(constant) == 1
            @test a._coefficients !== b._coefficients
            owned = copy(a)
            @test owned._coordinates !== a._coordinates && owned._exponents !== a._exponents
            if T === BigFloat
                @test IA.inf(owned._coordinates.box[1]) !== IA.inf(a._coordinates.box[1])
                @test IA.inf(owned._coefficients[1]) !== IA.inf(a._coefficients[1])
            end
            @test sameinterval(enclose(owned), enclose(a))
        end
        compatible = validated_adaptive_map(f, box; order = 1, atol = 1 // 16, splitter = :width)
        @test all(all(sameinterval.(domain(p), domain(q))) for (p, q) in zip(fit.patches, compatible.patches))
        @test all(sameinterval.(enclose(fit), enclose(compatible)))
        @test typeof(fit._estimator) === typeof(compatible._estimator)
        owned = copy(fit)
        @test owned.nodes !== fit.nodes
        @test first(owned.patches).models[1]._coordinates === first(owned.patches).models[2]._coordinates
        @test first(owned.patches).models[1]._coordinates !== first(fit.patches).models[1]._coordinates
        saved = enclose(fit)
        initialize!(1, 2)
        @test all(sameinterval.(saved, enclose(fit))) && all(sameinterval.(saved, enclose(owned)))
    end
end
