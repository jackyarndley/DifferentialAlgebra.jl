using Test, DifferentialAlgebra
using IntervalArithmetic

include("support/intervals.jl")
using .IntervalTestSupport: DA, IA, IF, IB, subset, interval_contains, sameinterval

function contract_contains(parent::NamedTuple, child::NamedTuple)
    return all(parent.lower .<= child.lower) && all(child.upper .<= parent.upper)
end
contract_contains(parent::ConvexPolygon, child::ConvexPolygon) =
    all(v -> DA.polygon_contains(parent, v), polygon_vertices(child))
function contract_contains(
        parent::AbstractVector{<:Interval}, child::AbstractVector{<:Interval}
    )
    return all(subset.(child, parent))
end
contract_value(actual::Interval, exact) = interval_contains(actual, exact)
contract_value(actual::Real, exact) = isapprox(actual, exact; atol = 1.0e-12)
contract_same(a::Interval, b::Interval) = sameinterval(a, b)
contract_same(a::Real, b::Real) = a == b
contract_same(a::AbstractVector, b::AbstractVector) = all(contract_same.(a, b))

@testset "Snapshot collections use owned vector storage" begin
    for T in (Float64, BigFloat), budget in (0, 32 * 1024^2)
        f(v) = [v[1]^2 + v[2]^2, v[1]]
        for kind in (:ordinary, :validated, :polygon)
            estimator = kind === :ordinary ? LastTerms() : IntervalBound()
            build(atol) = adaptive_map(
                f, T[-1, -1], T[1, 1]; estimator, order = 1, atol,
                splitter = kind === :polygon ? :oriented : :width,
                directions = kind === :polygon ? :axes : nothing,
                table_bytes = budget
            )
            fits = [build(atol) for atol in (1 // 2, 1 // 8)]
            @test typeof(fits[1]) === typeof(fits[2])
            @test length.(getproperty.(fits, :patches)) == [4, 16]
            for fit in fits
                @test isconcretetype(fieldtype(typeof(fit), :patches))
                @test fit.patches isa Vector
                @test @inferred(degree(fit)) == 1
                before = fit([1 // 8, -1 // 4])
                patches = fit.patches
                count = length(patches)
                empty!(patches)
                @test length(fit.patches) == count
                @test contract_same(fit([1 // 8, -1 // 4]), before)
                patch = first(fit.patches)
                raw = patch isa PolygonPatch ? patch._patch : patch
                if raw isa TaylorModelPatch
                    @test isconcretetype(fieldtype(typeof(raw), :models))
                    models = raw.models
                    @test models isa Vector && raw.error_bounds isa Vector
                    empty!(models)
                    @test noutputs(raw) == 2
                    bounds = raw.error_bounds
                    bounds[1] = interval(T, 999)
                    @test IA.sup(raw.error_bounds[1]) < 999
                end
                smooth = continuous_map(f, fit; table_bytes = budget)
                @test isconcretetype(fieldtype(typeof(smooth), :_patches))
                @test smooth._patches isa Vector
                saved = smooth([1 // 8, -1 // 4])
                @test copy(smooth)([1 // 8, -1 // 4]) == saved
            end
        end
        box = [interval(T, -1, 1), interval(T, 2)]
        x, y = taylor_models(box; order = 2, table_bytes = budget)
        for name in (:box, :center, :radius, :normalized)
            @test isconcretetype(fieldtype(typeof(x._coordinates), name))
            @test getproperty(x._coordinates, name) isa Vector{Interval{T}}
        end
        saved = domain(x)
        box[1] = interval(T, 999)
        @test all(sameinterval.(domain(x), saved))
        @test copy(x)._coordinates.box !== x._coordinates.box
    end
end

@testset "Refinement inherits owned construction options" begin
    square = ConvexPolygon([(-1, -1), (1, -1), (1, 1), (-1, 1)])
    for estimator in (GuardedTail(), LastTerms(), IntervalBound()),
            polygon in (false, true), budget in (0, 32 * 1024^2)
        guard = estimator isa GuardedTail ? 3 : 0
        names, atol = ["q", "r"], [1 // 8, 1 // 16]
        observed = Tuple{Int, Int, Tuple{String, String}}[]
        function f(v)
            p = v[1] isa TaylorPolynomial ? v[1] :
                v[1] isa TaylorModel ? getfield(v[1], :_polynomial) : nothing
            if p !== nothing
                push!(
                    observed, (
                        truncation_order(), p.algebra.basis.table_bytes,
                        Tuple(p.algebra.names),
                    )
                )
            end
            return [v[1]^2 + v[2], v[2]]
        end
        options = (;
            estimator, order = 2, atol, rtol = 0, guard_order = guard,
            max_depth = 0, max_patches = 17, check_points = false,
            strict = false, splitter = :width, names, table_bytes = budget,
        )
        source = polygon ? adaptive_map(f, square; directions = :axes, options...) :
            adaptive_map(f, [-1, -1], [1, 1]; options...)
        atol[1], names[1] = 99, "changed"
        owned = copy(source)
        @test owned._options.atol !== source._options.atol
        @test owned._options.names !== source._options.names
        empty!(observed)
        refined = adaptive_map(f, owned)
        @test !isempty(observed)
        @test all(==((2 + guard, budget, ("q", "r"))), observed)
        @test refined._options.atol == [1 // 8, 1 // 16]
        @test refined._options.max_depth == 0 && refined._options.max_patches == 17
        @test !refined._options.strict && !refined._options.check_points
        @test refined._options.splitter == :width && iszero(refined._options.rtol)
        @test length(refined.patches) == length(source.patches) == 1
        actual = refined([1 // 4, -1 // 2])
        @test all(contract_value.(actual, [-7 // 16, -1 // 2]))
        revised = adaptive_map(
            f, refined; atol = [1 // 4, 1 // 8], max_depth = 1, max_patches = 19
        )
        @test revised._options.atol == [1 // 4, 1 // 8]
        @test revised._options.max_depth == 1 && revised._options.max_patches == 19
        @test refined._options.atol == [1 // 8, 1 // 16]
    end
end

@testset "Refinement beyond the constructor's default patch limit" begin
    f(v) = v[1]^2
    source = adaptive_map(
        f, [-1.0], [1.0]; order = 1, estimator = LastTerms(),
        atol = 1 // 2^22, max_patches = 3000, max_depth = 12, table_bytes = 0
    )
    # A linear Taylor fit to x² has exact endpoint error r². This dyadic
    # tolerance requires radius 1/2048 and hence 2048 equal-width patches.
    @test length(source.patches) == 2048
    refined = adaptive_map(f, source)
    @test refined.converged && length(refined.patches) == 2048
    @test abs(refined([0.123]) - 0.123^2) <= 1 // 2^22
end

@testset "Changing estimators adopts compatible guard defaults" begin
    f(v) = v[1]^2
    source = adaptive_map(
        f, [-1.0], [1.0]; order = 2, guard_order = 3,
        atol = 1 // 8, max_depth = 0, strict = false
    )
    retained = adaptive_map(f, source; estimator = LastTerms())
    @test retained._options.guard_order == 0
    guarded = adaptive_map(f, retained; estimator = GuardedTail())
    @test guarded._options.guard_order == 2
    explicit = adaptive_map(f, retained; estimator = GuardedTail(), guard_order = 1)
    @test explicit._options.guard_order == 1
end

@testset "Invalid tolerances reject before callback evaluation" begin
    square = ConvexPolygon([(-1, -1), (1, -1), (1, 1), (-1, 1)])
    for geometry in (:box, :oriented, :polygon)
        build(f; kwargs...) = geometry === :polygon ?
            adaptive_map(f, square; kwargs...) :
            adaptive_map(
                f, [-1, -1], [1, 1];
                splitter = geometry === :oriented ? :oriented : :tail, kwargs...
            )
        for estimator in (GuardedTail(), ExtrapolatedTail(), LastTerms(), IntervalBound())
            for atol in (-1, 0, Inf, NaN, [-1], [0], [Inf], [NaN], Float64[])
                calls = Ref(0)
                f(v) = (calls[] += 1; v[1])
                @test_throws ArgumentError build(f; estimator, atol)
                @test iszero(calls[])
            end
        end
        for atol in (
                IA.interval(Float64, -1, 1), IA.interval(Float64, 0),
                convert(IF, 0.125),
            )
            calls = Ref(0)
            f(v) = (calls[] += 1; v[1])
            @test_throws ArgumentError build(f; estimator = IntervalBound(), atol)
            @test iszero(calls[])
        end
        for atol in (1 // 8, IA.interval(Float64, 1 // 8, 1 // 4), IA.exact(0.125))
            fit = build(v -> v[1]; estimator = IntervalBound(), atol)
            @test fit.converged
        end
        for estimator in (GuardedTail(), IntervalBound())
            fit = build(v -> v[1]; estimator, atol = 1 // 8)
            for atol in (-1, NaN, [0])
                calls = Ref(0)
                f(v) = (calls[] += 1; v[1])
                @test_throws ArgumentError continuous_map(f, fit; atol)
                @test iszero(calls[])
            end
        end
    end
    fit = adaptive_map(first, [-1.0], [1.0]; atol = 0, rtol = 0.001)
    @test fit.converged && fit([0.25]) == 0.25
end

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
            narrow(v) = exp((v[1] - 1) / IA.exact(eps()))
            roundoff = polygon ?
                adaptive_map(narrow, DA.polygon_box(lo, hi); estimator, directions = :axes, order = 2, atol = 1.0e-12, check_points = false, strict = false) :
                adaptive_map(narrow, lo, hi; estimator, order = 2, atol = 1.0e-12, check_points = false, strict = false)
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
            failure = estimator isa IntervalBound ? ErrorException : DomainError
            @test_throws failure build(
                v -> throw(DomainError(v, "callback failure"));
                max_depth = 0
            )
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
    explicit_fit = adaptive_map(v -> 2, [interval(-1, 1)]; estimator = IntervalBound(), order = 3, atol = 1.0e-6, splitter = :width)
    canonical = adaptive_map(v -> 2, [interval(-1, 1)]; estimator = IntervalBound())
    @test max_order(explicit_fit) == 3 && max_order(canonical) == 5
    @test degree(explicit_fit) == degree(canonical) == 0
    @test sameinterval(enclose(explicit_fit), enclose(canonical))
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
        owned = copy(fit)
        @test owned.nodes !== fit.nodes
        @test first(owned.patches).models[1]._coordinates === first(owned.patches).models[2]._coordinates
        @test first(owned.patches).models[1]._coordinates !== first(fit.patches).models[1]._coordinates
        saved = enclose(fit)
        initialize!(1, 2)
        @test all(sameinterval.(saved, enclose(fit))) && all(sameinterval.(saved, enclose(owned)))
    end
end

@testset "Patch bounds do not alias the map or neighboring patches" begin
    for T in (Float64, BigFloat)
        lower, upper = T[-1, -1], T[1, 1]
        fit = adaptive_map(
            v -> exp(v[1]), lower, upper;
            order = 3, atol = 1.0e-3, check_points = false
        )
        @test length(fit.patches) > 1
        @test all(p -> p.lower !== fit.lower && p.upper !== fit.upper, fit.patches)
        @test all(
            i == j || fit.patches[i].lower !== fit.patches[j].lower
                for i in eachindex(fit.patches), j in eachindex(fit.patches)
        )
        bounds = [deepcopy((p.lower, p.upper)) for p in fit.patches]
        saved = domain(fit)
        fit.patches[1].lower[1] = 42
        @test domain(fit) == saved
        @test all(
            i -> (fit.patches[i].lower, fit.patches[i].upper) == bounds[i],
            2:length(fit.patches)
        )
        @test lower == T[-1, -1] && upper == T[1, 1]
        if T === BigFloat
            @test all(p -> all(i -> p.upper[i] !== fit.upper[i], 1:2), fit.patches)
        end
    end
end
