# Included in DifferentialAlgebraIntervalArithmeticExt. Exact polygon geometry
# and its inverse frame live in DA; every coefficient/range operation uses IA.
DA.polygon_mid(x::Interval) = IA.mid(x)

function DA.ads_interval_domain(a::DA.PiecewisePolygonMap{T, P, DA.IntervalBound}) where {T, P}
    patch = first(DA.ads_patches(a))._patch
    return first(DA.model_snapshots(patch))._coordinates.box
end

function DA.polygon_candidate(
        f, polygon::DA.ConvexPolygon{T}, B, A, options, ctx,
        estimator::DA.IntervalBound; probe::Bool = false
    ) where {T}
    validated_rounding()
    T in (Float64, BigFloat) || throw(ArgumentError("IntervalBound supports Float64 or BigFloat polygon endpoints"))
    I = Interval{T}
    lo, hi = DA.polygon_projection_bounds(polygon, B)
    box = [checked_interval(IA.interval(T, lo[i], hi[i]); guaranteed = true) for i in 1:2]
    retain_decoration!(box, options.interval_domain)
    z = coordinate_models(box, options.order)
    x = [sum(IA.interval(T, A[i, j]) * z[j] for j in 1:2) for i in 1:2]
    candidate = validated_candidate(f, x, first(z), box, ctx, options)
    geometry = DA.ads_geometry(IA.inf.(box), IA.sup.(box))
    if candidate isa DA.ADSDomainFailure
        cover = (; geometry..., projection = (lo, hi), directions = nothing)
        return DA.ADSDomainFailure(cover, candidate.cause, candidate.zero_score)
    end
    values = candidate.payload.models
    directions = if probe
        polynomials = [m._polynomial for m in values]
        radii = [IA.inf(r) for r in first(z)._coordinates.radius]
        DA.polygon_sensitivity(polynomials, radii, options.scales, candidate.tolerance, T)
    else
        nothing
    end
    cover = (; geometry..., projection = (lo, hi), directions)
    return DA.ADSCandidate(cover, candidate.payload, candidate.errors, candidate.contributions, candidate.tolerance, candidate.scalar, candidate.accepted)
end

function polygon_interval_point(a, point)
    point isa Union{Tuple, AbstractVector} && length(point) == 2 || throw(DimensionMismatch("Provide two physical coordinates"))
    point isa AbstractVector && Base.require_one_based_indexing(point)
    if any(x -> x isa AbstractIrrational, point)
        domain = DA.ads_interval_domain(a)
        return [point_interval(eltype(domain), point[i], domain[i]) for i in 1:2]
    end
    x = map(point) do value
        if value isa Interval
            checked_interval(value; guaranteed = true)
            IA.isthin(value) || throw(ArgumentError("Use enclose for non-thin query intervals"))
            return DA.polygon_real(IA.inf(value))
        end
        return DA.polygon_real(value)
    end
    exact = (x[1], x[2])
    routing = DA.polygon_query(a, exact)
    routing.exact = exact
    return routing
end

function polygon_interval_constraints(a, query)
    if query isa DA.ConvexPolygon
        all(v -> DA.polygon_contains(a._domain, v), query._vertices) || throw(DomainError(query, "Query polygon lies outside the validity domain"))
        return [DA.polygon_halfspace(query, i) for i in eachindex(query._vertices)]
    end
    interval_box(query, 2; guaranteed = true)
    lo = DA.polygon_real.(IA.inf.(query))
    hi = DA.polygon_real.(IA.sup.(query))
    corners = ((lo[1], lo[2]), (lo[1], hi[2]), (hi[1], lo[2]), (hi[1], hi[2]))
    all(v -> DA.polygon_contains(a._domain, v), corners) || throw(DomainError(query, "Query box lies outside the polygon domain"))
    R = DA.PolygonReal
    return [((one(R), zero(R)), hi[1]), ((-one(R), zero(R)), -lo[1]), ((zero(R), one(R)), hi[2]), ((zero(R), -one(R)), -lo[2])]
end

function polygon_snapshot_bounds(
        models::AbstractVector{<:DA.CompiledTaylorModel{I}}, points, B
    ) where {I <: Interval}
    return polygon_projected_snapshot_bounds(
        models, map(x -> DA.polygon_project(B, x), points)
    )
end
function polygon_projected_snapshot_bounds(
        models::AbstractVector{<:DA.CompiledTaylorModel{I}}, projected
    ) where {I <: Interval}
    model = first(models)
    T = IA.numtype(I)
    box = map(1:2) do i
        value = IA.interval(T, minimum(x -> x[i], projected), maximum(x -> x[i], projected))
        # Both sets contain every exact projected query point. Intersection also
        # handles a caller reducing BigFloat precision since construction.
        IA.intersect_interval(value, model._coordinates.box[i]; dec = :auto)
    end
    return snapshot_bounds(models, box)
end

function DA.polygon_interval_query(a::DA.PiecewisePolygonMap, query = nothing; point = false)
    validated_rounding()
    a._estimator isa DA.IntervalBound || throw(ArgumentError("Select IntervalBound for polygon enclosures"))
    routing = point ? polygon_interval_point(a, query) : nothing
    # A transcendental coordinate has no exact rational representation. Query
    # its outward enclosure, admitting it only when the entire box is inside.
    routing isa AbstractVector && return DA.polygon_interval_query(a, routing)
    exact = point ? DA.polygon_exact!(routing) : nothing
    constraints = query === nothing || point ? nothing : polygon_interval_constraints(a, query)
    result = nothing
    query_points = if point
        routing
    elseif query === nothing
        a._domain._vertices
    elseif query isa DA.ConvexPolygon
        query._vertices
    else
        lo, hi = DA.polygon_real.(IA.inf.(query)), DA.polygon_real.(IA.sup.(query))
        ((lo[1], lo[2]), (lo[1], hi[2]), (hi[1], lo[2]), (hi[1], hi[2]))
    end
    projected = point ? (DA.polygon_project(a._directions, exact),) : nothing
    DA.ads_visit_intersections(a.nodes, query_points) do index
        patch = DA.ads_patches(a)[index]
        points = if point
            DA.polygon_contains(patch._domain, exact) ? (exact,) : ()
        elseif query === nothing
            patch._domain._vertices
        else
            remaining = patch._domain._vertices
            for (normal, cut) in constraints
                remaining = DA.polygon_clip_points(remaining, normal, cut)
                isempty(remaining) && break
            end
            remaining
        end
        isempty(points) && return nothing
        models = DA.model_snapshots(patch._patch)
        values = point ? polygon_projected_snapshot_bounds(models, projected) :
            polygon_snapshot_bounds(models, points, a._directions)
        result = hull_outputs!(result, values)
        return nothing
    end
    result === nothing && throw(ErrorException("Polygon partition has a coverage gap"))
    retain_decoration!(result, query)
    return a._scalar ? only(result) : result
end
