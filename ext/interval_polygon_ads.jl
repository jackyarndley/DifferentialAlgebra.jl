# Included in DifferentialAlgebraIntervalArithmeticExt. Exact polygon geometry
# and its inverse frame live in DA; every coefficient/range operation uses IA.
DA.polygon_mid(x::Interval) = IA.mid(x)

function DA.polygon_candidate(
        f, polygon::DA.ConvexPolygon{T}, B, A, options, ctx,
        estimator::DA.IntervalBound
    ) where {T}
    validated_rounding()
    T in (Float64, BigFloat) || throw(ArgumentError("IntervalBound supports Float64 or BigFloat polygon endpoints"))
    I = Interval{T}
    lo, hi = DA.polygon_projection_bounds(polygon, B)
    box = [checked_interval(IA.interval(T, lo[i], hi[i]); guaranteed = true) for i in 1:2]
    z = coordinate_models(box, options.order)
    x = [sum(IA.interval(T, A[i, j]) * z[j] for j in 1:2) for i in 1:2]
    candidate = validated_candidate(f, x, first(z), box, ctx, options)
    values = candidate.payload.models
    geometry = DA.ads_geometry(IA.inf.(box), IA.sup.(box))
    polynomials = [m._polynomial for m in values]
    radii = [IA.inf(r) for r in first(z)._coordinates.radius]
    directions = DA.polygon_sensitivity(polynomials, radii, options.scales, candidate.tolerance, T)
    cover = (; geometry..., projection = (lo, hi), directions)
    return DA.ADSCandidate(cover, candidate.payload, candidate.errors, candidate.contributions, candidate.tolerance, candidate.scalar, candidate.accepted)
end

function polygon_query_constraints(p::DA.ConvexPolygon)
    v = p._vertices
    return map(eachindex(v)) do i
        a, b = v[i], v[mod1(i + 1, length(v))]
        normal = (b[2] - a[2], a[1] - b[1])
        (normal, normal[1] * a[1] + normal[2] * a[2])
    end
end

function polygon_interval_point(a, point)
    point isa Union{Tuple, AbstractVector} && length(point) == 2 || throw(DimensionMismatch("Provide two physical coordinates"))
    point isa AbstractVector && Base.require_one_based_indexing(point)
    x = map(point) do value
        if value isa Interval
            checked_interval(value; guaranteed = true)
            IA.isthin(value) || throw(ArgumentError("Use enclose for non-thin query intervals"))
            return DA.polygon_real(IA.inf(value))
        end
        return DA.polygon_real(value)
    end
    exact = (x[1], x[2])
    DA.polygon_contains(a._domain, exact) || throw(DomainError(point, "Point lies outside the polygon domain"))
    return exact
end

function polygon_interval_constraints(a, query)
    if query isa DA.ConvexPolygon
        all(v -> DA.polygon_contains(a._domain, v), query._vertices) || throw(DomainError(query, "Query polygon lies outside the validity domain"))
        return polygon_query_constraints(query)
    end
    interval_box(query, 2; guaranteed = true)
    lo = DA.polygon_real.(IA.inf.(query))
    hi = DA.polygon_real.(IA.sup.(query))
    corners = ((lo[1], lo[2]), (lo[1], hi[2]), (hi[1], lo[2]), (hi[1], hi[2]))
    all(v -> DA.polygon_contains(a._domain, v), corners) || throw(DomainError(query, "Query box lies outside the polygon domain"))
    R = DA.PolygonReal
    return [((one(R), zero(R)), hi[1]), ((-one(R), zero(R)), -lo[1]), ((zero(R), one(R)), hi[2]), ((zero(R), -one(R)), -lo[2])]
end

function polygon_snapshot_bounds(models::Tuple{Vararg{DA.CompiledTaylorModel{I}}}, points, B) where {I <: Interval}
    model = first(models)
    T = IA.numtype(I)
    projected = map(x -> DA.polygon_project(B, x), points)
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
    exact = point ? polygon_interval_point(a, query) : nothing
    constraints = query === nothing || point ? nothing : polygon_interval_constraints(a, query)
    result = nothing
    query_points = if point
        (exact,)
    elseif query === nothing
        a._domain._vertices
    elseif query isa DA.ConvexPolygon
        query._vertices
    else
        lo, hi = DA.polygon_real.(IA.inf.(query)), DA.polygon_real.(IA.sup.(query))
        ((lo[1], lo[2]), (lo[1], hi[2]), (hi[1], lo[2]), (hi[1], hi[2]))
    end
    DA.ads_visit_intersections(a.nodes, query_points) do index
        patch = a.patches[index]
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
        values = polygon_snapshot_bounds(patch._patch.models, points, a._directions)
        result = result === nothing ? values : [IA.hull(x, y; dec = :auto) for (x, y) in zip(result, values)]
        return nothing
    end
    result === nothing && throw(ErrorException("Polygon partition has a coverage gap"))
    return a._scalar ? only(result) : result
end
