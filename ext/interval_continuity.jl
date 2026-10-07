# Included in the optional interval extension. Midpoint polynomials define only
# the explicitly requested smooth surrogate; certificates retain all widths/R.
function DA.continuity_layout(a::DA.PiecewiseTaylorModel)
    box = DA.domain(a)
    lo, hi = DA.polygon_real.(IA.inf.(box)), DA.polygon_real.(IA.sup.(box))
    B = Matrix{DA.PolygonReal}(DA.I, length(lo), length(lo))
    cores = [(DA.polygon_real.(IA.inf.(DA.domain(p))), DA.polygon_real.(IA.sup.(DA.domain(p)))) for p in a.patches]
    return (domain = box, bounds = (lo, hi), B, A = copy(B), cores, scalar = a._scalar)
end
DA.continuity_scalar_type(::DA.PiecewiseTaylorModel{I}) where {I <: Interval} = IA.numtype(I)

function DA.continuity_candidate(f, lo, hi, A, options, ctx, ::DA.IntervalBound)
    validated_rounding()
    T = options.T
    box = [checked_interval(IA.interval(T, a, b); guaranteed = true) for (a, b) in zip(lo, hi)]
    z = coordinate_models(box, options.order)
    x = [sum(IA.interval(T, A[i, j]) * z[j] for j in eachindex(z)) for i in eachindex(z)]
    fit = validated_fit(f, x, first(z), ctx, ctx.basis.order)
    values = fit.payload.models
    mids = [DA.TaylorPolynomial{T}(IA.mid.(m._polynomial.coeffs[1:m._polynomial.len]), m._polynomial.len, ctx) for m in values]
    map = DA.compile(mids)
    coordinates = first(z)._coordinates
    return (;
        center = IA.inf.(collect(coordinates.center)), radius = IA.inf.(collect(coordinates.radius)),
        map, certificates = compile_models(values),
        errors = fit.errors, scalar = fit.scalar,
    )
end
function DA.continuity_errors(patches::AbstractVector{<:DA.ContinuityPatch{T, C, V, E}}) where {T, C, V, E <: Tuple}
    return [foldl((a, p) -> IA.hull(a, p.errors[j]; dec = :auto), Iterators.drop(patches, 1); init = first(patches).errors[j]) for j in eachindex(first(patches).errors)]
end
function DA.continuity_check_tolerance(errors::AbstractVector{I}, atol) where {I <: Interval}
    atol === nothing && return nothing
    tol = validated_tolerances(I, atol, length(errors))
    all(j -> IA.sup(abs(errors[j])) <= tol[j], eachindex(tol)) || throw(ArgumentError("Overlap certificates exceed atol; refine the source map, increase order or reduce overlap"))
    return nothing
end

function continuity_query_bounds(a, query)
    if query === nothing
        return a._bounds
    elseif query isa DA.ConvexPolygon
        size(a._directions) == (2, 2) || throw(DimensionMismatch("Polygon queries require two variables"))
        v = query._vertices
        if a._domain isa DA.ConvexPolygon
            all(x -> DA.polygon_contains(a._domain, x), v) || throw(DomainError(query, "Query polygon lies outside the domain"))
        else
            all(x -> all(i -> a._bounds[1][i] <= x[i] <= a._bounds[2][i], eachindex(x)), v) || throw(DomainError(query, "Query polygon lies outside the domain"))
        end
        return DA.polygon_projection_bounds(query, a._directions)
    end
    query isa Union{Tuple, AbstractVector} && length(query) == DA.nvariables(a) || throw(DimensionMismatch("Provide one coordinate/interval per variable"))
    query isa AbstractVector && Base.require_one_based_indexing(query)
    if all(x -> x isa SupportedScalar, query)
        z = DA.continuity_point(a, query)
        return z, z
    end
    interval_box(query, length(query); guaranteed = true)
    lo, hi = DA.polygon_real.(IA.inf.(query)), DA.polygon_real.(IA.sup.(query))
    if a._domain isa DA.ConvexPolygon
        corners = ((lo[1], lo[2]), (lo[1], hi[2]), (hi[1], lo[2]), (hi[1], hi[2]))
        all(x -> DA.polygon_contains(a._domain, x), corners) || throw(DomainError(query, "Query box lies outside the polygon domain"))
    else
        all(i -> a._bounds[1][i] <= lo[i] <= hi[i] <= a._bounds[2][i], eachindex(lo)) || throw(DomainError(query, "Query box lies outside the domain"))
    end
    # Exact projection of an axis box under arbitrary signs; no corner explosion.
    B = a._directions
    lower = [sum(min(B[i, j] * lo[j], B[i, j] * hi[j]) for j in eachindex(lo)) for i in eachindex(lo)]
    upper = [sum(max(B[i, j] * lo[j], B[i, j] * hi[j]) for j in eachindex(lo)) for i in eachindex(lo)]
    return lower, upper
end

function DA.enclose(a::DA.ContinuousTaylorMap, query = nothing)
    validated_rounding()
    first(a._patches).certificates !== nothing || throw(ArgumentError("Use an IntervalBound source for certified function enclosures"))
    lo, hi = continuity_query_bounds(a, query)
    result = nothing
    for patch in a._patches
        lower, upper = max.(lo, patch.lower), min.(hi, patch.upper)
        any(lower .> upper) && continue
        # The query/support intersection is conservatively bounded in projected
        # coordinates. Certificate snapshots include the entire remainder.
        m = first(patch.certificates)
        T = IA.numtype(typeof(m._remainder))
        box = [IA.intersect_interval(IA.interval(T, l, h), b; dec = :auto) for (l, h, b) in zip(lower, upper, m._coordinates.box)]
        values = snapshot_bounds(patch.certificates, box)
        result = result === nothing ? values : map((x, y) -> IA.hull(x, y; dec = :auto), result, values)
    end
    result === nothing && throw(ErrorException("Overlap fits did not cover the query"))
    return a._scalar ? only(result) : collect(result)
end
