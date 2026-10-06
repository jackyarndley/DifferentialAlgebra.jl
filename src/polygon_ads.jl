"""
    PolygonPatch

One local map patch on a physical convex polygon subdomain of a `PiecewisePolygonMap`.
`domain(patch)` owns a copy of its exact vertices.
`error_estimate`, `depth` and `status` describe the
selected ADS error method. `error_bounds` is available only for `IntervalBound`.
The internal box snapshot is in projected coordinates, not physical coordinates.
"""
struct PolygonPatch{T <: AbstractFloat, P}
    _domain::ConvexPolygon{T}
    _patch::P
end
domain(p::PolygonPatch) = copy(p._domain)
function Base.getproperty(p::PolygonPatch, name::Symbol)
    raw = getfield(p, :_patch)
    name in (:depth, :status) && return getproperty(raw, name)
    name === :error_estimate && return raw isa TaylorModelPatch ? deepcopy(raw.error_bounds) : copy(raw.error_estimate)
    name === :error_bounds && return raw isa TaylorModelPatch ? deepcopy(raw.error_bounds) : throw(ArgumentError("This patch has heuristic estimates, not certified bounds"))
    return getfield(p, name)
end
Base.propertynames(::PolygonPatch) = (:_domain, :_patch, :depth, :status, :error_estimate, :error_bounds)

"""
    PiecewisePolygonMap

Static ADS patches on an exact convex polygon subdomain partition. Physical point
queries enforce the polygon domain; shared edges are included. With `IntervalBound`, evaluation
and enclosure include the remainders. Other estimators return ordinary numeric
approximations. Snapshots survive algebra reinitialization. Public geometry
accessors and `copy` own their data; internal arrays are read-only.
"""
struct PiecewisePolygonMap{T <: AbstractFloat, P, E <: ADSEstimator}
    _domain::ConvexPolygon{T}
    _directions::Matrix{PolygonReal}
    patches::Tuple{Vararg{PolygonPatch{T, P}}}
    _order::Int
    _scalar::Bool
    _estimator::E
    converged::Bool
end
domain(a::PiecewisePolygonMap) = copy(a._domain)
"""
    split_directions(map::PiecewisePolygonMap)

Return an independent exact 2x2 matrix whose rows project physical coordinates
onto the splitting frame. It need not be orthonormal.
"""
split_directions(a::PiecewisePolygonMap) = deepcopy(a._directions)
nvariables(::PiecewisePolygonMap) = 2
noutputs(a::PiecewisePolygonMap) = polygon_noutputs(first(a.patches)._patch)
polygon_noutputs(p::TaylorPatch) = noutputs(p.map)
polygon_noutputs(p::TaylorModelPatch) = length(p.models)
degree(a::PiecewisePolygonMap) = maximum(p -> polygon_degree(p._patch), a.patches)
polygon_degree(p::TaylorPatch) = degree(p.map)
polygon_degree(p::TaylorModelPatch) = maximum(degree, p.models)
max_order(a::PiecewisePolygonMap) = a._order
Base.copy(a::PiecewisePolygonMap) = PiecewisePolygonMap(copy(a._domain), deepcopy(a._directions), deepcopy(a.patches), a._order, a._scalar, a._estimator, a.converged)
Base.deepcopy_internal(a::PiecewisePolygonMap, copies::IdDict) = get!(() -> copy(a), copies, a)
(a::PiecewisePolygonMap)(point) = evaluate(a, point)
Base.show(io::IO, a::PiecewisePolygonMap) = print(io, "PiecewisePolygonMap(", length(a.patches), " polygons, ", a.converged ? "converged" : "unresolved", ")")

"""
    adaptive_map(
        f, polygon::ConvexPolygon; estimator = GuardedTail(),
        splitter = :oriented, directions = nothing, kwargs...
    )
    adaptive_map(f, lower, upper; splitter = :oriented, kwargs...)

Split a 2D convex domain by straight lines in a linear coordinate frame. Exact
rational clipping preserves coverage without gaps. Automatic directions use
gradient/Hessian sensitivity of an interior degree-two probe; they are a
heuristic choice of frame, not an optimality claim. Supply a nonsingular 2x2
matrix with projection directions in its rows to override this choice. Use
`directions=:axes` for an unrotated polygon baseline. Each child reevaluates `f`
over its bounding parallelogram. The callback must be valid over these covers.

All four ADS error methods are supported. Only `IntervalBound()` supplies
certified uniform errors and requires IntervalArithmetic. `:oriented` uses
retained-coefficient sensitivity to choose a projection axis; `:width` uses
relative projected widths. Domain queries retain the physical polygon.
"""
function adaptive_map(f, polygon::ConvexPolygon; kwargs...)
    return polygon_construct(f, polygon; kwargs...)
end
function adaptive_map(
        f, previous::PiecewisePolygonMap;
        estimator::ADSEstimator = previous._estimator,
        order::Integer = max(1, degree(previous)),
        directions = split_directions(previous), kwargs...
    )
    return polygon_construct(f, domain(previous); estimator, order, directions, partition = previous, kwargs...)
end

function polygon_check_context(ctx, order)
    CURRENT_ALGEBRA[] === ctx && ctx.cutoff == order && iszero(ctx.epsilon) && ctx.big_epsilon === nothing ||
        throw(ArgumentError("The callback changed the algebra configuration"))
    return nothing
end

# Sensitivity only selects a frame. Neither midpoints nor these nearest-rounded
# scores enter a validated coefficient, remainder or acceptance calculation.
function polygon_sensitivity(polynomials, radii, scales, tolerances, ::Type{T}) where {T}
    vectors = NTuple{2, T}[]
    for (p, tolerance) in zip(polynomials, tolerances)
        value(alpha) = T(polygon_mid(coefficient(p, alpha)))
        g = (value([1, 0]) / radii[1], value([0, 1]) / radii[2])
        h11, h12, h22 = 2value([2, 0]) / radii[1]^2, value([1, 1]) / prod(radii), 2value([0, 2]) / radii[2]^2
        push!(vectors, (g[1] / tolerance, g[2] / tolerance))
        push!(vectors, (h11 * scales[1] / tolerance, h12 * scales[1] / tolerance))
        push!(vectors, (h12 * scales[2] / tolerance, h22 * scales[2] / tolerance))
    end
    all(v -> all(isfinite, v), vectors) || return Matrix{T}(I, 2, 2)
    scale = maximum(v -> max(abs(v[1]), abs(v[2])), vectors)
    iszero(scale) && return Matrix{T}(I, 2, 2)
    a = b = d = zero(T)
    for (x, y) in vectors
        x /= scale; y /= scale
        a += x^2; b += x * y; d += y^2
    end
    if iszero(b)
        return a >= d ? Matrix{T}(I, 2, 2) : T[0 1; -1 0]
    end
    delta = (a - d) / 2
    t = abs(delta) + hypot(delta, b)
    x, y = delta >= 0 ? (t, b) : (b, t)
    scale = max(abs(x), abs(y))
    x /= scale; y /= scale
    return T[x y; -y x]
end
polygon_mid(x::Real) = x

# A centered, strictly interior axis box for the direction probe. Its vertices
# satisfy every physical halfspace exactly, including for very thin polygons.
function polygon_probe(p::ConvexPolygon{T}) where {T}
    v = p._vertices
    center = ntuple(i -> sum(x -> x[i], v) / length(v), 2)
    lo, hi = polygon_projection_bounds(p, Matrix{PolygonReal}(I, 2, 2))
    radii = ((hi[1] - lo[1]) / 2, (hi[2] - lo[2]) / 2)
    factor = one(PolygonReal)
    for i in eachindex(v)
        a, b = v[i], v[mod1(i + 1, length(v))]
        normal = (b[2] - a[2], a[1] - b[1])
        slack = normal[1] * (a[1] - center[1]) + normal[2] * (a[2] - center[2])
        factor = min(factor, slack / (abs(normal[1]) * radii[1] + abs(normal[2]) * radii[2]))
    end
    r = (factor * radii[1] / 2, factor * radii[2] / 2)
    points = ((center[1] - r[1], center[2] - r[2]), (center[1] + r[1], center[2] - r[2]), (center[1] + r[1], center[2] + r[2]), (center[1] - r[1], center[2] + r[2]))
    return ConvexPolygon{T}(points, Val(:owned)), T.(radii)
end

function polygon_candidate(f, polygon, B, A, options, ctx, estimator::ADSEstimator)
    T = eltype(options.scales)
    lo, hi = polygon_projection_bounds(polygon, B)
    geometry = ads_geometry(T.(collect(lo)), T.(collect(hi)))
    all(r -> isfinite(r) && r > 0, geometry.radius) || throw(ArgumentError("Polygon cover is not representable; use a wider scalar type"))
    x = [variable(i, T) for i in 1:2]
    z = geometry.center .+ geometry.radius .* x
    values, scalar = ads_outputs(f(T.(A) * z))
    polygon_check_context(ctx, options.working_order)
    C = foldl((R, p) -> promote_type(R, p isa TaylorPolynomial ? coefficient_type(p) : typeof(p)), values; init = T)
    C <: AbstractFloat || throw(ArgumentError("Ordinary polygon ADS requires floating coefficients; select IntervalBound for interval validation"))
    polynomials = TaylorPolynomial{C}.(values)
    all(p -> valid(p) === ctx && all(isfinite, @view(p.coeffs[1:p.len])), polynomials) || throw(ArgumentError("Invalid Taylor coefficients"))
    errors = [ads_error(p, options.order, estimator) for p in polynomials]
    contributions = zeros(C, length(values), 2)
    for (j, p) in enumerate(polynomials)
        ads_contributions!(contributions, p, j, options.order, estimator)
    end
    compiled = CompiledMap([trim(p, 0, options.order) for p in polynomials])
    p = (; geometry..., compiled, errors, contributions, scale = abs.(constant_term.(polynomials)), scalar)
    if options.check_points
        points = collect(polygon._vertices)
        push!(points, ntuple(i -> sum(v -> v[i], polygon._vertices) / length(polygon._vertices), 2))
        probes = [(T.(collect(v)), (T.(collect(polygon_project(B, v))) - geometry.center) ./ geometry.radius, 0) for v in points]
        ads_check_points!(p, probes, (f(point) for (point, _, _) in probes))
    end
    polygon_check_context(ctx, options.working_order)
    assessed = ads_assess(p, options)
    directions = polygon_sensitivity(polynomials, geometry.radius, options.scales, assessed.tolerance, T)
    return (; assessed..., projection = (lo, hi), directions)
end
function polygon_candidate(f, polygon, B, A, options, ctx, estimator::IntervalBound)
    throw(ArgumentError("Load IntervalArithmetic to use IntervalBound on polygons"))
end

function polygon_axis(p, original, options)
    best, axis = -Inf, 0
    widest = -Inf
    for i in 1:2
        polygon_splittable(p, i) || continue
        width = (p.projection[2][i] - p.projection[1][i]) / (original[2][i] - original[1][i])
        score = options.splitter == :width ? zero(best) : maximum(j -> ads_ratio(p.contributions[j, i], p.tolerance[j]), eachindex(p.errors))
        if score > best || (score == best && width > widest)
            axis, best, widest = i, score, width
        end
    end
    return axis
end
polygon_splittable(p, i) = p.lo[i] < p.center[i] < p.hi[i]
polygon_patch(p, depth, status) = TaylorPatch(p.lo, p.hi, p.center, p.radius, p.compiled, p.errors, depth, status)

function polygon_construct(
        f, polygon::ConvexPolygon{T}; order::Integer = 5, atol = 1.0e-8, rtol::Real = 0,
        estimator::ADSEstimator = GuardedTail(), splitter::Symbol = :oriented,
        guard_order::Integer = estimator isa GuardedTail ? 2 : 0,
        check_points::Bool = !(estimator isa IntervalBound), strict::Bool = true,
        directions = nothing, max_depth::Integer = 20, max_patches::Integer = 1024,
        names = nothing, table_bytes::Integer = 32 * 1024^2, partition = nothing
    ) where {T}
    1 <= order <= 65535 && 0 <= guard_order <= 65535 - order || throw(ArgumentError("Invalid order"))
    (estimator isa GuardedTail ? guard_order > 0 : guard_order == 0) || throw(ArgumentError("Invalid guard_order for estimator"))
    splitter in (:oriented, :tail, :width) || throw(ArgumentError("Unsupported polygon splitter"))
    isfinite(rtol) && rtol >= 0 || throw(ArgumentError("Invalid relative tolerance"))
    0 <= max_depth <= typemax(Int) && 1 <= max_patches <= typemax(Int) || throw(ArgumentError("Invalid splitting limits"))
    estimator isa IntervalBound && (rtol != 0 || check_points) && throw(ArgumentError("IntervalBound uses absolute bounds, not relative tolerances or samples"))
    partition !== nothing && (length(partition.patches) > max_patches || maximum(p -> p.depth, partition.patches) > max_depth) && throw(ArgumentError("Limits are smaller than the existing partition"))
    _, scales = polygon_probe(polygon)
    options = (; order = Int(order), working_order = Int(order + guard_order), atol, rtol, estimator, splitter, check_points, strict, scales)
    return with_algebra(max(2, options.working_order), 2; names, table_bytes) do ctx
        B, A = if directions === nothing || directions === :auto
            probe, _ = polygon_probe(polygon)
            probe_options = (; options..., order = 2, working_order = 2, check_points = false)
            candidate = with_order(2) do
                polygon_candidate(f, probe, Matrix{PolygonReal}(I, 2, 2), Matrix{PolygonReal}(I, 2, 2), probe_options, ctx, estimator)
            end
            polygon_frame(candidate.directions)
        else
            polygon_frame(directions === :axes ? Matrix{T}(I, 2, 2) : directions)
        end
        original = polygon_projection_bounds(polygon, B)
        pending = partition === nothing ? [(copy(polygon), 0)] : [(domain(p), p.depth) for p in reverse(partition.patches)]
        patch_count = length(pending)
        patches = nothing
        shape = nothing
        with_order(options.working_order) do
            while !isempty(pending)
                child, depth = pop!(pending)
                p = polygon_candidate(f, child, B, A, options, ctx, estimator)
                newshape = (p.scalar, length(p.errors))
                shape === nothing ? (shape = newshape) : shape == newshape || throw(DimensionMismatch("Callback changed output shape"))
                axis = p.accepted ? 0 : polygon_axis(p, original, options)
                status = p.accepted ? :converged : depth >= max_depth ? :max_depth : patch_count >= max_patches ? :max_patches : axis == 0 ? :roundoff : :split
                if status == :split
                    cut = (p.projection[1][axis] + p.projection[2][axis]) / 2
                    normal = (B[axis, 1], B[axis, 2])
                    left = polygon_clip(child, normal, cut)
                    right = polygon_clip(child, (-normal[1], -normal[2]), -cut)
                    left === nothing || right === nothing ? throw(ArgumentError("Degenerate polygon split")) : nothing
                    push!(pending, (right, depth + 1), (left, depth + 1))
                    patch_count += 1
                else
                    strict && !p.accepted && throw(ErrorException("Polygon ADS reached $status; increase limits or use strict=false"))
                    patch = PolygonPatch(child, polygon_patch(p, depth, status))
                    patches === nothing && (patches = typeof(patch)[])
                    patch isa eltype(patches) || throw(ArgumentError("Callback changed coefficient type"))
                    push!(patches, patch)
                end
            end
        end
        return PiecewisePolygonMap(copy(polygon), deepcopy(B), Tuple(patches), Int(order), shape[1], estimator, all(p -> p.status == :converged, patches))
    end
end

function polygon_point(a, point)
    point isa Union{Tuple, AbstractVector} && length(point) == 2 || throw(DimensionMismatch("Provide two physical coordinates"))
    x = (polygon_real(point[1]), polygon_real(point[2]))
    polygon_contains(a._domain, x) || throw(DomainError(point, "Point lies outside the polygon domain"))
    return x
end
function evaluate(a::PiecewisePolygonMap, point)
    a._estimator isa IntervalBound && return polygon_interval_query(a, point; point = true)
    x = polygon_point(a, point)
    patch = findfirst(p -> polygon_contains(p._domain, x), a.patches)
    patch === nothing && throw(ErrorException("Polygon partition has a coverage gap"))
    p = a.patches[patch]._patch
    T = eltype(p.center)
    z = (T.(collect(polygon_project(a._directions, x))) - p.center) ./ p.radius
    result = evaluate(p.map, z)
    return a._scalar ? only(result) : result
end
function enclose(a::PiecewisePolygonMap, args...)
    a._estimator isa IntervalBound || throw(ArgumentError("Select IntervalBound for a certified function enclosure"))
    return polygon_interval_query(a, args...)
end
function polygon_interval_query end
polygon_interval_query(args...; kwargs...) = throw(ArgumentError("Load IntervalArithmetic for polygon enclosures"))
