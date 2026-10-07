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
nvariables(::PolygonPatch) = 2
noutputs(p::PolygonPatch) = noutputs(p._patch)
degree(p::PolygonPatch) = degree(p._patch)
max_order(p::PolygonPatch) = max_order(p._patch)

# Cuts remain in physical coordinates, so refinement can choose a new frame
# while preserving all of the old partition's boundaries and lookup branches.
struct PolygonNode
    normal::PolygonPoint
    midpoint::PolygonReal
    left::Int
    right::Int
    patch::Int
end
PolygonNode() = PolygonNode((zero(PolygonReal), zero(PolygonReal)), zero(PolygonReal), 0, 0, 0)

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
    nodes::Vector{PolygonNode}
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
noutputs(a::PiecewisePolygonMap) = noutputs(first(a.patches))
degree(a::PiecewisePolygonMap) = maximum(degree, a.patches)
max_order(a::PiecewisePolygonMap) = a._order
Base.copy(a::PiecewisePolygonMap) = PiecewisePolygonMap(copy(a._domain), deepcopy(a._directions), deepcopy(a.patches), deepcopy(a.nodes), a._order, a._scalar, a._estimator, a.converged)
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
        order::Integer = max_order(previous),
        directions = split_directions(previous), kwargs...
    )
    return polygon_construct(f, domain(previous); estimator, order, directions, partition = previous, kwargs...)
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
    p = ads_analyze(f(T.(A) * z), geometry, options, ctx)
    if options.check_points
        points = collect(polygon._vertices)
        push!(points, ntuple(i -> sum(v -> v[i], polygon._vertices) / length(polygon._vertices), 2))
        probes = [(T.(collect(v)), (T.(collect(polygon_project(B, v))) - geometry.center) ./ geometry.radius, 0) for v in points]
        ads_check_points!(p, probes, (f(point) for (point, _, _) in probes))
    end
    ads_check_context(ctx, options.working_order)
    assessed = ads_assess(p, options)
    directions = polygon_sensitivity(p.polynomials, geometry.radius, options.scales, assessed.tolerance, T)
    cover = (; geometry..., projection = (lo, hi), directions)
    return ADSCandidate(cover, assessed.payload, assessed.errors, assessed.contributions, assessed.tolerance, assessed.scalar, assessed.accepted)
end
function polygon_candidate(f, polygon, B, A, options, ctx, estimator::IntervalBound)
    throw(ArgumentError("Load IntervalArithmetic to use IntervalBound on polygons"))
end

struct PolygonADS{T <: AbstractFloat}
    root::ConvexPolygon{T}
    B::Matrix{PolygonReal}
    A::Matrix{PolygonReal}
    original::Tuple{NTuple{2, PolygonReal}, NTuple{2, PolygonReal}}
end
ads_root_node(::PolygonADS) = PolygonNode()
ads_patch_domain(::PolygonADS, p::PolygonPatch) = domain(p)
ads_leaf(::PolygonNode, patch) = PolygonNode((zero(PolygonReal), zero(PolygonReal)), zero(PolygonReal), 0, 0, patch)
ads_node_value(node::PolygonNode, point) = node.normal[1] * point[1] + node.normal[2] * point[2]
function ads_query_bounds(node::PolygonNode, points)
    values = map(x -> ads_node_value(node, x), points)
    return minimum(values), maximum(values)
end
# Positive-area exact polygons can be clipped; representability of the fitted
# cover is checked separately before choosing a cut.
ads_has_split(::ConvexPolygon) = true
ads_fit(problem::StaticMap, g::PolygonADS, polygon, ctx, options, allow_split) =
    polygon_candidate(problem.f, polygon, g.B, g.A, options, ctx, options.estimator)
ads_patch_type(::PolygonADS{T}, p) where {T} = PolygonPatch{T, ads_payload_patch_type(p.payload, p.geometry)}
function ads_finalize(::PolygonADS, polygon, p, depth, status, order)
    return PolygonPatch(polygon, ads_local_patch(p.payload, p.geometry, p.errors, depth, status, order))
end
function ads_result(g::PolygonADS, patches, nodes, scalar, options)
    return PiecewisePolygonMap(copy(g.root), deepcopy(g.B), Tuple(patches), nodes, options.order, scalar, options.estimator, all(p -> p.status == :converged, patches))
end

function ads_relative_width(g::PolygonADS, cover, i)
    projection = cover.projection
    return (projection[2][i] - projection[1][i]) / (g.original[2][i] - g.original[1][i])
end
function ads_split(g::PolygonADS, polygon, p, axis, left, right)
    projection = p.geometry.projection
    cut = (projection[1][axis] + projection[2][axis]) / 2
    normal = (g.B[axis, 1], g.B[axis, 2])
    left_polygon = polygon_clip(polygon, normal, cut)
    right_polygon = polygon_clip(polygon, (-normal[1], -normal[2]), -cut)
    left_polygon === nothing || right_polygon === nothing ? throw(ArgumentError("Degenerate polygon split")) : nothing
    return left_polygon, right_polygon, PolygonNode(normal, cut, left, right, 0)
end

function polygon_construct(
        f, polygon::ConvexPolygon{T}; order::Integer = 5, atol = 1.0e-8, rtol::Real = 0,
        estimator::ADSEstimator = GuardedTail(), splitter::Symbol = :oriented,
        guard_order::Integer = estimator isa GuardedTail ? 2 : 0,
        check_points::Bool = !(estimator isa IntervalBound), strict::Bool = true,
        directions = nothing, max_depth::Integer = 20, max_patches::Integer = 1024,
        names = nothing, table_bytes::Integer = 32 * 1024^2, partition = nothing
    ) where {T}
    settings = ads_options(;
        order, atol, rtol, estimator, splitter, guard_order,
        max_depth, max_patches, check_points, strict, splitters = (:oriented, :tail, :width)
    )
    ads_check_partition(partition, settings)
    _, scales = polygon_probe(polygon)
    options = (; settings..., scales)
    return with_algebra(max(2, options.working_order), 2; names, table_bytes) do ctx
        B, A = if directions === nothing || directions === :auto
            probe, _ = polygon_probe(polygon)
            probe_options = (; options..., order = 2, working_order = 2, check_points = false)
            candidate = with_order(2) do
                polygon_candidate(f, probe, Matrix{PolygonReal}(I, 2, 2), Matrix{PolygonReal}(I, 2, 2), probe_options, ctx, estimator)
            end
            polygon_frame(candidate.geometry.directions)
        else
            polygon_frame(directions === :axes ? Matrix{T}(I, 2, 2) : directions)
        end
        original = polygon_projection_bounds(polygon, B)
        geometry = PolygonADS(copy(polygon), B, A, original)
        return with_order(options.working_order) do
            ads_build(StaticMap(f), geometry, ctx, options, partition)
        end
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
    patch = ads_point_patch(a.nodes, x)
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
