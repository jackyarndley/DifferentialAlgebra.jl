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
    if name in (:error_estimate, :error_bounds)
        raw isa TaylorModelPatch && return raw.error_bounds
        name === :error_bounds && throw(
            ArgumentError("This patch has heuristic estimates, not certified bounds")
        )
        return copy(raw.error_estimate)
    end
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
    axis::Int
    bounds::Tuple{Float64, Float64}
end
PolygonNode(normal, midpoint, left, right, patch) = PolygonNode(
    normal, midpoint, left, right, patch, 0,
    (Float64(midpoint, RoundDown), Float64(midpoint, RoundUp))
)
PolygonNode() = PolygonNode((zero(PolygonReal), zero(PolygonReal)), zero(PolygonReal), 0, 0, 0)

# Current-frame rows come first; refinement may retain cuts in older frames.
# Nodes store a row index, so each distinct projection is bounded once per query.
struct PolygonLookup
    directions::Matrix{Float64}
    bounds::Tuple{Matrix{Float64}, Matrix{Float64}}
    domain::Tuple{Matrix{Float64}, Matrix{Float64}}
end
function polygon_lookup(polygon, B, nodes)
    normals = PolygonPoint[(B[i, 1], B[i, 2]) for i in 1:2]
    rows = Dict(normal => i for (i, normal) in enumerate(normals))
    for (i, node) in enumerate(nodes)
        node.patch != 0 && continue
        axis = get(rows, node.normal, 0)
        if axis == 0
            push!(normals, node.normal)
            axis = length(normals)
            rows[node.normal] = axis
        end
        nodes[i] = PolygonNode(
            node.normal, node.midpoint, node.left, node.right, node.patch,
            axis, node.bounds
        )
    end
    projections = [normal[j] for normal in normals, j in 1:2]
    return PolygonLookup(
        Float64.(B), ads_float_bounds(projections), polygon_filter(polygon)
    )
end

mutable struct PolygonQuery{P}
    point::P
    bounds::Tuple{Vector{Float64}, Vector{Float64}}
    exact::Union{Nothing, PolygonPoint}
end
function polygon_exact!(query::PolygonQuery)
    query.exact === nothing && (query.exact = polygon_real.(query.point))
    return query.exact::PolygonPoint
end
ads_query_bounds(node::PolygonNode, query::PolygonQuery) =
    (query.bounds[1][node.axis], query.bounds[2][node.axis])
function ads_visit_left(node::PolygonNode, query::PolygonQuery, lo, hi)
    hi <= node.bounds[1] && return true
    lo > node.bounds[2] && return false
    return ads_node_value(node, polygon_exact!(query)) <= node.midpoint
end
function ads_visit_right(node::PolygonNode, query::PolygonQuery, lo, hi)
    lo >= node.bounds[2] && return true
    hi < node.bounds[1] && return false
    return ads_node_value(node, polygon_exact!(query)) >= node.midpoint
end
ads_point_left(node::PolygonNode, query::PolygonQuery) =
    ads_visit_left(node, query, ads_query_bounds(node, query)...)

"""
    PiecewisePolygonMap

Static ADS patches on an exact convex polygon subdomain partition. Physical point
queries enforce the polygon domain; shared edges are included. With `IntervalBound`, evaluation
and enclosure include the remainders. Other estimators return ordinary numeric
approximations. Snapshots survive algebra reinitialization. Public geometry
accessors and `copy` own their data; internal arrays are read-only.
"""
struct PiecewisePolygonMap{T <: AbstractFloat, P, E <: ADSEstimator, O <: NamedTuple}
    _domain::ConvexPolygon{T}
    _directions::Matrix{PolygonReal}
    _lookup::PolygonLookup
    patches::Vector{PolygonPatch{T, P}}
    nodes::Vector{PolygonNode}
    _order::Int
    _scalar::Bool
    _estimator::E
    _options::O
    converged::Bool
end
Base.getproperty(a::PiecewisePolygonMap, name::Symbol) =
    name === :patches ? copy(ads_patches(a)) : getfield(a, name)
domain(a::PiecewisePolygonMap) = copy(a._domain)
"""
    split_directions(map::PiecewisePolygonMap)

Return an independent exact 2x2 matrix whose rows project physical coordinates
onto the splitting frame. It need not be orthonormal.
"""
split_directions(a::PiecewisePolygonMap) = deepcopy(a._directions)
nvariables(::PiecewisePolygonMap) = 2
noutputs(a::PiecewisePolygonMap) = noutputs(first(ads_patches(a)))
degree(a::PiecewisePolygonMap) = maximum(degree, ads_patches(a))
max_order(a::PiecewisePolygonMap) = a._order
Base.copy(a::PiecewisePolygonMap) = PiecewisePolygonMap(
    copy(a._domain), deepcopy(a._directions), deepcopy(a._lookup),
    deepcopy(ads_patches(a)), deepcopy(a.nodes),
    a._order, a._scalar, a._estimator, deepcopy(a._options), a.converged
)
Base.deepcopy_internal(a::PiecewisePolygonMap, copies::IdDict) = get!(() -> copy(a), copies, a)
(a::PiecewisePolygonMap)(point) = evaluate(a, point)
Base.show(io::IO, a::PiecewisePolygonMap) = print(
    io, "PiecewisePolygonMap(", length(ads_patches(a)), " polygons, ",
    a.converged ? "converged" : "unresolved", ")"
)

"""
    adaptive_map(
        f, polygon::ConvexPolygon; estimator = GuardedTail(),
        splitter = :oriented, directions = nothing, kwargs...
    )
    adaptive_map(f, lower, upper; splitter = :oriented, kwargs...)

Split a 2D convex domain by straight lines in a linear coordinate frame. Exact
rational clipping preserves coverage without gaps. Automatic directions use
nonlinear sensitivity of one interior probe of degree three to five; they are a
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
        directions = split_directions(previous), kwargs...
    )
    options = ads_refinement_options(previous, estimator; kwargs...)
    return polygon_construct(
        f, domain(previous); estimator, directions, partition = previous, options...
    )
end

# Sensitivity only selects a frame. Neither midpoints nor these nearest-rounded
# scores enter a validated coefficient, remainder or acceptance calculation.
function polygon_sensitivity(polynomials, radii, scales, tolerances, ::Type{T}) where {T}
    vectors = NTuple{2, T}[]
    powers = zeros(Int, 2)
    for (p, tolerance) in zip(polynomials, tolerances)
        ratios = scales ./ radii
        # Group gradient coefficients by their residual multiindex. Excluding
        # the constant gradient removes affine variation. Scale the probe's
        # normalized coordinates to the root span before comparing outputs.
        for d in 1:(degree(p) - 1), j in 0:d
            powers .= (j + 1, d - j)
            x = (j + 1) * T(polygon_mid(coefficient(p, powers)))
            powers .= (j, d - j + 1)
            y = (d - j + 1) * T(polygon_mid(coefficient(p, powers)))
            iszero(x) && iszero(y) && continue
            weight = ratios[1]^j * ratios[2]^(d - j) / tolerance
            push!(vectors, (weight * (x / radii[1]), weight * (y / radii[2])))
        end
    end
    isempty(vectors) && return Matrix{T}(I, 2, 2)
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
        normal, cut = polygon_halfspace(p, i)
        slack = cut - normal[1] * center[1] - normal[2] * center[2]
        factor = min(factor, slack / (abs(normal[1]) * radii[1] + abs(normal[2]) * radii[2]))
    end
    r = (factor * radii[1] / 2, factor * radii[2] / 2)
    points = ((center[1] - r[1], center[2] - r[2]), (center[1] + r[1], center[2] - r[2]), (center[1] + r[1], center[2] + r[2]), (center[1] - r[1], center[2] + r[2]))
    return ConvexPolygon{T}(points, Val(:owned)), T.(radii)
end

function polygon_candidate(f, polygon, B, A, options, ctx, estimator::ADSEstimator; probe::Bool = false)
    T = eltype(options.scales)
    lo, hi = polygon_projection_bounds(polygon, B)
    lower = [T(v, RoundDown) for v in lo]
    upper = [T(v, RoundUp) for v in hi]
    geometry = ads_geometry(lower, upper)
    all(r -> isfinite(r) && r >= 0, geometry.radius) || throw(
        ArgumentError(
            "Polygon cover is not representable; use a wider scalar type"
        )
    )
    x = [variable(i, T) for i in 1:2]
    z = geometry.center .+ geometry.radius .* x
    p = ads_analyze(f(T.(A) * z), geometry, options, ctx)
    if options.check_points
        points = collect(polygon._vertices)
        push!(points, ntuple(i -> sum(v -> v[i], polygon._vertices) / length(polygon._vertices), 2))
        probes = map(points) do v
            projected = T.(collect(polygon_project(B, v)))
            normalized = ads_normalize.(projected, geometry.center, geometry.radius)
            (T.(collect(v)), normalized, 0)
        end
        ads_check_points!(p, probes, (f(point) for (point, _, _) in probes))
    end
    ads_check_context(ctx, options.working_order)
    assessed = ads_assess(p, options)
    directions = probe ? polygon_sensitivity(p.polynomials, geometry.radius, options.scales, assessed.tolerance, T) : nothing
    cover = (; geometry..., projection = (lo, hi), directions)
    return ADSCandidate(cover, assessed.payload, assessed.errors, assessed.contributions, assessed.tolerance, assessed.scalar, assessed.accepted)
end
function polygon_candidate(f, polygon, B, A, options, ctx, estimator::IntervalBound; probe::Bool = false)
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
    lookup = polygon_lookup(g.root, g.B, nodes)
    return PiecewisePolygonMap(
        copy(g.root), deepcopy(g.B), lookup, patches, nodes, options.order, scalar,
        options.estimator, ads_saved_options(options),
        all(p -> p.status == :converged, patches)
    )
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
        names = nothing, table_bytes::Integer = 32 * 1024^2, partition = nothing,
        interval_domain = ads_interval_domain(partition)
    ) where {T}
    settings = ads_options(;
        order, atol, rtol, estimator, splitter, guard_order,
        max_depth, max_patches, check_points, strict,
        splitters = (:oriented, :tail, :width), names, table_bytes
    )
    ads_check_partition(partition, settings)
    probe, scales = polygon_probe(polygon)
    options = (; settings..., scales, interval_domain)
    automatic = directions === nothing || directions === :auto
    probe_order = clamp(options.working_order, 3, 5)
    basis_order = max(options.working_order, automatic ? probe_order : 2)
    return with_algebra(basis_order, 2; names, table_bytes) do ctx
        B, A = if automatic
            probe_options = (;
                options..., order = probe_order, working_order = probe_order,
                check_points = false,
            )
            candidate = with_order(probe_order) do
                polygon_candidate(f, probe, Matrix{PolygonReal}(I, 2, 2), Matrix{PolygonReal}(I, 2, 2), probe_options, ctx, estimator; probe = true)
            end
            frame = candidate isa ADSDomainFailure ? Matrix{T}(I, 2, 2) :
                candidate.geometry.directions
            polygon_frame(frame)
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

function polygon_query(a, point)
    point isa Union{Tuple, AbstractVector} && length(point) == 2 || throw(DimensionMismatch("Provide two physical coordinates"))
    point isa AbstractVector && Base.require_one_based_indexing(point)
    x = (
        ads_primal_coordinate(point[1]), ads_primal_coordinate(point[2]),
    )
    all(isfinite, x) || throw(ArgumentError("Physical coordinates must be finite"))
    bounds = ads_float_bounds(x)
    polygon_contains_filtered(a._domain, x, bounds, a._lookup.domain) ||
        throw(DomainError(point, "Point lies outside the polygon domain"))
    return PolygonQuery(x, ads_project_bounds(a._lookup.bounds, bounds), nothing)
end
function evaluate(a::PiecewisePolygonMap, point)
    a._estimator isa IntervalBound && return polygon_interval_query(a, point; point = true)
    x = polygon_query(a, point)
    patch = ads_point_patch(a.nodes, x)
    p = ads_patches(a)[patch]._patch
    T = evaluation_type(eltype(p.center), point)
    coordinates = T[point[1], point[2]]
    B = coefficient_primal(first(coordinates)) isa Float64 ?
        a._lookup.directions : T.(a._directions)
    projected = B * coordinates
    z = ads_normalize.(projected, p.center, p.radius)
    result = evaluate(p.map, z)
    return a._scalar ? only(result) : result
end
function enclose(a::PiecewisePolygonMap, args...)
    a._estimator isa IntervalBound || throw(ArgumentError("Select IntervalBound for a certified function enclosure"))
    return polygon_interval_query(a, args...)
end
function polygon_interval_query end
polygon_interval_query(args...; kwargs...) = throw(ArgumentError("Load IntervalArithmetic for polygon enclosures"))
