"""
    TaylorPatch

One box in a [`PiecewiseTaylorMap`](@ref). `lower` and `upper` are physical
coordinate bounds; `map` is a [`CompiledMap`](@ref) on normalized coordinates
`(point - center) ./ radius`. A fixed coordinate has normalized value zero.
`error_estimate` contains one absolute error indicator per output (the maximum
over checkpoints for a flow). Its meaning depends on the selected estimator.
`status` is `:converged`, `:max_depth`, `:max_patches`, or `:roundoff`.
The estimates are heuristic, not certified error bounds. Treat patch data as read-only.
"""
struct TaylorPatch{T <: AbstractFloat, C <: Real}
    lower::Vector{T}
    upper::Vector{T}
    center::Vector{T}
    radius::Vector{T}
    map::CompiledMap{C}
    error_estimate::Vector{C}
    depth::Int
    status::Symbol
end

struct DomainNode{T}
    axis::Int
    midpoint::T
    left::Int
    right::Int
    patch::Int
end
DomainNode(::Type{T}) where {T} = DomainNode(0, zero(T), 0, 0, 0)

"""
    PiecewiseTaylorMap

A piecewise Taylor approximation built by [`adaptive_map`](@ref).
Call `map(point)` with physical coordinates to select and evaluate a patch.
`patches` holds its [`TaylorPatch`](@ref)s; `converged` is true only if every
patch meets the requested error estimate. `lower` and `upper` bound the domain.
Numeric evaluation remains valid after the global algebra is reinitialized.
Shared faces belong to the lower child of a split. Extrapolation is rejected.
"""
struct PiecewiseTaylorMap{T <: AbstractFloat, C <: Real}
    lower::Vector{T}
    upper::Vector{T}
    patches::Vector{TaylorPatch{T, C}}
    nodes::Vector{DomainNode{T}}
    max_degree::Int
    scalar::Bool
    converged::Bool
end
nvariables(map::PiecewiseTaylorMap) = length(map.lower)
noutputs(map::PiecewiseTaylorMap) = noutputs(first(map.patches).map)
degree(map::PiecewiseTaylorMap) = map.max_degree

function Base.show(io::IO, map::PiecewiseTaylorMap)
    return print(
        io, "PiecewiseTaylorMap(", length(map.patches), " patches, ",
        noutputs(map), " outputs, ", nvariables(map), " variables, ",
        map.converged ? "converged" : "unresolved", ")"
    )
end

"""
    adaptive_map(
        f, lower, upper; order = 5, atol = 1.0e-8, rtol = 0,
        estimator = GuardedTail(), splitter = :tail,
        guard_order = estimator isa GuardedTail ? 2 : 0,
        max_depth = 20, max_patches = 1024,
        check_points = !(estimator isa IntervalBound), strict = true, names = nothing,
        directions = nothing,
        table_bytes = 32 * 1024^2
    )

Approximate `f` on a box by automatic domain splitting (ADS). `f` accepts a
vector of coordinates and returns a real scalar or a nonempty vector of reals.
It must support both ordinary numbers and [`TaylorPolynomial`](@ref)s and be
deterministic. Bounds are finite, ordered vectors; equal bounds fix a coordinate.
Their promoted floating-point type determines the input coefficient type.

Choose [`GuardedTail`](@ref) (extra degrees), [`ExtrapolatedTail`](@ref)
(coefficient decay), or [`LastTerms`](@ref) (retained tail). Only `GuardedTail`
uses a positive `guard_order`. Point checks compare the map with `f` at axis
endpoints and selected corners. `splitter = :tail` maximizes tolerance-scaled
tail reduction, with point checks and relative widths as fallbacks;
`:width` bisects the longest side relative to the original box.
Both children recompute `f` on their own domains.

Each output must satisfy `error_estimate ≤ atol + rtol * scale`, where `scale`
is the largest absolute value at the center and checked points. `atol` may be
a scalar or one tolerance per output. These are **heuristic estimates**, not
uniform error guarantees; validate on independent points. Disabling
`check_points` can miss terms beyond the computed order entirely.

The returned [`PiecewiseTaylorMap`](@ref) covers the whole input box. By default,
resource limits or unsplittable boxes throw an error. With `strict = false`,
unresolved patches are retained with their status and `converged = false`.
`max_depth` counts splits along a path; `max_patches` bounds the number of leaves.

Construction temporarily uses a separate global algebra and restores the
caller's configuration and polynomials, including on failure. Do not change
the algebra or its settings inside `f`, or construct maps concurrently with
other polynomial calculations. Evaluation of completed maps is independent.

`estimator=IntervalBound()` selects the optional validated static-map method.
Its callback accepts Taylor models, and its returned `PiecewiseTaylorModel`
retains interval coefficients, remainder and domain. It requires positive
absolute tolerances, zero relative tolerance, no guard degrees or point checks.
An explicit interval box can replace the two endpoint vectors. The other
estimators preserve their ordinary polynomial behavior.

`splitter=:oriented` selects exact convex polygon geometry in two dimensions,
with automatic or supplied 2x2 projection rows in `directions`. It works with
all four estimators and returns `PiecewisePolygonMap`; see the polygon overload.
"""
function adaptive_map(f, lower::AbstractVector{<:Real}, upper::AbstractVector{<:Real}; kwargs...)
    return ads_construct(StaticMap(f), lower, upper; kwargs...)
end

"""
    adaptive_map(f, previous::PiecewiseTaylorMap; order = max(1, degree(previous)), kwargs...)

Refine an existing partition using the same options as `adaptive_map(f, lower,
upper)`. Each leaf is reevaluated with `f`; the input map is unchanged. Existing
boundaries are retained, and `max_depth` and `max_patches` apply to the entire
tree. Supply the original function, not `previous` as a surrogate: subdividing
an already truncated polynomial cannot recover missing information.
"""
function adaptive_map(f, previous::PiecewiseTaylorMap; order::Integer = max(1, degree(previous)), kwargs...)
    return ads_construct(StaticMap(f), previous.lower, previous.upper, previous; order, kwargs...)
end

struct StaticMap{F}
    f::F
end

function ads_construct(
        problem, lower::AbstractVector{<:Real}, upper::AbstractVector{<:Real}, partition = nothing;
        order::Integer = 5, atol = 1.0e-8, rtol::Real = 0,
        estimator::ADSEstimator = GuardedTail(), splitter::Symbol = :tail,
        guard_order::Integer = estimator isa GuardedTail ? 2 : 0,
        max_depth::Integer = 20, max_patches::Integer = 1024,
        check_points::Bool = !(estimator isa IntervalBound), strict::Bool = true, names = nothing,
        directions = nothing,
        table_bytes::Integer = 32 * 1024^2
    )
    Base.require_one_based_indexing(lower, upper)
    length(lower) == length(upper) || throw(DimensionMismatch("Domain bounds differ in length"))
    isempty(lower) && throw(ArgumentError("The domain must have at least one coordinate"))
    if estimator isa IntervalBound
        problem isa StaticMap || throw(ArgumentError("IntervalBound supports static maps, not time integration"))
        partition === nothing || throw(ArgumentError("Construct a fresh interval map from the original function"))
        return interval_ads(
            problem.f, lower, upper; order, atol, rtol, guard_order, splitter,
            check_points, strict, names, table_bytes, max_depth, max_patches, directions
        )
    elseif splitter == :oriented
        problem isa StaticMap || throw(ArgumentError("Polygonal splitting supports static maps"))
        partition === nothing || throw(ArgumentError("Construct a fresh polygon map or refine a PiecewisePolygonMap"))
        length(lower) == 2 || throw(ArgumentError("Polygonal ADS currently supports two dimensions"))
        return polygon_construct(
            problem.f, polygon_box(lower, upper); order, atol, rtol, guard_order,
            estimator, splitter, check_points, strict, names, table_bytes,
            max_depth, max_patches, directions
        )
    end
    directions === nothing || throw(ArgumentError("Supply directions with splitter=:oriented"))
    1 <= order <= 65535 && 0 <= guard_order <= 65535 - order || throw(ArgumentError("Invalid expansion or guard order"))
    (estimator isa GuardedTail ? guard_order > 0 : guard_order == 0) ||
        throw(ArgumentError("guard_order must be positive for GuardedTail and zero for other estimators"))
    splitter in (:tail, :width) || throw(ArgumentError("splitter must be :tail or :width"))
    0 <= max_depth <= typemax(Int) && 1 <= max_patches <= typemax(Int) || throw(ArgumentError("Invalid splitting limits"))
    isfinite(rtol) && rtol >= 0 || throw(ArgumentError("rtol must be finite and nonnegative"))
    T = mapreduce(x -> typeof(float(x)), promote_type, Iterators.flatten((lower, upper)))
    T <: AbstractFloat && isconcretetype(T) || throw(ArgumentError("Bounds must convert to a concrete floating-point type"))
    lo, hi = T.(lower), T.(upper)
    all(i -> isfinite(lo[i]) && isfinite(hi[i]) && lo[i] <= hi[i], eachindex(lo)) ||
        throw(ArgumentError("Bounds must be finite and ordered"))
    if partition !== nothing
        length(partition.patches) <= max_patches || throw(ArgumentError("max_patches is smaller than the existing partition"))
        maximum(p -> p.depth, partition.patches) <= max_depth || throw(ArgumentError("max_depth is smaller than the existing partition"))
    end
    options = (;
        order = Int(order), atol, rtol, estimator, splitter, check_points,
        max_depth = Int(max_depth), max_patches = Int(max_patches), strict,
    )
    return with_algebra(Int(order + guard_order), length(lo); names, table_bytes) do ctx
        x = [variable(i, T) for i in eachindex(lo)]
        return ads_build(problem, lo, hi, x, ctx, options, partition)
    end
end

function interval_ads end
interval_ads(args...; kwargs...) = throw(ArgumentError("Load IntervalArithmetic to select IntervalBound"))

"Select certified ADS on an explicit interval box with `estimator=IntervalBound()`."
function adaptive_map(f, box::AbstractVector{<:Real}; estimator::ADSEstimator = IntervalBound(), kwargs...)
    estimator isa IntervalBound || throw(ArgumentError("Supply lower and upper bounds for ordinary box ADS"))
    return interval_ads(f, box; kwargs...)
end

function ads_tolerances(::Type{C}, atol, rtol, n) where {C}
    atol isa Union{Real, AbstractVector{<:Real}} || throw(ArgumentError("atol must be a real scalar or vector"))
    a = atol isa Real ? fill(C(atol), n) : C.(collect(atol))
    length(a) == n || throw(DimensionMismatch("Provide one absolute tolerance per output"))
    all(t -> isfinite(t) && t >= 0, a) || throw(ArgumentError("atol must be finite and nonnegative"))
    all(t -> t > 0 || C(rtol) > 0, a) || throw(ArgumentError("Each output needs a positive tolerance"))
    isfinite(C(rtol)) || throw(ArgumentError("rtol is not representable in the coefficient type"))
    return a
end

function ads_outputs(value)
    result = value isa Real ? [value] : value
    result isa AbstractVector && !isempty(result) && all(x -> x isa Real, result) ||
        throw(ArgumentError("The map must return a real scalar or a nonempty vector of reals"))
    Base.require_one_based_indexing(result)
    return result, value isa Real
end

function ads_check_context(ctx)
    CURRENT_ALGEBRA[] === ctx && ctx.cutoff == ctx.basis.order && iszero(ctx.epsilon) &&
        ctx.big_epsilon === nothing || throw(ArgumentError("The callback changed the algebra configuration"))
    return nothing
end

function ads_geometry(lo, hi)
    center = lo ./ 2 .+ hi ./ 2
    # Enclose both endpoints even when the midpoint rounds toward one of them.
    radius = max.(center - lo, hi - center)
    return (; lo, hi, center, radius)
end

function ads_analyze(value, box, options, ctx)
    values, scalar = ads_outputs(value)
    ads_check_context(ctx)
    C = foldl((R, p) -> promote_type(R, p isa TaylorPolynomial ? coefficient_type(p) : typeof(p)), values; init = eltype(box.lo))
    polynomials = TaylorPolynomial{C}.(values)
    all(p -> valid(p) === ctx && all(isfinite, @view(p.coeffs[1:p.len])), polynomials) ||
        throw(ArgumentError("The map returned invalid or nonfinite Taylor coefficients"))
    compiled = CompiledMap([trim(p, 0, options.order) for p in polynomials])
    errors = zeros(C, length(values))
    contributions = zeros(C, length(values), length(box.lo))
    for (j, p) in enumerate(polynomials)
        errors[j] = ads_error(p, options.order, options.estimator)
        ads_contributions!(contributions, p, j, options.order, options.estimator)
    end
    scale = abs.(constant_term.(polynomials))
    return (; box..., compiled, errors, contributions, scale, scalar)
end

function ads_probes(box)
    (; lo, hi, center, radius) = box
    n = length(lo)
    probes = Tuple{Vector{eltype(lo)}, Vector{eltype(lo)}, Int}[]
    # Axis endpoints, two opposite corners and two alternating corners. Linear
    # growth in dimension avoids an exponential number of callback evaluations.
    for probe in 1:(2n + 4)
        axis = probe <= 2n ? (probe + 1) ÷ 2 : 0
        point, normalized = similar(lo), similar(lo)
        for i in 1:n
            sign = axis > 0 ? (i == axis ? (isodd(probe) ? -1 : 1) : 0) :
                (isodd(probe) ? -1 : 1) * (probe <= 2n + 2 || isodd(i) ? 1 : -1)
            point[i] = sign == 0 ? center[i] : sign < 0 ? lo[i] : hi[i]
            normalized[i] = iszero(radius[i]) ? 0 : (point[i] - center[i]) / radius[i]
        end
        any(p -> p[1] == point, probes) || push!(probes, (point, normalized, axis))
    end
    return probes
end

function ads_check_points!(p, probes, values)
    (; errors, contributions, scale, compiled, scalar) = p
    predicted = similar(errors)
    work = zeros(eltype(errors), degree(compiled) + 1)
    for ((_, normalized, axis), value) in zip(probes, values)
        actual, isscalar = ads_outputs(value)
        length(actual) == length(errors) && isscalar == scalar || throw(DimensionMismatch("The map changed output shape"))
        all(y -> !(y isa TaylorPolynomial) && isfinite(y), actual) ||
            throw(ArgumentError("Numeric evaluation must return finite numeric values"))
        evaluate!(predicted, compiled, normalized, work)
        for j in eachindex(errors)
            discrepancy = abs(actual[j] - predicted[j])
            errors[j] = max(errors[j], discrepancy)
            scale[j] = max(scale[j], abs(actual[j]))
            axis > 0 && (contributions[j, axis] = max(contributions[j, axis], discrepancy))
        end
    end
    return nothing
end

function ads_assess(p, options)
    C = eltype(p.errors)
    tolerance = ads_tolerances(C, options.atol, options.rtol, length(p.errors)) .+ C(options.rtol) .* p.scale
    accepted = all(j -> isfinite(p.errors[j]) && p.errors[j] <= tolerance[j], eachindex(tolerance))
    return (; p..., tolerance, accepted)
end

function ads_candidate(problem::StaticMap, box, x, ctx, options, allow_split)
    p = ads_analyze(problem.f(box.center .+ box.radius .* x), box, options, ctx)
    if options.check_points
        probes = ads_probes(box)
        ads_check_points!(p, probes, (problem.f(point) for (point, _, _) in probes))
    end
    ads_check_context(ctx)
    return ads_assess(p, options)
end

ads_ratio(error, tolerance) = iszero(tolerance) ? (iszero(error) ? zero(error) : oftype(error, Inf)) : error / tolerance

function ads_split_axis(candidate, initial_radius, splitter)
    axis = 0; best = -one(eltype(candidate.errors)); widest = zero(eltype(initial_radius))
    for i in eachindex(initial_radius)
        candidate.lo[i] < candidate.center[i] < candidate.hi[i] || continue
        score = splitter == :width ? zero(best) : maximum(j -> ads_ratio(candidate.contributions[j, i], candidate.tolerance[j]), eachindex(candidate.tolerance))
        width = candidate.radius[i] / initial_radius[i]
        if score > best || (score == best && width > widest)
            axis, best, widest = i, score, width
        end
    end
    return axis
end

function ads_queue(lo::Vector{T}, hi, partition) where {T}
    partition === nothing && return [DomainNode(T)], [(lo, hi, 0, 1)]
    nodes = copy(partition.nodes)
    pending = Tuple{Vector{T}, Vector{T}, Int, Int}[]
    for (i, node) in enumerate(nodes)
        node.patch == 0 && continue
        p = partition.patches[node.patch]
        push!(pending, (copy(p.lower), copy(p.upper), p.depth, i))
    end
    reverse!(pending)
    return nodes, pending
end

function ads_can_split(box, depth, leaves, options)
    return depth < options.max_depth && leaves < options.max_patches &&
        any(i -> box.lo[i] < box.center[i] < box.hi[i], eachindex(box.lo))
end

function ads_build(problem, lo::Vector{T}, hi, x, ctx, options, partition) where {T}
    nodes, pending = ads_queue(lo, hi, partition)
    leaves = length(pending)
    lower, upper, depth, _ = last(pending)
    box = ads_geometry(lower, upper)
    first_patch = ads_candidate(problem, box, x, ctx, options, ads_can_split(box, depth, leaves, options))
    return ads_build!(problem, lo, hi, x, ctx, options, first_patch, nodes, pending)
end

function ads_build!(problem, lo::Vector{T}, hi, x, ctx, options, first_patch, nodes, pending) where {T}
    C = eltype(first_patch.errors)
    patches = TaylorPatch{T, C}[]
    leaves = length(pending)
    initial_radius = ads_geometry(lo, hi).radius
    first_node = last(pending)[4]
    while !isempty(pending)
        lower, upper, depth, node = pop!(pending)
        box = ads_geometry(lower, upper)
        p = node == first_node ? first_patch :
            ads_candidate(problem, box, x, ctx, options, ads_can_split(box, depth, leaves, options))
        p.scalar == first_patch.scalar && length(p.errors) == length(first_patch.errors) || throw(DimensionMismatch("The map changed output shape"))
        eltype(p.errors) === C || throw(ArgumentError("The map changed coefficient type between patches"))
        axis = p.accepted ? 0 : ads_split_axis(p, initial_radius, options.splitter)
        status = p.accepted ? :converged : depth >= options.max_depth ? :max_depth :
            leaves >= options.max_patches ? :max_patches : axis == 0 ? :roundoff : :split
        if status == :split
            midpoint = p.center[axis]
            left, right = length(nodes) + 1, length(nodes) + 2
            nodes[node] = DomainNode(axis, midpoint, left, right, 0)
            push!(nodes, DomainNode(T), DomainNode(T))
            left_upper, right_lower = copy(upper), copy(lower)
            left_upper[axis] = right_lower[axis] = midpoint
            push!(pending, (right_lower, upper, depth + 1, right), (lower, left_upper, depth + 1, left))
            leaves += 1
        else
            options.strict && !p.accepted && throw(ErrorException("Automatic domain splitting reached $status at depth $depth; increase the limit or use strict=false to inspect unresolved patches"))
            push!(patches, TaylorPatch(lower, upper, p.center, p.radius, p.compiled, p.errors, depth, status))
            nodes[node] = DomainNode(0, zero(T), 0, 0, length(patches))
        end
    end
    max_degree = maximum(p -> degree(p.map), patches)
    return PiecewiseTaylorMap(lo, hi, patches, nodes, max_degree, first_patch.scalar, all(p -> p.status == :converged, patches))
end

function ads_patch(map::PiecewiseTaylorMap, point)
    Base.require_one_based_indexing(point)
    length(point) == nvariables(map) || throw(DimensionMismatch("Provide one coordinate per domain dimension"))
    all(i -> !(point[i] isa TaylorPolynomial) && map.lower[i] <= point[i] <= map.upper[i], eachindex(point)) ||
        throw(DomainError(point, "Point lies outside the map domain"))
    node = map.nodes[1]
    while node.patch == 0
        node = map.nodes[point[node.axis] <= node.midpoint ? node.left : node.right]
    end
    return map.patches[node.patch]
end

"""
    evaluate!(result, map::PiecewiseTaylorMap, point, normalized, work)

Evaluate a piecewise map using reusable buffers. `normalized` has
`nvariables(map)` entries; `work` has at least `degree(map) + 1` entries.
All three buffers have the same real element type and must not alias each
other or `point`. This method always returns the output vector, even for a
scalar-valued map. The three-argument form allocates the two work buffers.
"""
function evaluate!(
        out::AbstractVector{R}, map::PiecewiseTaylorMap, point::AbstractVector{<:Real},
        normalized::AbstractVector{R}, work::AbstractVector{R}
    ) where {R <: Real}
    Base.require_one_based_indexing(out, normalized, work)
    length(normalized) == nvariables(map) || throw(DimensionMismatch("Invalid coordinate workspace"))
    buffers = (out, point, normalized, work)
    for i in 1:3, j in (i + 1):4
        Base.mightalias(buffers[i], buffers[j]) && throw(ArgumentError("Evaluation buffers must not alias"))
    end
    patch = ads_patch(map, point)
    for i in eachindex(normalized)
        normalized[i] = iszero(patch.radius[i]) ? zero(R) : (point[i] - patch.center[i]) / patch.radius[i]
    end
    return evaluate!(out, patch.map, normalized, work)
end
function evaluate!(out::AbstractVector{R}, map::PiecewiseTaylorMap, point::AbstractVector{<:Real}) where {R <: Real}
    return evaluate!(out, map, point, zeros(R, nvariables(map)), zeros(R, degree(map) + 1))
end
function evaluate(map::PiecewiseTaylorMap{T, C}, point::AbstractVector{<:Real}) where {T, C}
    R = evaluation_type(promote_type(T, C), point)
    out = evaluate!(zeros(R, noutputs(map)), map, point)
    return map.scalar ? only(out) : out
end
(map::PiecewiseTaylorMap)(point::AbstractVector{<:Real}) = evaluate(map, point)
