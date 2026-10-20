# One construction lifecycle for static boxes/polygons and ordinary flows.
# Geometry supplies domains and splits; local fitting supplies payloads and errors.
struct ADSCandidate{G, P, E, C, T}
    geometry::G
    payload::P
    errors::E
    contributions::C
    tolerance::T
    scalar::Bool
    accepted::Bool
end

# A failed analytic fit has neither a payload nor an output shape. Subdivide it
# geometrically until a valid model establishes the concrete patch storage type.
struct ADSDomainFailure{G, E, T}
    geometry::G
    cause::E
    zero_score::T
end
ads_score_zero(p::ADSCandidate) = zero(eltype(p.contributions))
ads_score_zero(p::ADSDomainFailure) = p.zero_score
ads_direction_score(p::ADSDomainFailure, axis) = p.zero_score

struct BoxADS{D, X}
    root::D
    variables::X
end

function ads_check_order(order, guard_order, estimator)
    1 <= order <= 65535 && 0 <= guard_order <= 65535 - order || throw(ArgumentError("Invalid expansion or guard order"))
    (estimator isa GuardedTail ? guard_order > 0 : guard_order == 0) ||
        throw(ArgumentError("guard_order must be positive for GuardedTail and zero for other estimators"))
    return nothing
end

function ads_check_tolerance(t::Real, rtol, estimator)
    t isa TaylorPolynomial &&
        throw(ArgumentError("Absolute tolerances must be numeric scalars, not polynomials"))
    isfinite(t) && t >= 0 ||
        throw(ArgumentError("atol must be finite and nonnegative"))
    (t > 0 || (!(estimator isa IntervalBound) && rtol > 0)) ||
        throw(ArgumentError("Each output needs a positive tolerance"))
    return nothing
end
function ads_check_atol(atol, rtol, estimator)
    atol isa Real && return ads_check_tolerance(atol, rtol, estimator)
    atol isa AbstractVector{<:Real} ||
        throw(ArgumentError("atol must be a real scalar or vector"))
    Base.require_one_based_indexing(atol)
    isempty(atol) && throw(ArgumentError("Provide an absolute tolerance for each output"))
    foreach(t -> ads_check_tolerance(t, rtol, estimator), atol)
    return nothing
end

function ads_options(;
        order, atol, rtol, estimator, splitter, guard_order, max_depth, max_patches,
        check_points, strict, splitters = (:tail, :width), names = nothing,
        table_bytes::Integer = 32 * 1024^2
    )
    ads_check_order(order, guard_order, estimator)
    splitter in splitters || throw(ArgumentError("Unsupported ADS splitter: $splitter"))
    0 <= max_depth <= typemax(Int) && 1 <= max_patches <= typemax(Int) || throw(ArgumentError("Invalid splitting limits"))
    isfinite(rtol) && rtol >= 0 || throw(ArgumentError("rtol must be finite and nonnegative"))
    estimator isa IntervalBound && (!iszero(rtol) || check_points) &&
        throw(ArgumentError("IntervalBound uses positive absolute bounds, no guard degrees and no sampled acceptance"))
    ads_check_atol(atol, rtol, estimator)
    0 <= table_bytes <= typemax(Int) ||
        throw(ArgumentError("Invalid lookup-table budget"))
    return (;
        order = Int(order), working_order = Int(order + guard_order),
        atol = deepcopy(atol), rtol = deepcopy(rtol),
        estimator, splitter, max_depth = Int(max_depth), max_patches = Int(max_patches),
        check_points, strict, names = deepcopy(names), table_bytes = Int(table_bytes),
    )
end

function ads_saved_options(options)
    return deepcopy(
        (;
            order = options.order, atol = options.atol, rtol = options.rtol,
            splitter = options.splitter,
            guard_order = options.working_order - options.order,
            max_depth = options.max_depth, max_patches = options.max_patches,
            check_points = options.check_points, strict = options.strict,
            names = options.names, table_bytes = options.table_bytes,
        )
    )
end
function ads_refinement_options(previous, estimator; kwargs...)
    saved = deepcopy(previous._options)
    guard_order = 0
    if estimator isa GuardedTail
        guard_order = previous._estimator isa GuardedTail ? saved.guard_order : 2
    end
    check_points = saved.check_points
    if estimator isa IntervalBound
        check_points = false
    elseif previous._estimator isa IntervalBound
        check_points = true
    end
    return (; saved..., guard_order, check_points, kwargs...)
end

function ads_check_partition(partition, options)
    partition === nothing && return nothing
    length(ads_patches(partition)) <= options.max_patches ||
        throw(ArgumentError("max_patches is smaller than the existing partition"))
    maximum(p -> p.depth, ads_patches(partition)) <= options.max_depth ||
        throw(ArgumentError("max_depth is smaller than the existing partition"))
    return nothing
end

function ads_queue(geometry, partition)
    root = geometry.root
    nodes = partition === nothing ? [ads_root_node(geometry)] : deepcopy(partition.nodes)
    pending = Tuple{typeof(root), Int, Int}[]
    if partition === nothing
        push!(pending, (root, 0, 1))
    else
        for (i, node) in enumerate(nodes)
            node.patch == 0 && continue
            p = ads_patches(partition)[node.patch]
            push!(pending, (ads_patch_domain(geometry, p), p.depth, i))
        end
        reverse!(pending)
    end
    return nodes, pending
end

ads_root_node(g::BoxADS) = DomainNode(eltype(g.root.lo))
ads_patch_domain(::BoxADS, p::TaylorPatch) = ads_geometry(deepcopy(p.lower), deepcopy(p.upper))
ads_splittable(box, i) = box.lo[i] < box.center[i] < box.hi[i]
ads_dimensions(box) = length(box.lo)
ads_relative_width(g::BoxADS, box, i) = box.radius[i] / g.root.radius[i]
ads_has_split(box) = any(i -> ads_splittable(box, i), 1:ads_dimensions(box))

function ads_can_split(box, depth, patch_count, options)
    return depth < options.max_depth && patch_count < options.max_patches &&
        ads_has_split(box)
end

function ads_direction_score(p, axis)
    score = zero(eltype(p.contributions))
    for j in eachindex(p.errors)
        score = max(score, ads_ratio(p.contributions[j, axis], p.tolerance[j]))
    end
    return score
end

function ads_split_axis(p, geometry, options)
    axis = 0
    best = -one(ads_score_zero(p))
    widest = zero(best)
    for i in 1:ads_dimensions(p.geometry)
        ads_splittable(p.geometry, i) || continue
        score = options.splitter == :width ? zero(best) : ads_direction_score(p, i)
        width = ads_relative_width(geometry, p.geometry, i)
        if score > best || (score == best && width > widest)
            axis, best, widest = i, score, width
        end
    end
    return axis
end

function ads_split(::BoxADS, box, p, axis, left, right)
    midpoint = box.center[axis]
    left_upper, right_lower = copy(box.hi), copy(box.lo)
    left_upper[axis] = right_lower[axis] = midpoint
    return ads_geometry(box.lo, left_upper), ads_geometry(right_lower, box.hi),
        DomainNode(axis, midpoint, left, right, 0)
end

ads_fit(problem, g::BoxADS, box, ctx, options, allow_split) =
    ads_candidate(problem, box, g.variables, ctx, options, allow_split)

ads_payload_patch_type(::CompiledMap{C}, box) where {C} = TaylorPatch{eltype(box.lo), C}
ads_patch_type(::BoxADS, p) = ads_payload_patch_type(p.payload, p.geometry)
function ads_local_patch(payload::CompiledMap, box, errors, depth, status, order)
    return TaylorPatch(
        deepcopy(box.lo), deepcopy(box.hi), deepcopy(box.center), deepcopy(box.radius),
        payload, copy(errors), depth, status, order
    )
end
ads_finalize(::BoxADS, box, p, depth, status, order) =
    ads_local_patch(p.payload, p.geometry, p.errors, depth, status, order)

function ads_result(g::BoxADS, patches::Vector{TaylorPatch{T, C}}, nodes, scalar, options) where {T, C}
    max_degree = maximum(p -> degree(p.map), patches)
    return PiecewiseTaylorMap(
        deepcopy(g.root.lo), deepcopy(g.root.hi), patches, nodes, max_degree, options.order,
        scalar,
        options.estimator, ads_saved_options(options),
        all(p -> p.status == :converged, patches)
    )
end

function ads_candidate_status(p, geometry, depth, patch_count, options)
    accepted = p isa ADSCandidate && p.accepted
    axis = accepted ? 0 : ads_split_axis(p, geometry, options)
    status = accepted ? :converged : depth >= options.max_depth ? :max_depth :
        patch_count >= options.max_patches ? :max_patches : axis == 0 ? :roundoff : :split
    return axis, status
end

function ads_enqueue_split!(geometry, box, p, depth, node, axis, nodes, pending)
    left, right = length(nodes) + 1, length(nodes) + 2
    left_box, right_box, branch = ads_split(geometry, box, p, axis, left, right)
    nodes[node] = branch
    push!(nodes, ads_leaf(branch, 0), ads_leaf(branch, 0))
    push!(pending, (right_box, depth + 1, right), (left_box, depth + 1, left))
    return nothing
end

function ads_retry_domain!(
        geometry, box, p, depth, node, nodes, pending,
        patch_count, options
    )
    axis, status = ads_candidate_status(p, geometry, depth, patch_count, options)
    if status != :split
        cause = sprint(showerror, p.cause)
        throw(
            ErrorException(
                "Cannot construct a valid interval model: $status at depth $depth. " *
                    "Check the function domain or refine the input box. Cause: $cause"
            )
        )
    end
    return ads_enqueue_split!(geometry, box, p, depth, node, axis, nodes, pending)
end

function ads_build(problem, geometry, ctx, options, partition = nothing)
    ads_check_partition(partition, options)
    nodes, pending = ads_queue(geometry, partition)
    while true
        box, depth, node = last(pending)
        first_candidate = ads_fit(
            problem, geometry, box, ctx, options,
            ads_can_split(box, depth, length(pending), options)
        )
        if !(first_candidate isa ADSDomainFailure)
            return ads_build!(problem, geometry, ctx, options, first_candidate, nodes, pending)
        end
        pop!(pending)
        ads_retry_domain!(
            geometry, box, first_candidate, depth, node, nodes, pending,
            length(pending) + 1, options
        )
    end
    return
end

# The first candidate is the function barrier: no runtime patch counts in types,
# and a concrete patch vector for every estimator/geometry combination.
function ads_build!(problem, geometry, ctx, options, first_candidate, nodes, pending)
    patches = ads_patch_type(geometry, first_candidate)[]
    patch_count = length(pending)
    first_node = last(pending)[3]
    while !isempty(pending)
        box, depth, node = pop!(pending)
        p = node == first_node ? first_candidate :
            ads_fit(problem, geometry, box, ctx, options, ads_can_split(box, depth, patch_count, options))
        if p isa ADSDomainFailure
            ads_retry_domain!(
                geometry, box, p, depth, node, nodes, pending,
                patch_count, options
            )
            patch_count += 1
            continue
        end
        p.scalar == first_candidate.scalar && length(p.errors) == length(first_candidate.errors) ||
            throw(DimensionMismatch("The map changed output shape"))
        eltype(p.errors) === eltype(first_candidate.errors) || throw(ArgumentError("The map changed coefficient type between patches"))
        axis, status = ads_candidate_status(p, geometry, depth, patch_count, options)
        if status == :split
            ads_enqueue_split!(geometry, box, p, depth, node, axis, nodes, pending)
            patch_count += 1
        else
            if options.strict && !p.accepted
                hint = status == :roundoff ? "Use wider precision or relax the tolerance" :
                    "Increase the limit or use strict=false to inspect unresolved patches"
                throw(
                    ErrorException(
                        "Automatic domain splitting reached $status at depth $depth. " * hint
                    )
                )
            end
            push!(patches, ads_finalize(geometry, box, p, depth, status, options.order))
            nodes[node] = ads_leaf(nodes[node], length(patches))
        end
    end
    return ads_result(geometry, patches, nodes, first_candidate.scalar, options)
end

ads_leaf(::DomainNode{T}, patch) where {T} = DomainNode(0, zero(T), 0, 0, patch)
ads_node_value(node::DomainNode, point) = point[node.axis]
ads_query_bounds(node::DomainNode, query) = (query[1][node.axis], query[2][node.axis])
function ads_point_patch(nodes, point)
    node = first(nodes)
    while node.patch == 0
        node = nodes[ads_node_value(node, point) <= node.midpoint ? node.left : node.right]
    end
    return node.patch
end

# Closed-set queries visit both sides of a shared face. In particular a subbox
# cannot be routed using its midpoint. Geometry supplies bounds on each cut.
function ads_visit_intersections(f, nodes, query, index = 1)
    node = nodes[index]
    if node.patch != 0
        f(node.patch)
    else
        lo, hi = ads_query_bounds(node, query)
        lo <= node.midpoint && ads_visit_intersections(f, nodes, query, node.left)
        hi >= node.midpoint && ads_visit_intersections(f, nodes, query, node.right)
    end
    return nothing
end
