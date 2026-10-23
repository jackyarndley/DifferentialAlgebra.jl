# Overlapping fits use the native polynomial engine. Blending does not alter
# coefficients or pretend that independently computed jets match at a face.
struct ContinuityPatch{T, C, V, E}
    core_lower::Vector{PolygonReal}
    core_upper::Vector{PolygonReal}
    lower::Vector{PolygonReal}
    upper::Vector{PolygonReal}
    center::Vector{T}
    radius::Vector{T}
    map::CompiledMap{C}
    certificates::V
    errors::E
    # Rows are support/core lower and core/support upper faces, respectively.
    bounds::Tuple{Matrix{Float64}, Matrix{Float64}}
    tapers::Matrix{Float64}
end

# Unlike a partition tree, sibling supports overlap. Each cut bounds the entire
# corresponding subtree outward; visiting both sides preserves closed faces.
struct ContinuityNode
    axis::Int
    left_upper::Float64
    right_lower::Float64
    left::Int
    right::Int
    patch::Int
end
ads_query_bounds(node::ContinuityNode, query) =
    (query[1][node.axis], query[2][node.axis])
ads_visit_left(node::ContinuityNode, query, lo, hi) = lo <= node.left_upper
ads_visit_right(node::ContinuityNode, query, lo, hi) = hi >= node.right_lower
function continuity_index!(nodes, patches, indices)
    index = length(nodes) + 1
    push!(nodes, ContinuityNode(0, 0, 0, 0, 0, first(indices)))
    length(indices) == 1 && return index
    widths = [
        maximum(i -> patches[i].center[j], indices) -
            minimum(i -> patches[i].center[j], indices)
            for j in eachindex(first(patches).center)
    ]
    axis = argmax(widths)
    sort!(indices; by = i -> patches[i].center[axis])
    middle = length(indices) ÷ 2
    left, right = indices[1:middle], indices[(middle + 1):end]
    upper = maximum(i -> patches[i].bounds[2][4, axis], left)
    lower = minimum(i -> patches[i].bounds[1][1, axis], right)
    l = continuity_index!(nodes, patches, left)
    r = continuity_index!(nodes, patches, right)
    nodes[index] = ContinuityNode(axis, upper, lower, l, r, 0)
    return index
end
function continuity_index(patches)
    nodes = ContinuityNode[]
    sizehint!(nodes, 2length(patches) - 1)
    continuity_index!(nodes, patches, collect(eachindex(patches)))
    return nodes
end
function continuity_indices(a, bounds)
    indices = Int[]
    ads_visit_intersections(i -> push!(indices, i), a._nodes, bounds)
    # Retain the original summation and hulling order, independent of the tree.
    return sort!(indices)
end

"""
    ContinuousTaylorMap

A partition-of-unity surrogate constructed by `continuous_map`. Calling it or
`evaluate` returns an ordinary smooth approximation, also with ForwardDiff
dual-number inputs when ForwardDiff is loaded. The original ADS map is unchanged.
`continuity` is `:c0`, `:c1` or `:c2`; this is regularity of the surrogate, not
accuracy of its derivatives relative to the original function. The ideal real
surrogate uses exact rational geometry and stored polynomial coefficients;
floating numeric evaluation has ordinary rounding error.

For an interval source, `enclose` retains fresh Taylor-model certificates for
the original function. `error_bounds` bounds that function minus the ideal real
blended surrogate uniformly, including coefficient widths and remainders.
For an ordinary source, only heuristic `error_estimate` is available. No original
ADS tolerance is inherited: overlap fits are reevaluated. Supply `atol` to
`continuous_map` to check their recomputed bounds/estimates explicitly.

Domain/accessor/copy storage is owned. Snapshots survive algebra reinitialization.
Internal data are read-only. A spatial index selects overlapping supports;
conservative floating bounds fall back to exact geometry near faces.
"""
struct ContinuousTaylorMap{D, P, E}
    _domain::D
    _bounds::Tuple{Vector{PolygonReal}, Vector{PolygonReal}}
    _directions::Matrix{PolygonReal}
    _float_directions::Matrix{Float64}
    _direction_bounds::Tuple{Matrix{Float64}, Matrix{Float64}}
    _root_bounds::Tuple{Matrix{Float64}, Matrix{Float64}}
    _domain_filter::Tuple{Matrix{Float64}, Matrix{Float64}}
    _nodes::Vector{ContinuityNode}
    _patches::Vector{P}
    _errors::E
    _order::Int
    _scalar::Bool
    continuity::Symbol
end
domain(a::ContinuousTaylorMap) = deepcopy(a._domain)
nvariables(a::ContinuousTaylorMap) = size(a._directions, 2)
noutputs(a::ContinuousTaylorMap) = noutputs(first(a._patches).map)
max_order(a::ContinuousTaylorMap) = a._order
degree(::ContinuousTaylorMap) = throw(ArgumentError("A blended map is not a polynomial; max_order gives its retained local order"))
(a::ContinuousTaylorMap)(point) = evaluate(a, point)
function Base.copy(a::ContinuousTaylorMap)
    return ContinuousTaylorMap(
        deepcopy(a._domain), deepcopy(a._bounds), deepcopy(a._directions),
        copy(a._float_directions), deepcopy(a._direction_bounds),
        deepcopy(a._root_bounds), deepcopy(a._domain_filter), copy(a._nodes),
        deepcopy(a._patches),
        deepcopy(a._errors), a._order, a._scalar, a.continuity
    )
end
Base.deepcopy_internal(a::ContinuousTaylorMap, copies::IdDict) = get!(() -> copy(a), copies, a)
function Base.getproperty(a::ContinuousTaylorMap, name::Symbol)
    if name in (:error_bounds, :error_estimate)
        validated = first(getfield(a, :_patches)).certificates !== nothing
        name === :error_bounds && !validated && throw(ArgumentError("Ordinary overlap fits have estimates, not certified bounds"))
        return deepcopy(getfield(a, :_errors))
    end
    return getfield(a, name)
end
Base.propertynames(::ContinuousTaylorMap) = (:continuity, :error_bounds, :error_estimate)
Base.show(io::IO, a::ContinuousTaylorMap) = print(io, "ContinuousTaylorMap(", a.continuity, ", ", length(a._patches), " overlapping fits)")

function continuity_layout(a::PiecewiseTaylorMap)
    lo, hi = polygon_real.(a.lower), polygon_real.(a.upper)
    B = Matrix{PolygonReal}(I, length(lo), length(lo))
    cores = [(polygon_real.(p.lower), polygon_real.(p.upper)) for p in ads_patches(a)]
    return (domain = (lower = copy(a.lower), upper = copy(a.upper)), bounds = (lo, hi), B, A = copy(B), cores, scalar = a.scalar)
end
function continuity_layout(a::PiecewisePolygonMap)
    B, A = polygon_frame(split_directions(a))
    lo, hi = polygon_projection_bounds(a._domain, B)
    cores = [
        (collect(x), collect(y)) for p in ads_patches(a)
            for (x, y) in (polygon_projection_bounds(p._domain, B),)
    ]
    return (domain = domain(a), bounds = (collect(lo), collect(hi)), B, A, cores, scalar = a._scalar)
end
continuity_scalar_type(a::PiecewiseTaylorMap{T}) where {T} = T
continuity_scalar_type(a::PiecewisePolygonMap{T}) where {T} = T

# Exact widening is clipped to the projected root cover. Fixed coordinates
# have weight one and never divide by a zero radius.
function continuity_support(core, root, overlap)
    lo, hi = core
    margin = overlap .* (hi - lo)
    return max.(root[1], lo - margin), min.(root[2], hi + margin)
end
function continuity_check_support(core, lo, hi, T)
    for i in eachindex(lo)
        for gap in (core[1][i] - lo[i], hi[i] - core[2][i])
            gap == 0 && continue
            value = T(gap)
            isfinite(value) && value > 0 || throw(ArgumentError("Overlap taper width is not representable in the source scalar type"))
        end
    end
    return nothing
end
function continuity_candidate(f, lo, hi, A, options, ctx, estimator::ADSEstimator)
    estimator isa IntervalBound && throw(ArgumentError("Load IntervalArithmetic for validated overlap fits"))
    T = options.T
    box = ads_geometry(T.(lo), T.(hi))
    all(i -> isfinite(box.radius[i]) && (lo[i] == hi[i] || box.radius[i] > 0), eachindex(lo)) || throw(ArgumentError("Unrepresentable overlap domain"))
    x = [variable(i, T) for i in eachindex(lo)]
    center, linear = polygon_affine_frame(A, box.center, box.radius, T)
    callback(z) = f(center + linear * ads_normalize.(z, box.center, box.radius))
    candidate = ads_static_fit(StaticMap(callback), box, x, ctx, options)
    C = eltype(candidate.errors)
    C <: AbstractFloat || throw(ArgumentError("Ordinary continuity fits require floating coefficients"))
    return (;
        center = box.center, radius = box.radius, map = candidate.compiled,
        certificates = nothing, errors = candidate.errors, scalar = candidate.scalar,
    )
end
continuity_errors(patches) = map(j -> maximum(p -> p.errors[j], patches), eachindex(first(patches).errors))
function continuity_check_tolerance(errors, atol)
    atol === nothing && return nothing
    tolerance = ads_tolerances(eltype(errors), atol, 0, length(errors))
    all(errors .<= tolerance) || throw(ArgumentError("Overlap fits exceed atol; refine the source map, increase order or reduce overlap"))
    return nothing
end

"""
    continuous_map(
        f, fit; continuity = :c2, overlap = 1 // 4,
        order = max_order(fit), atol = nothing, kwargs...
    )

Reevaluate the original static function on overlapping covers of the box or
polygon ADS subdomains and blend their local polynomials. `overlap` is a positive
fraction (at most one) of each core subdomain's full width added on either side, clipped
to the projected root cover. `:c0`, `:c1` and `:c2` select linear, cubic and
quintic compact tapers. Each weight is one on its original core subdomain's cover and
zero outside its support; normalized weights are nonnegative and sum to one.
The denominator is at least one throughout the original physical domain.

The taper and its first k derivatives join at both endpoints for `:ck`.
Products and normalization therefore give a Ck surrogate across faces and
junctions, including oriented polygons. It generally is not a polynomial and
need not interpolate the original function or match individual patch jets.

For ordinary fits, `estimator` and retained `order` default to the source map's
saved settings, `guard_order` to two for GuardedTail, and `check_points` to true.
An interval source requires IntervalBound and retains its certificates while
using stored midpoint coefficients for the explicitly requested surrogate.
There is no midpoint Taylor-model arithmetic backend. `atol`, if supplied,
checks the newly computed uniform errors (heuristic for ordinary fits) and
throws on failure; the source tolerance alone says nothing about enlarged fits.
The callback must be valid on every overlapping cover, which can extend beyond
a physical polygon, and must not change algebra settings. No time/derivative
error certification is inferred. Construction restores the caller's algebra.
"""
function continuous_map(
        f, fit::Union{PiecewiseTaylorMap, PiecewiseTaylorModel, PiecewisePolygonMap};
        continuity::Symbol = :c2, overlap::Real = 1 // 4,
        order::Integer = max_order(fit), atol = nothing,
        estimator::ADSEstimator = fit._estimator,
        guard_order::Integer = estimator isa GuardedTail ? 2 : 0,
        check_points::Bool = !(estimator isa IntervalBound), names = nothing,
        table_bytes::Integer = 32 * 1024^2
    )
    continuity in (:c0, :c1, :c2) || throw(ArgumentError("continuity must be :c0, :c1 or :c2"))
    fraction = polygon_real(overlap)
    0 < fraction <= 1 || throw(ArgumentError("overlap must lie in (0,1]"))
    ads_check_order(order, guard_order, estimator)
    validated = fit isa PiecewiseTaylorModel || fit isa PiecewisePolygonMap && fit._estimator isa IntervalBound
    validated == (estimator isa IntervalBound) || throw(ArgumentError("Use IntervalBound only with an interval source fit"))
    validated && check_points && throw(ArgumentError("IntervalBound does not use sampled acceptance"))
    atol === nothing || ads_check_atol(atol, 0, estimator)
    layout = continuity_layout(fit)
    T = continuity_scalar_type(fit)
    options = (;
        order = Int(order), working_order = Int(order + guard_order),
        atol, rtol = 0, estimator, T, check_points, names, table_bytes,
        interval_domain = ads_interval_domain(fit),
    )
    return with_algebra(Int(order + guard_order), length(layout.bounds[1]); names, table_bytes) do ctx
        patches = nothing
        for core in layout.cores
            lo, hi = continuity_support(core, layout.bounds, fraction)
            continuity_check_support(core, lo, hi, T)
            p = continuity_candidate(f, lo, hi, layout.A, options, ctx, estimator)
            p.scalar == layout.scalar && noutputs(p.map) == noutputs(fit) || throw(DimensionMismatch("Overlap callback changed output shape"))
            faces = permutedims(hcat(lo, core[1], core[2], hi))
            tapers = permutedims(hcat(core[1], core[2], core[1] - lo, hi - core[2]))
            patch = ContinuityPatch(
                deepcopy(core[1]), deepcopy(core[2]), lo, hi, p.center, p.radius,
                p.map, p.certificates, p.errors, ads_float_bounds(faces), Float64.(tapers)
            )
            patches === nothing && (patches = typeof(patch)[])
            patch isa eltype(patches) || throw(ArgumentError("Overlap callback changed coefficient type"))
            push!(patches, patch)
        end
        errors = continuity_errors(patches)
        continuity_check_tolerance(errors, atol)
        return ContinuousTaylorMap(
            deepcopy(layout.domain), deepcopy(layout.bounds), deepcopy(layout.B),
            Float64.(layout.B), ads_float_bounds(layout.B),
            ads_float_bounds(hcat(layout.bounds...)),
            layout.domain isa ConvexPolygon ? polygon_filter(layout.domain) :
                (zeros(0, 0), zeros(0, 0)),
            continuity_index(patches),
            patches, errors, Int(order), layout.scalar, continuity
        )
    end
end

function continuity_point(a, point)
    point isa Union{Tuple, AbstractVector} && length(point) == nvariables(a) || throw(DimensionMismatch("Provide exactly one physical coordinate per variable"))
    point isa AbstractVector && Base.require_one_based_indexing(point)
    x = polygon_real.(ads_primal_coordinate.(point))
    inside = a._domain isa ConvexPolygon ? polygon_contains(a._domain, x) : all(i -> a._bounds[1][i] <= x[i] <= a._bounds[2][i], eachindex(x))
    inside || throw(DomainError(point, "Point lies outside the continuity map domain"))
    return [sum(a._directions[i, j] * x[j] for j in eachindex(x)) for i in eachindex(x)]
end
function continuity_taper(t, order)
    u = one(t) - t
    return order == 0 ? u : order == 1 ? u^2 * (1 + 2t) : u^3 * (1 + 3t + 6t^2)
end
function continuity_query(a, point)
    point isa Union{Tuple, AbstractVector} && length(point) == nvariables(a) ||
        throw(DimensionMismatch("Provide exactly one physical coordinate per variable"))
    point isa AbstractVector && Base.require_one_based_indexing(point)
    x = ads_primal_coordinate.(point)
    all(isfinite, x) || throw(ArgumentError("Physical coordinates must be finite"))
    bounds = ads_float_bounds(x)
    lower, upper = a._root_bounds
    inside = if a._domain isa ConvexPolygon
        polygon_contains_filtered(a._domain, x, bounds, a._domain_filter)
    elseif all(
            i -> bounds[1][i] >= upper[i, 1] && bounds[2][i] <= lower[i, 2],
            eachindex(x)
        )
        true
    elseif any(
            i -> bounds[2][i] < lower[i, 1] || bounds[1][i] > upper[i, 2],
            eachindex(x)
        )
        false
    else
        all(
            i -> a._bounds[1][i] <= polygon_real(x[i]) <= a._bounds[2][i],
            eachindex(x)
        )
    end
    inside || throw(DomainError(point, "Point lies outside the continuity map domain"))
    return x, ads_project_bounds(a._direction_bounds, bounds)
end
function continuity_side(a, p, point, bounds, i)
    lo, hi = bounds[1][i], bounds[2][i]
    lower, upper = p.bounds
    lo >= upper[2, i] && hi <= lower[3, i] && return 0
    (hi <= lower[1, i] || lo >= upper[4, i]) && return 2
    lo > upper[1, i] && hi < lower[2, i] && return -1
    lo > upper[3, i] && hi < lower[4, i] && return 1
    # Only an inconclusive filter reaches exact projection. This includes core
    # faces, closed support endpoints, fixed axes and sub-Float64 geometry.
    x = sum(a._directions[i, j] * polygon_real(point[j]) for j in eachindex(point))
    p.core_lower[i] <= x <= p.core_upper[i] && return 0
    p.lower[i] < x < p.upper[i] || return 2
    return x < p.core_lower[i] ? -1 : 1
end
function continuity_weight(a, p, z, point, bounds, order)
    T = eltype(z)
    cached = coefficient_primal(first(z)) isa Float64
    weight = one(first(z))
    for i in eachindex(z)
        side = continuity_side(a, p, point, bounds, i)
        side == 0 && continue
        side == 2 && return zero(weight)
        if side == -1
            face = cached ? T(p.tapers[1, i]) : T(p.core_lower[i])
            gap = cached ? T(p.tapers[3, i]) : T(p.core_lower[i] - p.lower[i])
            t = (face - z[i]) / gap
        else
            face = cached ? T(p.tapers[2, i]) : T(p.core_upper[i])
            gap = cached ? T(p.tapers[4, i]) : T(p.upper[i] - p.core_upper[i])
            t = (z[i] - face) / gap
        end
        weight *= continuity_taper(clamp(t, zero(t), one(t)), order)
    end
    return weight
end
function continuity_coordinates(a, point)
    primal, bounds = continuity_query(a, point)
    T = evaluation_type(eltype(first(a._patches).center), point)
    coordinates = T[point...]
    B = coefficient_primal(first(coordinates)) isa Float64 ?
        a._float_directions : T.(a._directions)
    z = B * coordinates
    indices = continuity_indices(a, bounds)
    order = a.continuity == :c0 ? 0 : a.continuity == :c1 ? 1 : 2
    weights = [
        continuity_weight(a, a._patches[i], z, primal, bounds, order) for i in indices
    ]
    return z, indices, weights
end
"""
    blend_weights(map::ContinuousTaylorMap, physical_point)

Return normalized compact-support weights. They are nonnegative and sum to
one in exact arithmetic, with ordinary rounding in numeric/dual evaluation.
"""
function blend_weights(a::ContinuousTaylorMap, point)
    _, indices, active = continuity_coordinates(a, point)
    weights = zeros(eltype(active), length(a._patches))
    weights[indices] .= active ./ sum(active)
    return weights
end
function evaluate(a::ContinuousTaylorMap, point)
    z, indices, weights = continuity_coordinates(a, point)
    R = promote_type(eltype(first(a._patches).map.coefficients), eltype(z))
    result = zeros(R, noutputs(a))
    for (index, w) in zip(indices, weights)
        p = a._patches[index]
        iszero(coefficient_primal(w)) && continue
        normalized = ads_normalize.(z, p.center, p.radius)
        result .+= w .* evaluate(p.map, normalized)
    end
    result ./= sum(weights)
    return a._scalar ? only(result) : result
end
enclose(a::ContinuousTaylorMap, args...) = throw(ArgumentError("Certified continuity enclosures require an IntervalBound source and IntervalArithmetic"))
