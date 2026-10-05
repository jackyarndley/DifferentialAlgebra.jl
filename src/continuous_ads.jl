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
Internal data are read-only. Evaluation scans the overlapping supports.
"""
struct ContinuousTaylorMap{D, P, E}
    _domain::D
    _bounds::Tuple{Vector{PolygonReal}, Vector{PolygonReal}}
    _directions::Matrix{PolygonReal}
    _patches::Tuple{Vararg{P}}
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
Base.copy(a::ContinuousTaylorMap) = ContinuousTaylorMap(deepcopy(a._domain), deepcopy(a._bounds), deepcopy(a._directions), deepcopy(a._patches), deepcopy(a._errors), a._order, a._scalar, a.continuity)
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
    cores = [(polygon_real.(p.lower), polygon_real.(p.upper)) for p in a.patches]
    return (domain = (lower = copy(a.lower), upper = copy(a.upper)), bounds = (lo, hi), B, A = copy(B), cores, scalar = a.scalar)
end
function continuity_layout(a::PiecewisePolygonMap)
    B, A = polygon_frame(split_directions(a))
    lo, hi = polygon_projection_bounds(a._domain, B)
    cores = [(collect(x), collect(y)) for p in a.patches for (x, y) in (polygon_projection_bounds(p._domain, B),)]
    return (domain = domain(a), bounds = (collect(lo), collect(hi)), B, A, cores, scalar = a._scalar)
end
continuity_estimator(::PiecewiseTaylorMap) = GuardedTail()
continuity_estimator(::PiecewiseTaylorModel) = IntervalBound()
continuity_estimator(a::PiecewisePolygonMap) = a._estimator
continuity_scalar_type(a::PiecewiseTaylorMap{T}) where {T} = T
continuity_scalar_type(a::PiecewisePolygonMap{T}) where {T} = T
function continuity_layout end
function continuity_scalar_type end

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
    callback(z) = f(T.(A) * z)
    settings = (; order = options.order, atol = 1, rtol = 0, estimator, check_points = options.check_points)
    candidate = ads_candidate(StaticMap(callback), box, x, ctx, settings, false)
    C = eltype(candidate.errors)
    C <: AbstractFloat || throw(ArgumentError("Ordinary continuity fits require floating coefficients"))
    return (; center = candidate.center, radius = candidate.radius, map = candidate.compiled, certificates = nothing, errors = candidate.errors, scalar = candidate.scalar)
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
        order = max(1, degree(fit)), atol = nothing, kwargs...
    )

Reevaluate the original static function on overlapping covers of the box or
polygon ADS leaves and blend their local polynomials. `overlap` is a positive
fraction (at most one) of each leaf's full width added on either side, clipped
to the projected root cover. `:c0`, `:c1` and `:c2` select linear, cubic and
quintic compact tapers. Each weight is one on its original leaf's cover and
zero outside its support; normalized weights are nonnegative and sum to one.
The denominator is at least one throughout the original physical domain.

The taper and its first k derivatives join at both endpoints for `:ck`.
Products and normalization therefore give a Ck surrogate across faces and
junctions, including oriented polygons. It generally is not a polynomial and
need not interpolate the original function or match individual patch jets.

For ordinary fits, `estimator` defaults to GuardedTail (or the polygon source's
estimator), `guard_order` to two for GuardedTail, and `check_points` to true.
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
        order::Integer = max(1, degree(fit)), atol = nothing,
        estimator::ADSEstimator = continuity_estimator(fit),
        guard_order::Integer = estimator isa GuardedTail ? 2 : 0,
        check_points::Bool = !(estimator isa IntervalBound), names = nothing,
        table_bytes::Integer = 32 * 1024^2
    )
    continuity in (:c0, :c1, :c2) || throw(ArgumentError("continuity must be :c0, :c1 or :c2"))
    fraction = polygon_real(overlap)
    0 < fraction <= 1 || throw(ArgumentError("overlap must lie in (0,1]"))
    1 <= order <= 65535 && 0 <= guard_order <= 65535 - order || throw(ArgumentError("Invalid retained/guard order"))
    (estimator isa GuardedTail ? guard_order > 0 : guard_order == 0) || throw(ArgumentError("Invalid guard order for estimator"))
    validated = fit isa PiecewiseTaylorModel || fit isa PiecewisePolygonMap && fit._estimator isa IntervalBound
    validated == (estimator isa IntervalBound) || throw(ArgumentError("Use IntervalBound only with an interval source fit"))
    validated && check_points && throw(ArgumentError("IntervalBound does not use sampled acceptance"))
    layout = continuity_layout(fit)
    T = continuity_scalar_type(fit)
    options = (; order = Int(order), T, check_points)
    return with_algebra(Int(order + guard_order), length(layout.bounds[1]); names, table_bytes) do ctx
        patches = nothing
        for core in layout.cores
            lo, hi = continuity_support(core, layout.bounds, fraction)
            continuity_check_support(core, lo, hi, T)
            p = continuity_candidate(f, lo, hi, layout.A, options, ctx, estimator)
            p.scalar == layout.scalar && noutputs(p.map) == noutputs(fit) || throw(DimensionMismatch("Overlap callback changed output shape"))
            patch = ContinuityPatch(deepcopy(core[1]), deepcopy(core[2]), lo, hi, p.center, p.radius, p.map, p.certificates, p.errors)
            patches === nothing && (patches = typeof(patch)[])
            patch isa eltype(patches) || throw(ArgumentError("Overlap callback changed coefficient type"))
            push!(patches, patch)
        end
        errors = continuity_errors(patches)
        continuity_check_tolerance(errors, atol)
        return ContinuousTaylorMap(deepcopy(layout.domain), deepcopy(layout.bounds), deepcopy(layout.B), Tuple(patches), errors, Int(order), layout.scalar, continuity)
    end
end

continuity_primal(x::Union{Integer, Rational, AbstractFloat}) = x
continuity_primal(x) = throw(ArgumentError("Provide scalar physical coordinates; use enclose for intervals and load ForwardDiff for dual queries"))
function continuity_point(a, point)
    point isa Union{Tuple, AbstractVector} && length(point) == nvariables(a) || throw(DimensionMismatch("Provide exactly one physical coordinate per variable"))
    point isa AbstractVector && Base.require_one_based_indexing(point)
    x = polygon_real.(continuity_primal.(point))
    inside = a._domain isa ConvexPolygon ? polygon_contains(a._domain, x) : all(i -> a._bounds[1][i] <= x[i] <= a._bounds[2][i], eachindex(x))
    inside || throw(DomainError(point, "Point lies outside the continuity map domain"))
    return [sum(a._directions[i, j] * x[j] for j in eachindex(x)) for i in eachindex(x)]
end
function continuity_taper(t, order)
    u = one(t) - t
    return order == 0 ? u : order == 1 ? u^2 * (1 + 2t) : u^3 * (1 + 3t + 6t^2)
end
function continuity_weight(p, z, exact, order)
    T = eltype(p.center)
    weight = one(first(z))
    for i in eachindex(z)
        x = exact[i]
        p.core_lower[i] <= x <= p.core_upper[i] && continue
        p.lower[i] < x < p.upper[i] || return zero(weight)
        t = x < p.core_lower[i] ? (T(p.core_lower[i]) - z[i]) / T(p.core_lower[i] - p.lower[i]) : (z[i] - T(p.core_upper[i])) / T(p.upper[i] - p.core_upper[i])
        weight *= continuity_taper(clamp(t, zero(t), one(t)), order)
    end
    return weight
end
function continuity_coordinates(a, point)
    exact = continuity_point(a, point)
    T = eltype(first(a._patches).center)
    z = T.(a._directions) * collect(point)
    order = a.continuity == :c0 ? 0 : a.continuity == :c1 ? 1 : 2
    weights = [continuity_weight(p, z, exact, order) for p in a._patches]
    return z, weights
end
"""
    blend_weights(map::ContinuousTaylorMap, physical_point)

Return normalized compact-support weights. They are nonnegative and sum to
one in exact arithmetic, with ordinary rounding in numeric/dual evaluation.
"""
function blend_weights(a::ContinuousTaylorMap, point)
    _, weights = continuity_coordinates(a, point)
    return weights ./ sum(weights)
end
function evaluate(a::ContinuousTaylorMap, point)
    z, weights = continuity_coordinates(a, point)
    R = promote_type(eltype(first(a._patches).map.coefficients), eltype(z))
    result = zeros(R, noutputs(a))
    for (p, w) in zip(a._patches, weights)
        iszero(continuity_primal(w)) && continue
        normalized = [iszero(p.radius[i]) ? zero(z[i]) : (z[i] - p.center[i]) / p.radius[i] for i in eachindex(z)]
        result .+= w .* evaluate(p.map, normalized)
    end
    result ./= sum(weights)
    return a._scalar ? only(result) : result
end
enclose(a::ContinuousTaylorMap, args...) = throw(ArgumentError("Certified continuity enclosures require an IntervalBound source and IntervalArithmetic"))
