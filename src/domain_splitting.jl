"""
    TaylorPatch

One local Taylor map with its box subdomain in a [`PiecewiseTaylorMap`](@ref).
`lower` and `upper` are physical coordinate bounds; `map` is a [`CompiledMap`](@ref)
on normalized coordinates
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
    _order::Int
end
domain(p::TaylorPatch) = (lower = deepcopy(p.lower), upper = deepcopy(p.upper))
nvariables(p::TaylorPatch) = length(p.lower)
noutputs(p::TaylorPatch) = noutputs(p.map)
degree(p::TaylorPatch) = degree(p.map)
max_order(p::TaylorPatch) = p._order

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
Its input subdomains form a partition, with one local Taylor map per patch.
Call `map(point)` with physical coordinates to select and evaluate a patch.
`patches` holds its [`TaylorPatch`](@ref)s; `converged` is true only if every
patch meets the requested error estimate. `lower` and `upper` bound the domain.
Numeric evaluation remains valid after the global algebra is reinitialized.
Construction and refinement own the domain endpoints, including BigFloat data;
`domain(map)` returns independent bounds.
Shared faces belong to the lower child of a split. Extrapolation is rejected.
`max_order(map)` is the requested retained order; `degree(map)` is the largest
actual stored polynomial degree. Refinement defaults to the saved order and estimator.
"""
struct PiecewiseTaylorMap{T <: AbstractFloat, C <: Real, E <: ADSEstimator}
    lower::Vector{T}
    upper::Vector{T}
    patches::Vector{TaylorPatch{T, C}}
    nodes::Vector{DomainNode{T}}
    max_degree::Int
    _order::Int
    scalar::Bool
    _estimator::E
    converged::Bool
end
nvariables(map::PiecewiseTaylorMap) = length(map.lower)
noutputs(map::PiecewiseTaylorMap) = noutputs(first(map.patches).map)
degree(map::PiecewiseTaylorMap) = map.max_degree
max_order(map::PiecewiseTaylorMap) = map._order
domain(map::PiecewiseTaylorMap) = (lower = deepcopy(map.lower), upper = deepcopy(map.upper))
Base.copy(map::PiecewiseTaylorMap) = PiecewiseTaylorMap(
    deepcopy(map.lower), deepcopy(map.upper), deepcopy(map.patches), deepcopy(map.nodes),
    map.max_degree, map._order, map.scalar, map._estimator, map.converged
)

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
`max_depth` counts splits along a path; `max_patches` bounds the number of patches.

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
    adaptive_map(f, previous::PiecewiseTaylorMap; order = max_order(previous), estimator = previous._estimator, kwargs...)

Refine an existing partition using the same options as `adaptive_map(f, lower,
upper)`. Each patch is reevaluated with `f`; the input map is unchanged. Existing
boundaries are retained, and `max_depth` and `max_patches` apply to the entire
tree. Supply the original function, not `previous` as a surrogate: subdividing
an already truncated polynomial cannot recover missing information.
The retained order defaults to the originally requested order, even if an output
has accidentally lower degree; the estimator defaults to the saved selection.
"""
function adaptive_map(
        f, previous::PiecewiseTaylorMap; order::Integer = max_order(previous),
        estimator::ADSEstimator = previous._estimator, kwargs...
    )
    return ads_construct(StaticMap(f), previous.lower, previous.upper, previous; order, estimator, kwargs...)
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
        return interval_ads(
            problem.f, lower, upper; partition, order, atol, rtol, guard_order, splitter,
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
    options = ads_options(;
        order, atol, rtol, estimator, splitter, guard_order,
        max_depth, max_patches, check_points, strict
    )
    T = mapreduce(x -> typeof(float(x)), promote_type, Iterators.flatten((lower, upper)))
    T <: AbstractFloat && isconcretetype(T) || throw(ArgumentError("Bounds must convert to a concrete floating-point type"))
    # Conversion alone can reuse BigFloat endpoint storage, including on refinement.
    lo, hi = deepcopy.(T.(lower)), deepcopy.(T.(upper))
    all(i -> isfinite(lo[i]) && isfinite(hi[i]) && lo[i] <= hi[i], eachindex(lo)) ||
        throw(ArgumentError("Bounds must be finite and ordered"))
    return with_algebra(Int(order + guard_order), length(lo); names, table_bytes) do ctx
        x = [variable(i, T) for i in eachindex(lo)]
        return ads_build(problem, BoxADS(ads_geometry(lo, hi), x), ctx, options, partition)
    end
end

function interval_ads end
interval_ads(args...; kwargs...) = throw(ArgumentError("Load IntervalArithmetic to select IntervalBound"))

"Select certified ADS on an explicit interval box with `estimator=IntervalBound()`."
function adaptive_map(f, box::Union{AbstractVector{<:Real}, Tuple{Vararg{Real}}}; estimator::ADSEstimator = IntervalBound(), kwargs...)
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

function ads_check_context(ctx, order = ctx.basis.order)
    CURRENT_ALGEBRA[] === ctx && ctx.cutoff == order && iszero(ctx.epsilon) &&
        ctx.big_epsilon === nothing || throw(ArgumentError("The callback changed the algebra configuration"))
    return nothing
end

function ads_geometry(lo, hi)
    # Preserve fixed subnormal coordinates: halving each endpoint can underflow.
    center = map((a, b) -> a == b ? copy(a) : a / 2 + b / 2, lo, hi)
    # Enclose both endpoints even when the midpoint rounds toward one of them.
    radius = max.(center - lo, hi - center)
    return (; lo, hi, center, radius)
end

function ads_analyze(value, box, options, ctx)
    values, scalar = ads_outputs(value)
    ads_check_context(ctx, options.working_order)
    C = foldl((R, p) -> promote_type(R, p isa TaylorPolynomial ? coefficient_type(p) : typeof(p)), values; init = eltype(box.lo))
    C <: AbstractFloat || throw(ArgumentError("Ordinary ADS requires floating coefficients; select IntervalBound for interval validation"))
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
    return (; geometry = box, polynomials, compiled, errors, contributions, scale, scalar)
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
    return ADSCandidate(p.geometry, p.compiled, p.errors, p.contributions, tolerance, p.scalar, accepted)
end

function ads_candidate(problem::StaticMap, box, x, ctx, options, allow_split)
    p = ads_analyze(problem.f(box.center .+ box.radius .* x), box, options, ctx)
    if options.check_points
        probes = ads_probes(box)
        ads_check_points!(p, probes, (problem.f(point) for (point, _, _) in probes))
    end
    ads_check_context(ctx, options.working_order)
    return ads_assess(p, options)
end

ads_ratio(error, tolerance) = iszero(tolerance) ? (iszero(error) ? zero(error) : oftype(error, Inf)) : error / tolerance

function ads_patch(map::PiecewiseTaylorMap, point)
    Base.require_one_based_indexing(point)
    length(point) == nvariables(map) || throw(DimensionMismatch("Provide one coordinate per domain dimension"))
    all(i -> !(point[i] isa TaylorPolynomial) && map.lower[i] <= point[i] <= map.upper[i], eachindex(point)) ||
        throw(DomainError(point, "Point lies outside the map domain"))
    return map.patches[ads_point_patch(map.nodes, point)]
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
