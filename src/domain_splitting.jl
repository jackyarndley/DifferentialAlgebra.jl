"""
    TaylorPatch

One box in a [`PiecewiseTaylorMap`](@ref). `lower` and `upper` are physical
coordinate bounds; `map` is a [`CompiledMap`](@ref) on normalized coordinates
`(point - center) ./ radius`. A fixed coordinate has normalized value zero.
`error_estimate` contains one absolute truncation-error estimate per output.
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
        guard_order = 2, max_depth = 20, max_patches = 1024,
        check_points = true, strict = true, names = nothing,
        table_bytes = 32 * 1024^2
    )

Approximate `f` on a box by automatic domain splitting (ADS). `f` accepts a
vector of coordinates and returns a real scalar or a nonempty vector of reals.
It must support both ordinary numbers and [`TaylorPolynomial`](@ref)s and be
deterministic. Bounds are finite, ordered vectors; equal bounds fix a coordinate.
Their promoted floating-point type determines the input coefficient type.

Each box is expanded to `order + guard_order`, then truncated to `order`.
The error estimate combines the L1 norm of discarded coefficients, an
exponential extrapolation to the next degree, and (by default) discrepancies
at axis endpoints and selected corners. The split direction maximizes the
estimated reduction in discarded coefficients, with point checks and relative
box widths as fallbacks. Both children recompute `f` on their own domains.

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
"""
function adaptive_map(
        f, lower::AbstractVector{<:Real}, upper::AbstractVector{<:Real};
        order::Integer = 5, atol = 1.0e-8, rtol::Real = 0, guard_order::Integer = 2,
        max_depth::Integer = 20, max_patches::Integer = 1024,
        check_points::Bool = true, strict::Bool = true, names = nothing,
        table_bytes::Integer = 32 * 1024^2
    )
    Base.require_one_based_indexing(lower, upper)
    length(lower) == length(upper) || throw(DimensionMismatch("Domain bounds differ in length"))
    isempty(lower) && throw(ArgumentError("The domain must have at least one coordinate"))
    1 <= order <= 65534 && 1 <= guard_order <= 65535 - order || throw(ArgumentError("Invalid expansion or guard order"))
    0 <= max_depth <= typemax(Int) && 1 <= max_patches <= typemax(Int) || throw(ArgumentError("Invalid splitting limits"))
    isfinite(rtol) && rtol >= 0 || throw(ArgumentError("rtol must be finite and nonnegative"))
    T = mapreduce(x -> typeof(float(x)), promote_type, Iterators.flatten((lower, upper)))
    T <: AbstractFloat && isconcretetype(T) || throw(ArgumentError("Bounds must convert to a concrete floating-point type"))
    lo, hi = T.(lower), T.(upper)
    all(i -> isfinite(lo[i]) && isfinite(hi[i]) && lo[i] <= hi[i], eachindex(lo)) ||
        throw(ArgumentError("Bounds must be finite and ordered"))
    return with_algebra(Int(order + guard_order), length(lo); names, table_bytes) do ctx
        x = [variable(i, T) for i in eachindex(lo)]
        first_patch = ads_candidate(f, lo, hi, x, Int(order), check_points, ctx)
        C = eltype(first_patch.errors)
        tolerances = ads_tolerances(C, atol, rtol, length(first_patch.errors))
        return ads_build(
            f, lo, hi, x, Int(order), check_points, ctx, first_patch,
            tolerances, C(rtol), Int(max_depth), Int(max_patches), strict
        )
    end
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

# Fit log(norm) against degree, ignoring zero orders. Only extrapolate when
# discarded terms exist: an exactly represented polynomial needs no splitting.
function ads_next_norm(norms)
    R = eltype(norms)
    count = 0; sx = sy = sxx = sxy = zero(R)
    for d in 1:(length(norms) - 1)
        c = norms[d + 1]
        iszero(c) && continue
        count += 1
        y = log(c)
        sx += d; sy += y; sxx += R(d)^2; sxy += d * y
    end
    count < 2 && return zero(R)
    slope = (count * sxy - sx * sy) / (count * sxx - sx * sx)
    return exp((sy - slope * sx) / count + slope * length(norms))
end

function ads_check_context(ctx)
    CURRENT_ALGEBRA[] === ctx && ctx.cutoff == ctx.basis.order && iszero(ctx.epsilon) &&
        ctx.big_epsilon === nothing || throw(ArgumentError("The callback changed the algebra configuration"))
    return nothing
end

function ads_candidate(f, lo::Vector{T}, hi, x, order, check_points, ctx) where {T}
    center = lo ./ 2 .+ hi ./ 2
    # Enclose both endpoints even when the midpoint rounds toward one of them.
    radius = max.(center - lo, hi - center)
    values, scalar = ads_outputs(f(center .+ radius .* x))
    ads_check_context(ctx)
    C = foldl((R, p) -> promote_type(R, p isa TaylorPolynomial ? coefficient_type(p) : typeof(p)), values; init = T)
    polynomials = TaylorPolynomial{C}.(values)
    all(p -> valid(p) === ctx && all(isfinite, @view(p.coeffs[1:p.len])), polynomials) ||
        throw(ArgumentError("The map returned invalid or nonfinite Taylor coefficients"))
    compiled = CompiledMap([trim(p, 0, order) for p in polynomials])
    errors = zeros(C, length(values))
    contributions = zeros(C, length(values), length(lo))
    basis = ctx.basis
    for (j, p) in enumerate(polynomials)
        norms = degree_norms(p, 0, 1)
        tail = sum(@view norms[(order + 2):end])
        errors[j] = iszero(tail) ? tail : tail + ads_next_norm(norms)
        for k in (basis.ends[order + 1] + 1):p.len
            c = abs(p.coeffs[k])
            iszero(c) && continue
            for i in eachindex(lo)
                # Halving coordinate i reduces this monomial by 2^(-exponent).
                contributions[j, i] += c * (1 - exp2(-C(basis.exponents[i, k])))
            end
        end
    end
    scale = abs.(constant_term.(polynomials))
    check_points && ads_check_points!(errors, contributions, scale, f, compiled, lo, hi, center, radius, scalar)
    ads_check_context(ctx)
    return (; lo, hi, center, radius, compiled, errors, contributions, scale, scalar)
end

function ads_check_points!(errors, contributions, scale, f, compiled, lo, hi, center, radius, scalar)
    C = eltype(errors); n = length(lo)
    normalized = zeros(eltype(lo), n)
    point = copy(center)
    predicted = similar(errors)
    work = zeros(C, degree(compiled) + 1)
    # Axis endpoints, two opposite corners and two alternating corners. Linear
    # growth in dimension avoids an exponential number of callback evaluations.
    for probe in 1:(2n + 4)
        axis = probe <= 2n ? (probe + 1) ÷ 2 : 0
        for i in 1:n
            sign = axis > 0 ? (i == axis ? (isodd(probe) ? -1 : 1) : 0) :
                (isodd(probe) ? -1 : 1) * (probe <= 2n + 2 || isodd(i) ? 1 : -1)
            point[i] = sign == 0 ? center[i] : sign < 0 ? lo[i] : hi[i]
            normalized[i] = iszero(radius[i]) ? 0 : (point[i] - center[i]) / radius[i]
        end
        actual, isscalar = ads_outputs(f(point))
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

ads_ratio(error, tolerance) = iszero(tolerance) ? (iszero(error) ? zero(error) : oftype(error, Inf)) : error / tolerance

function ads_split_axis(candidate, tolerance, initial_radius)
    axis = 0; best = -one(eltype(candidate.errors)); widest = zero(eltype(initial_radius))
    for i in eachindex(initial_radius)
        candidate.lo[i] < candidate.center[i] < candidate.hi[i] || continue
        score = maximum(j -> ads_ratio(candidate.contributions[j, i], tolerance[j]), eachindex(tolerance))
        width = candidate.radius[i] / initial_radius[i]
        if score > best || (score == best && width > widest)
            axis, best, widest = i, score, width
        end
    end
    return axis
end

function ads_build(
        f, lo::Vector{T}, hi, x, order, check_points, ctx, first_patch,
        atol::Vector{C}, rtol, max_depth, max_patches, strict
    ) where {T, C}
    patches = TaylorPatch{T, C}[]
    nodes = [DomainNode(T)]
    pending = [(lo, hi, 0, 1)]
    leaves = 1
    initial_radius = first_patch.radius
    while !isempty(pending)
        lower, upper, depth, node = pop!(pending)
        p = node == 1 ? first_patch : ads_candidate(f, lower, upper, x, order, check_points, ctx)
        p.scalar == first_patch.scalar && length(p.errors) == length(atol) || throw(DimensionMismatch("The map changed output shape"))
        eltype(p.errors) === C || throw(ArgumentError("The map changed coefficient type between patches"))
        tolerance = atol .+ rtol .* p.scale
        accepted = all(j -> isfinite(p.errors[j]) && p.errors[j] <= tolerance[j], eachindex(atol))
        axis = accepted ? 0 : ads_split_axis(p, tolerance, initial_radius)
        status = accepted ? :converged : depth >= max_depth ? :max_depth :
            leaves >= max_patches ? :max_patches : axis == 0 ? :roundoff : :split
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
            strict && !accepted && throw(ErrorException("Automatic domain splitting reached $status at depth $depth; increase the limit or use strict=false to inspect unresolved patches"))
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
