# Local validated fitting and outward-rounded box geometry. The ADS lifecycle
# itself is in DA.ads_build, shared with ordinary boxes, polygons and flows.
struct ValidatedFit{I <: Interval}
    models::Vector{TaylorModel{I}}
end

# Bound the original function minus its stored midpoint-coefficient polynomial.
# All retained widths and the entire absolute remainder enter acceptance.
function fit_error(a::TaylorModel{I}, powers) where {I <: Interval}
    p = DA.model_polycopy(a._polynomial)
    for k in 1:p.len
        p.coeffs[k] -= IA.interval(IA.numtype(I), IA.mid(p.coeffs[k]))
    end
    return polynomial_bound(p, powers, I) + a._remainder
end

function validated_outputs(value, reference)
    scalar = value isa Union{Real, TaylorModel}
    values = scalar ? [value] : value
    values isa AbstractVector && !isempty(values) || throw(ArgumentError("Return a model/real scalar or nonempty output vector"))
    Base.require_one_based_indexing(values)
    models = map(values) do x
        if x isa TaylorModel
            compatible_models(x, reference)
            return x
        elseif x isa Real && !(x isa TaylorPolynomial)
            return TaylorModel(x, reference)
        end
        throw(ArgumentError("Certified ADS requires outputs computed from its Taylor-model inputs"))
    end
    return models, scalar
end

# Boxes, projected polygons and overlapping continuity fits use the same
# callback normalization, model validation and full fitting-error calculation.
function validated_fit(f, inputs, reference, ctx, order)
    models, scalar = validated_outputs(f(inputs), reference)
    DA.ads_check_context(ctx, order)
    I = DA.coefficient_type(reference)
    powers = power_cache(reference._coordinates.normalized, maximum(DA.degree, models), I)
    errors = Tuple(fit_error(m, powers) for m in models)
    return (; payload = ValidatedFit(models), errors, scalar)
end

function validated_tolerances(::Type{I}, atol, n) where {I <: Interval}
    values = atol isa Real ? fill(atol, n) : atol
    values isa AbstractVector && length(values) == n || throw(DimensionMismatch("Provide one absolute tolerance per output"))
    Base.require_one_based_indexing(values)
    return map(values) do t
        v = checked_interval(asinterval(I, t); guaranteed = true)
        IA.inf(v) > 0 || throw(ArgumentError("Absolute tolerances must be finite and strictly positive"))
        IA.inf(v) # Conservative comparison, including for nonbinary rationals.
    end
end

# A scalar remainder has no exact directional attribution. Use the last two
# nonlinear retained degrees as a heuristic, plus coefficient uncertainty at
# lower degrees. Exact affine variation contributes nothing. If these contain
# no directional information the driver falls back to relative widths.
function interval_contribution(m, axis)
    p = m._polynomial
    b = p.algebra.basis
    T = IA.numtype(typeof(m._remainder))
    value = zero(T)
    tail_start = max(2, m._order - 1)
    for k in 2:p.len
        power = b.exponents[axis, k]
        iszero(power) && continue
        c = p.coeffs[k]
        magnitude = b.degrees[k] >= tail_start ? IA.sup(abs(c)) :
            IA.sup(abs(c - IA.interval(T, IA.mid(c))))
        value += magnitude * (1 - exp2(-power))
    end
    return value
end
fit_error_size(error::Interval) = IA.sup(abs(error))

function validated_candidate(f, inputs, reference, geometry, ctx, options)
    fit = validated_fit(f, inputs, reference, ctx, options.working_order)
    models = fit.payload.models
    I = DA.coefficient_type(reference)
    tolerance = validated_tolerances(I, options.atol, length(models))
    contributions = zeros(IA.numtype(I), length(models), length(reference._coordinates.box))
    for j in eachindex(models)
        fit_error_size(fit.errors[j]) <= tolerance[j] && continue
        for i in axes(contributions, 2)
            contributions[j, i] = interval_contribution(models[j], i)
        end
    end
    accepted = all(j -> fit_error_size(fit.errors[j]) <= tolerance[j], eachindex(tolerance))
    return DA.ADSCandidate(geometry, fit.payload, fit.errors, contributions, tolerance, fit.scalar, accepted)
end

DA.ads_root_node(::DA.BoxADS{<:Tuple{Vararg{Interval{T}}}}) where {T} = DA.DomainNode(T)
DA.ads_dimensions(box::Tuple{Vararg{Interval}}) = length(box)
DA.ads_splittable(box::Tuple{Vararg{Interval}}, i) = IA.inf(box[i]) < IA.mid(box[i]) < IA.sup(box[i])
function DA.ads_relative_width(g::DA.BoxADS, box::Tuple{Vararg{Interval}}, i)
    original = g.root
    T = IA.numtype(box[i])
    return IA.sup(
        (IA.interval(T, IA.sup(box[i])) - IA.interval(T, IA.inf(box[i]))) /
            (IA.interval(T, IA.sup(original[i])) - IA.interval(T, IA.inf(original[i])))
    )
end
DA.ads_patch_domain(::DA.BoxADS{<:Tuple{Vararg{Interval}}}, p::DA.TaylorModelPatch) = Tuple(DA.domain(p))
function DA.ads_patch_domain(::DA.BoxADS{<:Tuple{Vararg{Interval{T}}}}, p::DA.TaylorPatch) where {T}
    return Tuple(IA.interval(T, a, b) for (a, b) in zip(p.lower, p.upper))
end
function DA.ads_split(::DA.BoxADS, box::Tuple{Vararg{Interval{T}}}, p, axis, left, right) where {T}
    midpoint = IA.mid(box[axis])
    left_box = ntuple(i -> i == axis ? IA.intersect_interval(box[i], IA.interval(T, IA.inf(box[i]), midpoint); dec = :auto) : box[i], length(box))
    right_box = ntuple(i -> i == axis ? IA.intersect_interval(box[i], IA.interval(T, midpoint, IA.sup(box[i])); dec = :auto) : box[i], length(box))
    return left_box, right_box, DA.DomainNode(axis, midpoint, left, right, 0)
end
function DA.ads_candidate(problem::DA.StaticMap, box::Tuple{Vararg{Interval}}, x, ctx, options, allow_split)
    inputs = coordinate_models(box, options.order)
    return validated_candidate(problem.f, inputs, first(inputs), box, ctx, options)
end
DA.ads_payload_patch_type(::ValidatedFit{I}, box) where {I} = DA.TaylorModelPatch{I}
function DA.ads_local_patch(fit::ValidatedFit, box, errors, depth, status, order)
    return DA.TaylorModelPatch(compile_models(fit.models), deepcopy(errors), depth, status)
end
function DA.ads_result(g::DA.BoxADS, patches::Vector{DA.TaylorModelPatch{I}}, nodes, scalar, options) where {I}
    return DA.PiecewiseTaylorModel(g.root, Tuple(patches), nodes, options.order, scalar, options.estimator, all(p -> p.status == :converged, patches))
end

function DA.interval_ads(
        f, box; order::Integer = 5, atol = 1.0e-8, rtol::Real = 0,
        guard_order::Integer = 0, check_points::Bool = false, splitter::Symbol = :tail,
        directions = nothing, max_depth::Integer = 20, max_patches::Integer = 1024,
        strict::Bool = true, names = nothing, table_bytes::Integer = 32 * 1024^2,
        partition = nothing
    )
    validated_rounding()
    box isa Union{Tuple, AbstractVector} && !isempty(box) || throw(ArgumentError("Supply a nonempty interval box"))
    interval_box(box, length(box); guaranteed = true)
    T = mapreduce(IA.numtype, promote_type, box)
    T in (Float64, BigFloat) || throw(ArgumentError("Certified ADS supports Float64 or BigFloat endpoints"))
    options = DA.ads_options(;
        order, atol, rtol, estimator = DA.IntervalBound(), splitter,
        guard_order, check_points, max_depth, max_patches, strict, splitters = (:tail, :width, :oriented)
    )
    original = Tuple(deepcopy(asinterval(Interval{T}, x)) for x in box)
    if splitter == :oriented
        partition === nothing || throw(ArgumentError("Refine a box partition with splitter=:tail or :width"))
        length(box) == 2 || throw(ArgumentError("Polygonal ADS currently supports two dimensions"))
        polygon = DA.polygon_box([IA.inf(x) for x in original], [IA.sup(x) for x in original])
        return DA.polygon_construct(f, polygon; order, atol, estimator = DA.IntervalBound(), splitter, directions, max_depth, max_patches, strict, names, table_bytes)
    end
    directions === nothing || throw(ArgumentError("Supply directions with splitter=:oriented"))
    return DA.with_algebra(options.working_order, length(box); names, table_bytes) do ctx
        DA.ads_build(DA.StaticMap(f), DA.BoxADS(original, nothing), ctx, options, partition)
    end
end
function DA.interval_ads(f, lower::AbstractVector, upper::AbstractVector; kwargs...)
    Base.require_one_based_indexing(lower, upper)
    length(lower) == length(upper) && !isempty(lower) || throw(DimensionMismatch("Supply equal nonempty coordinate bounds"))
    all(x -> x isa SupportedScalar, Iterators.flatten((lower, upper))) || throw(ArgumentError("Supply scalar endpoint bounds"))
    T = mapreduce(x -> typeof(float(x)), promote_type, Iterators.flatten((lower, upper)))
    T in (Float64, BigFloat) || throw(ArgumentError("IntervalBound supports Float64 or BigFloat endpoints"))
    all(i -> isfinite(lower[i]) && isfinite(upper[i]) && lower[i] <= upper[i], eachindex(lower)) ||
        throw(ArgumentError("Bounds must be finite and ordered"))
    box = [checked_interval(IA.interval(T, a, b); guaranteed = true) for (a, b) in zip(lower, upper)]
    return DA.interval_ads(f, box; kwargs...)
end

function piecewise_query(a::DA.PiecewiseTaylorModel{I}, box; point = false) where {I <: Interval}
    validated_rounding()
    box isa Union{Tuple, AbstractVector} || throw(ArgumentError("Supply one physical coordinate per dimension"))
    box isa AbstractVector && Base.require_one_based_indexing(box)
    length(box) == DA.nvariables(a) || throw(DimensionMismatch("Provide exactly one physical coordinate per dimension"))
    !point && interval_box(box, length(box); guaranteed = true)
    query = Tuple(checked_interval(asinterval(I, x); guaranteed = true) for x in box)
    point && any(x -> x isa Interval && !IA.isthin(x), box) && throw(ArgumentError("Use enclose for interval subboxes"))
    all(i -> IA.issubset_interval(query[i], a._domain[i]), eachindex(query)) || throw(DomainError(box, "Query lies outside the ADS validity domain"))
    result = nothing
    DA.ads_visit_intersections(a.nodes, (IA.inf.(query), IA.sup.(query))) do index
        patch = a.patches[index]
        intersection = Tuple(IA.intersect_interval(query[i], first(patch.models)._coordinates.box[i]; dec = :auto) for i in eachindex(query))
        any(IA.isempty_interval, intersection) && return nothing
        bounds = snapshot_bounds(patch.models, intersection)
        if result === nothing
            result = bounds
        else
            for j in eachindex(bounds)
                result[j] = IA.hull(result[j], bounds[j]; dec = :auto)
            end
        end
        return nothing
    end
    result === nothing && throw(ErrorException("ADS partition did not cover the query"))
    return a._scalar ? only(result) : result
end
DA.enclose(a::DA.PiecewiseTaylorModel) = piecewise_query(a, a._domain)
DA.enclose(a::DA.PiecewiseTaylorModel, box) = piecewise_query(a, box)
DA.evaluate(a::DA.PiecewiseTaylorModel, point) = piecewise_query(a, point; point = true)
