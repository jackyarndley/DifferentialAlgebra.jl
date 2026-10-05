module DifferentialAlgebraIntervalArithmeticExt

import DifferentialAlgebra as DA
import IntervalArithmetic as IA
import DifferentialAlgebra: TaylorPolynomial, TaylorModel
using IntervalArithmetic: Interval

# IA 1.0.12 also offers non-rigorous rounding=:none. Its configuration record
# is inspected here, without calling configure or changing process settings.
# Keep this version-specific access isolated in the optional extension.
function validated_rounding()
    IA.configuration_options.rounding === :correct ||
        throw(ArgumentError("Validated enclosure and Taylor-model arithmetic require IntervalArithmetic.configure(rounding=:correct)"))
    return nothing
end
model_context(a) = (validated_rounding(); DA.model_arithmetic(a))
compatible_models(a, b) = (validated_rounding(); DA.model_compatible(a, b))

# Skip only a valid, guaranteed exact zero. Even a thin NG/degraded zero must
# participate in arithmetic so that its metadata is not silently erased.
DA.coefficient_iszero(c::Interval) = IA.isthinzero(c) && IA.isguaranteed(c) && IA.decoration(c) == IA.com
DA.coefficient_convert(::Type{Interval{T}}, c::Union{Integer, Rational}) where {T <: IA.NumTypes} = IA.interval(T, c)
DA.coefficient_convert(::Type{Interval{T}}, c::IA.ExactReal) where {T <: IA.NumTypes} = asinterval(Interval{T}, c)
DA.coefficient_operand(::Type{Interval{T}}, c::Real) where {T <: IA.NumTypes} = DA.coefficient_convert(Interval{T}, c)
DA.degree_factor(c::Interval{T}, n::Integer) where {T <: IA.NumTypes} = IA.interval(T, n)
DA.coefficient_muladd(a::Real, b::Real, c::Interval) = DA.coefficient_convert(typeof(c), a) * DA.coefficient_convert(typeof(c), b) + c
DA.invertible_constant(c::Interval) = !IA.isnai(c) && !IA.isempty_interval(c) && !IA.in_interval(0, c)

# IA has symmetric Real promotion/arithmetic rules. Resolve their intersections
# using methods that own a DifferentialAlgebra type, never methods on IA alone.
Base.promote_rule(::Type{TaylorPolynomial{T}}, ::Type{Interval{S}}) where {T <: Real, S <: IA.NumTypes} = TaylorPolynomial{promote_type(T, Interval{S})}
Base.promote_rule(::Type{Interval{S}}, ::Type{TaylorPolynomial{T}}) where {T <: Real, S <: IA.NumTypes} = TaylorPolynomial{promote_type(T, Interval{S})}
Base.promote_rule(::Type{TaylorPolynomial}, ::Type{Interval{S}}) where {S <: IA.NumTypes} = TaylorPolynomial
Base.promote_rule(::Type{Interval{S}}, ::Type{TaylorPolynomial}) where {S <: IA.NumTypes} = TaylorPolynomial
Base.promote_rule(::Type{TaylorPolynomial{T}}, ::Type{IA.ExactReal{S}}) where {T <: Real, S <: Real} = TaylorPolynomial{promote_type(T, S)}
Base.promote_rule(::Type{TaylorPolynomial}, ::Type{IA.ExactReal{S}}) where {S <: Real} = TaylorPolynomial
TaylorPolynomial{T}(x::IA.ExactReal) where {T <: Real} = invoke(TaylorPolynomial{T}, Tuple{Real}, x)
TaylorPolynomial(x::IA.ExactReal) = TaylorPolynomial(x.value)
Base.convert(::Type{TaylorPolynomial{T}}, x::IA.ExactReal) where {T <: Real} = TaylorPolynomial{T}(x)
Base.convert(::Type{TaylorPolynomial}, x::IA.ExactReal) = TaylorPolynomial(x)
for op in (:+, :-, :*, :/)
    @eval begin
        Base.$op(a::TaylorPolynomial{T}, b::Interval) where {T <: Real} = invoke(Base.$op, Tuple{TaylorPolynomial{T}, Real}, a, b)
        Base.$op(a::Interval, b::TaylorPolynomial) = invoke(Base.$op, Tuple{Real, TaylorPolynomial}, a, b)
    end
end

function checked_interval(x::Interval; guaranteed::Bool = false)
    !IA.isnai(x) && !IA.isempty_interval(x) && IA.isbounded(x) && IA.decoration(x) >= IA.def ||
        throw(ArgumentError("Require finite, nonempty intervals without loose-evaluation decorations"))
    guaranteed && !IA.isguaranteed(x) && throw(ArgumentError("Taylor models require guaranteed intervals; an NG input was supplied"))
    return x
end
const SupportedScalar = Union{Integer, Rational, AbstractFloat, AbstractIrrational}
asinterval(::Type{Interval{T}}, x::Interval) where {T <: IA.NumTypes} = IA.interval(T, x) # preserves flags
asinterval(::Type{Interval{T}}, x::SupportedScalar) where {T <: IA.NumTypes} = IA.interval(T, x)
asinterval(::Type{Interval{T}}, x::IA.ExactReal) where {T <: IA.NumTypes} = asinterval(Interval{T}, x.value)
asinterval(::Type{<:Interval}, x) = throw(ArgumentError("Unsupported scalar; supply a real value or decorated Interval explicitly"))

function interval_box(box, n; guaranteed = false)
    box isa Union{Tuple, AbstractVector} || throw(ArgumentError("Supply a tuple or vector of decorated intervals"))
    box isa AbstractVector && Base.require_one_based_indexing(box)
    length(box) == n || throw(DimensionMismatch("Provide exactly one interval per coordinate"))
    all(x -> x isa Interval, box) || throw(ArgumentError("Domains require decorated Interval values; BareInterval is unsupported"))
    foreach(x -> checked_interval(x; guaranteed), box)
    return box
end

# Integer powers, rather than repeated interval products, preserve even powers.
# The matrix caches all powers needed by either a range or a discarded product.
function power_cache(box, order, ::Type{I}) where {I <: Interval}
    powers = Matrix{I}(undef, length(box), order + 1)
    for v in eachindex(box)
        x = asinterval(I, box[v])
        powers[v, 1] = one(I)
        for k in 1:order
            powers[v, k + 1] = IA.pown(x, k)
        end
    end
    return powers
end
function polynomial_bound(p, powers, ::Type{I}) where {I <: Interval}
    b = DA.valid(p).basis
    result = zero(I)
    @inbounds for k in 1:p.len
        c = asinterval(I, p.coeffs[k])
        checked_interval(c)
        DA.coefficient_iszero(c) && continue
        term = c
        for v in 1:b.variables
            term *= powers[v, b.exponents[v, k] + 1]
        end
        result += term
    end
    return result
end
function DA.enclose(p::TaylorPolynomial, box)
    validated_rounding()
    b = DA.valid(p).basis
    interval_box(box, b.variables)
    T = mapreduce(IA.numtype, promote_type, box)
    C = DA.coefficient_type(p)
    C <: Interval && (T = promote_type(T, IA.numtype(C)))
    C <: AbstractFloat && (T = promote_type(T, C))
    T in (Float64, BigFloat) || throw(ArgumentError("Range bounding supports Float64 or BigFloat interval endpoints"))
    I = Interval{T}
    return polynomial_bound(p, power_cache(box, DA.degree(p), I), I)
end

function model_result(p::TaylorPolynomial{I}, r::I, coordinates, order) where {I <: Interval}
    checked_interval(r; guaranteed = true)
    for k in 1:p.len
        checked_interval(p.coeffs[k]; guaranteed = true)
    end
    return TaylorModel(p, r, coordinates, order, Val(:owned))
end
function coordinate_models(box, order)
    validated_rounding()
    box isa Union{Tuple, AbstractVector} && !isempty(box) || throw(ArgumentError("Supply a nonempty interval box"))
    interval_box(box, length(box); guaranteed = true)
    T = mapreduce(IA.numtype, promote_type, box)
    T in (Float64, BigFloat) || throw(ArgumentError("Taylor models support Float64 or BigFloat endpoints"))
    I = Interval{T}
    physical = Tuple(deepcopy(asinterval(I, x)) for x in box)
    centers = Tuple(IA.interval(T, IA.mid(x)) for x in physical)
    radii = Tuple(IA.interval(T, IA.sup(abs(x - c))) for (x, c) in zip(physical, centers))
    foreach(x -> checked_interval(x; guaranteed = true), radii)
    normalized = Tuple(IA.isthinzero(r) ? zero(I) : IA.interval(T, -1, 1) for r in radii)
    coordinates = DA.ModelCoordinates(physical, centers, radii, normalized, Ref(nothing))
    return [model_result(centers[i] + radii[i] * DA.variable(i, I), zero(I), coordinates, Int(order)) for i in eachindex(physical)]
end
function DA.taylor_models(box; order::Integer, names = nothing, table_bytes::Integer = 32 * 1024^2)
    validated_rounding()
    box isa Union{Tuple, AbstractVector} && !isempty(box) || throw(ArgumentError("Supply a nonempty interval box"))
    interval_box(box, length(box); guaranteed = true)
    T = mapreduce(IA.numtype, promote_type, box)
    T in (Float64, BigFloat) || throw(ArgumentError("Taylor models support Float64 or BigFloat endpoints"))
    DA.initialize!(order, length(box); names, table_bytes)
    return coordinate_models(box, order)
end
function TaylorModel(value::Real, reference::TaylorModel{I}) where {I <: Interval}
    model_context(reference)
    value isa TaylorPolynomial && throw(ArgumentError("Supply an explicit absolute remainder when constructing a model from a polynomial"))
    c = checked_interval(asinterval(I, value); guaranteed = true)
    p = TaylorPolynomial{I}(deepcopy(c))
    return model_result(p, zero(I), reference._coordinates, reference._order)
end
function TaylorModel(p::TaylorPolynomial, r::Real, reference::TaylorModel{I}) where {I <: Interval}
    model_context(reference)
    DA.compatible(p, reference._polynomial)
    DA.degree(p) <= reference._order || throw(ArgumentError("Polynomial exceeds the retained model order; no implicit truncation is allowed"))
    coeffs = [deepcopy(checked_interval(asinterval(I, p.coeffs[k]); guaranteed = true)) for k in 1:p.len]
    owned = TaylorPolynomial{I}(coeffs, p.len, p.algebra)
    return model_result(owned, deepcopy(checked_interval(asinterval(I, r); guaranteed = true)), reference._coordinates, reference._order)
end

function normalized_box(a::Union{TaylorModel{I}, DA.CompiledTaylorModel{I}}, box; point = false) where {I <: Interval}
    validated_rounding()
    a isa TaylorModel && DA.model_valid(a)
    c = a._coordinates
    box isa Union{Tuple, AbstractVector} || throw(ArgumentError("Supply one physical coordinate per dimension"))
    box isa AbstractVector && Base.require_one_based_indexing(box)
    length(box) == length(c.box) || throw(DimensionMismatch("Provide exactly one physical coordinate per dimension"))
    !point && interval_box(box, length(c.box); guaranteed = true)
    result = Vector{I}(undef, length(box))
    for i in eachindex(box)
        x = checked_interval(asinterval(I, box[i]); guaranteed = true)
        point && box[i] isa Interval && !IA.isthin(box[i]) && throw(ArgumentError("Point coordinates must be scalar values or thin intervals; use enclose for boxes"))
        IA.issubset_interval(x, c.box[i]) || throw(DomainError(box, "Query lies outside the Taylor-model validity domain"))
        # Both sets enclose the same exact normalized coordinate. The public
        # :auto policy preserves the minimum input decoration and all NG flags;
        # the default set-operation decoration would be trv even for valid sets.
        result[i] = IA.isthinzero(c.radius[i]) ? zero(I) : IA.intersect_interval((x - c.center[i]) / c.radius[i], c.normalized[i]; dec = :auto)
    end
    return result
end
function DA.enclose(a::TaylorModel{I}) where {I <: Interval}
    validated_rounding()
    DA.model_valid(a)
    return polynomial_bound(a._polynomial, power_cache(a._coordinates.normalized, DA.degree(a), I), I) + a._remainder
end
function DA.enclose(a::TaylorModel{I}, physical_subbox) where {I <: Interval}
    normalized = normalized_box(a, physical_subbox)
    return polynomial_bound(a._polynomial, power_cache(normalized, DA.degree(a), I), I) + a._remainder
end
function DA.evaluate(a::TaylorModel{I}, physical_point) where {I <: Interval}
    normalized = normalized_box(a, physical_point; point = true)
    return polynomial_bound(a._polynomial, power_cache(normalized, DA.degree(a), I), I) + a._remainder
end

Base.zero(a::TaylorModel) = TaylorModel(0, a)
Base.one(a::TaylorModel) = TaylorModel(1, a)
function Base.:-(a::TaylorModel)
    model_context(a)
    return model_result(-a._polynomial, -a._remainder, a._coordinates, a._order)
end
for op in (:+, :-)
    @eval function Base.$op(a::TaylorModel{I}, b::TaylorModel{I}) where {I <: Interval}
        compatible_models(a, b)
        return model_result($op(a._polynomial, b._polynomial), $op(a._remainder, b._remainder), a._coordinates, a._order)
    end
end

# Bound every ordered discarded product directly in the original basis. No
# order-2n algebra, exponent-vector allocation, or nearest-rounded sum is used.
function discarded_bound(p, q, order, powers, ::Type{I}) where {I <: Interval}
    basis = p.algebra.basis
    result = zero(I)
    @inbounds for i in 1:p.len
        ac = p.coeffs[i]
        DA.coefficient_iszero(ac) && continue
        for j in 1:q.len
            basis.degrees[i] + basis.degrees[j] > order || continue
            bc = q.coeffs[j]
            DA.coefficient_iszero(bc) && continue
            term = ac * bc
            for v in 1:basis.variables
                term *= powers[v, basis.exponents[v, i] + basis.exponents[v, j] + 1]
            end
            result += term
        end
    end
    return result
end
function Base.:*(a::TaylorModel{I}, b::TaylorModel{I}) where {I <: Interval}
    ctx = compatible_models(a, b)
    p, q = a._polynomial, b._polynomial
    powers = power_cache(a._coordinates.normalized, DA.degree(p) + DA.degree(q), I)
    tail = discarded_bound(p, q, a._order, powers, I)
    rp, rq = polynomial_bound(p, powers, I), polynomial_bound(q, powers, I)
    r = tail + rp * b._remainder + rq * a._remainder + a._remainder * b._remainder
    retained = DA.allocate(ctx, I, ctx.basis.ends[a._order + 1])
    DA.multiply!(retained, p, q, a._order)
    return model_result(retained, r, a._coordinates, a._order)
end
for op in (:+, :-, :*, :/)
    @eval begin
        Base.$op(a::TaylorModel, b::Real) = $op(a, TaylorModel(b, a))
        Base.$op(a::Real, b::TaylorModel) = $op(TaylorModel(a, b), b)
    end
end
Base.:/(a::TaylorModel, b::TaylorModel) = (compatible_models(a, b); a * inv(b))
function Base.:^(a::TaylorModel, n::Integer)
    model_context(a)
    typemin(Int32) < n <= typemax(Int32) || throw(ArgumentError("Exponent out of range"))
    n < 0 && return inv(a)^(-n)
    n == 1 && return copy(a)
    result, factor = one(a), a
    while n > 0
        isodd(n) && (result = result * factor)
        n >>= 1
        n > 0 && (factor = factor * factor)
    end
    return result
end

# f^(k)(X)/k!, with interval constants and no floating factorial approximation.
# These are valid on the whole hull from the expansion point to the input range.
function normalized_derivative(fn, x::I, k::Int) where {I <: Interval}
    T = IA.numtype(I)
    if fn === inv
        return (isodd(k) ? -one(I) : one(I)) / IA.pown(x, k + 1)
    elseif fn === log
        k == 0 && return log(x)
        return (iseven(k) ? -one(I) : one(I)) / (IA.interval(T, k) * IA.pown(x, k))
    elseif fn === sqrt
        factor = one(I)
        for j in 1:k
            factor *= (IA.interval(T, 1 // 2) - IA.interval(T, j - 1)) / IA.interval(T, j)
        end
        return factor * sqrt(x) / IA.pown(x, k)
    end
    factor = one(I)
    for j in 1:k
        factor /= IA.interval(T, j)
    end
    fn === exp && return exp(x) * factor
    values = fn === sin ? (sin(x), cos(x), -sin(x), -cos(x)) : (cos(x), -sin(x), -cos(x), sin(x))
    return values[mod(k, 4) + 1] * factor
end
function elementary_model(fn, a::TaylorModel{I}) where {I <: Interval}
    model_context(a)
    x = checked_interval(DA.enclose(a); guaranteed = true)
    fn === sqrt && IA.isthinzero(x) && return zero(a)
    fn === inv && IA.in_interval(0, x) && throw(DomainError(x, "Reciprocal requires exclusion of zero on the whole model enclosure"))
    fn in (log, sqrt) && IA.inf(x) <= 0 && throw(DomainError(x, "This expansion requires strict positivity on the whole model enclosure"))
    center = IA.interval(IA.numtype(I), IA.mid(x))
    IA.issubset_interval(center, x) || throw(ArgumentError("Could not construct an expansion point inside the input enclosure"))
    h = a - center
    n = a._order
    result = TaylorModel(normalized_derivative(fn, center, n), a)
    for k in (n - 1):-1:0
        result = result * h + normalized_derivative(fn, center, k)
    end
    # Taylor's theorem: the (n+1)th normalized derivative is bounded on x;
    # h^(n+1) here is an interval power of the whole input displacement.
    tail = normalized_derivative(fn, x, n + 1) * IA.pown(x - center, n + 1)
    return model_result(result._polynomial, result._remainder + tail, a._coordinates, n)
end
for fn in (:inv, :exp, :log, :sin, :cos, :sqrt)
    @eval Base.$fn(a::TaylorModel) = elementary_model($fn, a)
end

include("interval_polygon_ads.jl")
Base.sincos(a::TaylorModel) = (sin(a), cos(a))

function DA.compile(a::TaylorModel{I}) where {I <: Interval}
    validated_rounding()
    b = DA.model_valid(a).basis
    return DA.CompiledTaylorModel(
        deepcopy(a._polynomial.coeffs[1:a._polynomial.len]),
        Matrix{Int}(b.exponents[:, 1:a._polynomial.len]), deepcopy(a._remainder),
        deepcopy(a._coordinates), a._order
    )
end

function snapshot_bound(a::DA.CompiledTaylorModel{I}, normalized) where {I <: Interval}
    powers = power_cache(normalized, DA.degree(a), I)
    result = zero(I)
    @inbounds for k in eachindex(a._coefficients)
        term = a._coefficients[k]
        DA.coefficient_iszero(term) && continue
        for v in axes(a._exponents, 1)
            term *= powers[v, a._exponents[v, k] + 1]
        end
        result += term
    end
    return result + a._remainder
end
function DA.enclose(a::DA.CompiledTaylorModel)
    validated_rounding()
    return snapshot_bound(a, a._coordinates.normalized)
end
DA.enclose(a::DA.CompiledTaylorModel, box) = snapshot_bound(a, normalized_box(a, box))
DA.evaluate(a::DA.CompiledTaylorModel, point) = snapshot_bound(a, normalized_box(a, point; point = true))

# This is a split criterion, not a midpoint-coefficient arithmetic backend.
# Bound f minus the exact stored midpoint polynomial, including every retained
# coefficient width. Subtraction and summation are outward-rounded throughout.
function fit_error(a::TaylorModel{I}) where {I <: Interval}
    p = DA.model_polycopy(a._polynomial)
    for k in 1:p.len
        p.coeffs[k] -= IA.interval(IA.numtype(I), IA.mid(p.coeffs[k]))
    end
    return polynomial_bound(p, power_cache(a._coordinates.normalized, DA.degree(p), I), I) + a._remainder
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

function validated_split_axis(box, original, values = nothing, tolerances = nothing, splitter = :width)
    axis = 0
    score = -Inf
    widest = -Inf
    for i in eachindex(box)
        IA.inf(box[i]) < IA.mid(box[i]) < IA.sup(box[i]) || continue
        # Choose a direction using enclosed endpoint spans. The choice affects
        # efficiency only; acceptance always uses the uniform error bound.
        relative = IA.sup(
            (IA.interval(IA.numtype(box[i]), IA.sup(box[i])) - IA.interval(IA.numtype(box[i]), IA.inf(box[i]))) /
                (IA.interval(IA.numtype(original[i]), IA.sup(original[i])) - IA.interval(IA.numtype(original[i]), IA.inf(original[i])))
        )
        sensitivity = values === nothing || splitter == :width ? zero(relative) : maximum(j -> interval_contribution(values[j], i) / tolerances[j], eachindex(values))
        if sensitivity > score || (sensitivity == score && relative > widest)
            axis, score, widest = i, sensitivity, relative
        end
    end
    return axis
end

# Direction scores are heuristic. Acceptance always uses fit_error, never these
# nearest-rounded scalar scores or a sampled discrepancy.
function interval_contribution(m, axis)
    p = m._polynomial
    b = p.algebra.basis
    value = zero(IA.numtype(typeof(m._remainder)))
    for k in 2:p.len
        power = b.exponents[axis, k]
        iszero(power) && continue
        value += IA.sup(abs(p.coeffs[k])) * (1 - exp2(-power))
    end
    return value
end

function DA.interval_ads(
        f, box; order = 5, atol = 1.0e-8, rtol = 0,
        guard_order = 0, check_points = false, splitter = :tail,
        directions = nothing, kwargs...
    )
    iszero(rtol) && iszero(guard_order) && !check_points ||
        throw(ArgumentError("IntervalBound uses positive absolute bounds, no guard degrees and no sampled acceptance"))
    return DA.validated_adaptive_map(f, box; order, atol, splitter, directions, kwargs...)
end
function DA.interval_ads(f, lower::AbstractVector, upper::AbstractVector; kwargs...)
    length(lower) == length(upper) && !isempty(lower) || throw(DimensionMismatch("Supply equal nonempty coordinate bounds"))
    all(x -> x isa SupportedScalar, Iterators.flatten((lower, upper))) || throw(ArgumentError("Supply scalar endpoint bounds"))
    T = mapreduce(x -> typeof(float(x)), promote_type, Iterators.flatten((lower, upper)))
    T in (Float64, BigFloat) || throw(ArgumentError("IntervalBound supports Float64 or BigFloat endpoints"))
    box = [checked_interval(IA.interval(T, a, b); guaranteed = true) for (a, b) in zip(lower, upper)]
    return DA.interval_ads(f, box; kwargs...)
end

function DA.validated_adaptive_map(
        f, box; order::Integer = 3, atol = 1.0e-6,
        max_depth::Integer = 20, max_patches::Integer = 1024, strict::Bool = true,
        names = nothing, table_bytes::Integer = 32 * 1024^2,
        splitter::Symbol = :width, directions = nothing
    )
    validated_rounding()
    box isa Union{Tuple, AbstractVector} && !isempty(box) || throw(ArgumentError("Supply a nonempty interval box"))
    interval_box(box, length(box); guaranteed = true)
    T = mapreduce(IA.numtype, promote_type, box)
    T in (Float64, BigFloat) || throw(ArgumentError("Certified ADS supports Float64 or BigFloat endpoints"))
    1 <= order <= 65535 || throw(ArgumentError("Invalid retained order"))
    0 <= max_depth <= typemax(Int) && 1 <= max_patches <= typemax(Int) || throw(ArgumentError("Invalid splitting limits"))
    I = Interval{T}
    original = Tuple(deepcopy(asinterval(I, x)) for x in box)
    splitter in (:width, :tail, :oriented) || throw(ArgumentError("Unsupported interval ADS splitter"))
    if splitter == :oriented
        length(box) == 2 || throw(ArgumentError("Polygonal ADS currently supports two dimensions"))
        polygon = DA.polygon_box([IA.inf(x) for x in original], [IA.sup(x) for x in original])
        return DA.polygon_construct(f, polygon; order, atol, estimator = DA.IntervalBound(), splitter, directions, max_depth, max_patches, strict, names, table_bytes)
    end
    directions === nothing || throw(ArgumentError("Supply directions with splitter=:oriented"))
    return DA.with_algebra(Int(order), length(box); names, table_bytes) do ctx
        pending = [(original, 0)]
        patches = DA.TaylorModelPatch{I}[]
        leaves = 1
        scalar = nothing
        count = 0
        tolerances = T[]
        while !isempty(pending)
            child, depth = pop!(pending)
            x = coordinate_models(child, order)
            values, isscalar = validated_outputs(f(x), first(x))
            DA.ads_check_context(ctx)
            if scalar === nothing
                scalar, count = isscalar, length(values)
                tolerances = validated_tolerances(I, atol, count)
            end
            scalar == isscalar && count == length(values) || throw(DimensionMismatch("The function changed output shape between patches"))
            errors = Tuple(fit_error(m) for m in values)
            accepted = all(j -> IA.sup(abs(errors[j])) <= tolerances[j], eachindex(errors))
            axis = accepted ? 0 : validated_split_axis(child, original, values, tolerances, splitter)
            status = accepted ? :converged : depth >= max_depth ? :max_depth :
                leaves >= max_patches ? :max_patches : axis == 0 ? :roundoff : :split
            if status == :split
                midpoint = IA.mid(child[axis])
                left = ntuple(i -> i == axis ? IA.intersect_interval(child[i], IA.interval(T, IA.inf(child[i]), midpoint); dec = :auto) : child[i], length(child))
                right = ntuple(i -> i == axis ? IA.intersect_interval(child[i], IA.interval(T, midpoint, IA.sup(child[i])); dec = :auto) : child[i], length(child))
                push!(pending, (right, depth + 1), (left, depth + 1))
                leaves += 1
            else
                strict && !accepted && throw(ErrorException("Certified box ADS reached $status at depth $depth; increase limits or use strict=false"))
                push!(patches, DA.TaylorModelPatch(Tuple(DA.compile(m) for m in values), deepcopy(errors), depth, status))
            end
        end
        return DA.PiecewiseTaylorModel(original, Tuple(patches), Int(order), scalar, all(p -> p.status == :converged, patches))
    end
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
    for patch in a.patches
        intersection = Tuple(IA.intersect_interval(query[i], first(patch.models)._coordinates.box[i]; dec = :auto) for i in eachindex(query))
        any(IA.isempty_interval, intersection) && continue
        bounds = [DA.enclose(m, intersection) for m in patch.models]
        result = result === nothing ? bounds : [IA.hull(result[j], bounds[j]; dec = :auto) for j in eachindex(bounds)]
    end
    result === nothing && throw(ErrorException("ADS partition did not cover the query"))
    return a._scalar ? only(result) : result
end
DA.enclose(a::DA.PiecewiseTaylorModel) = piecewise_query(a, a._domain)
DA.enclose(a::DA.PiecewiseTaylorModel, box) = piecewise_query(a, box)
DA.evaluate(a::DA.PiecewiseTaylorModel, point) = piecewise_query(a, point; point = true)

# Polynomial functions enclose coefficient arithmetic only. They do not acquire
# a Taylor-model remainder. For interval coefficients sqrt avoids interval order
# comparisons and uses an exact enclosed half rather than a floating exponent.
function Base.sqrt(a::TaylorPolynomial{I}) where {I <: Interval}
    x = DA.constant_term(a)
    DA.coefficient_iszero(a) && return zero(a)
    IA.inf(checked_interval(x)) > 0 || throw(DomainError(x, "Polynomial square-root expansion requires a strictly positive constant interval"))
    c0 = sqrt(x)
    return isempty(a.algebra.basis.products) ? DA.power_series(a, IA.interval(IA.numtype(I), 1 // 2), c0) : DA.square_root(a, c0)
end

end
