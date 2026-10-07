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
function monomial_bound(coefficients, exponents, powers, ::Type{I}) where {I <: Interval}
    result = zero(I)
    @inbounds for k in eachindex(coefficients)
        c = asinterval(I, coefficients[k])
        checked_interval(c)
        DA.coefficient_iszero(c) && continue
        term = c
        for v in axes(exponents, 1)
            term *= powers[v, exponents[v, k] + 1]
        end
        result += term
    end
    return result
end
function polynomial_bound(p, powers, ::Type{I}) where {I <: Interval}
    b = DA.valid(p).basis
    return monomial_bound(@view(p.coeffs[1:p.len]), b.exponents, powers, I)
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


Base.sincos(a::TaylorModel) = (sin(a), cos(a))

function compile_models(models::AbstractVector{TaylorModel{I}}) where {I <: Interval}
    validated_rounding()
    reference = first(models)
    b = DA.model_valid(reference).basis
    foreach(m -> DA.model_valid(m), models)
    # Each patch owns one coordinate family and exponent matrix. Outputs share
    # these read-only snapshots, but own their coefficients and remainders.
    coordinates = deepcopy(reference._coordinates)
    exponents = Matrix{Int}(b.exponents[:, 1:maximum(m -> m._polynomial.len, models)])
    return Tuple(
        DA.CompiledTaylorModel(
            deepcopy(m._polynomial.coeffs[1:m._polynomial.len]), exponents,
            deepcopy(m._remainder), coordinates, m._order, DA.degree(m)
        ) for m in models
    )
end
DA.compile(a::TaylorModel{I}) where {I <: Interval} = only(compile_models([a]))

function snapshot_bound(a::DA.CompiledTaylorModel{I}, normalized) where {I <: Interval}
    powers = power_cache(normalized, DA.degree(a), I)
    return monomial_bound(a._coefficients, a._exponents, powers, I) + a._remainder
end
function snapshot_bounds(models::Tuple{Vararg{DA.CompiledTaylorModel{I}}}, box) where {I <: Interval}
    normalized = normalized_box(first(models), box)
    powers = power_cache(normalized, maximum(DA.degree, models), I)
    return [monomial_bound(m._coefficients, m._exponents, powers, I) + m._remainder for m in models]
end
function DA.enclose(a::DA.CompiledTaylorModel)
    validated_rounding()
    return snapshot_bound(a, a._coordinates.normalized)
end
DA.enclose(a::DA.CompiledTaylorModel, box) = snapshot_bound(a, normalized_box(a, box))
DA.evaluate(a::DA.CompiledTaylorModel, point) = snapshot_bound(a, normalized_box(a, point; point = true))

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
