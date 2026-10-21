function interval_box(box, n; guaranteed = false)
    box isa Union{Tuple, AbstractVector} || throw(ArgumentError("Supply a tuple or vector of decorated intervals"))
    box isa AbstractVector && Base.require_one_based_indexing(box)
    length(box) == n || throw(DimensionMismatch("Provide exactly one interval per coordinate"))
    all(x -> x isa Interval, box) || throw(ArgumentError("Domains require decorated Interval values; BareInterval is unsupported"))
    foreach(x -> checked_interval(x; guaranteed), box)
    return box
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

function model_result(
        p::TaylorPolynomial{I}, r::I, coordinates, order,
        enclosure::Union{Nothing, I} = nothing
    ) where {I <: Interval}
    checked_interval(r; guaranteed = true)
    for k in 1:p.len
        checked_interval(p.coeffs[k]; guaranteed = true)
    end
    if enclosure !== nothing
        !IA.isnai(enclosure) && !IA.isempty_interval(enclosure) &&
            IA.decoration(enclosure) >= IA.def && IA.isguaranteed(enclosure) ||
            throw(ArgumentError("A model range certificate must be valid and guaranteed"))
    end
    return TaylorModel(p, r, deepcopy(enclosure), coordinates, order, Val(:owned))
end
function coordinate_models(box, order)
    validated_rounding()
    box isa Union{Tuple, AbstractVector} && !isempty(box) || throw(ArgumentError("Supply a nonempty interval box"))
    interval_box(box, length(box); guaranteed = true)
    T = mapreduce(IA.numtype, promote_type, box)
    T in (Float64, BigFloat) || throw(ArgumentError("Taylor models support Float64 or BigFloat endpoints"))
    I = Interval{T}
    physical = [deepcopy(asinterval(I, x)) for x in box]
    centers = [
        IA.setdecoration(IA.interval(T, IA.mid(x)), IA.decoration(x))
            for x in physical
    ]
    radii = [
        IA.setdecoration(IA.interval(T, IA.sup(abs(x - c))), IA.decoration(x))
            for (x, c) in zip(physical, centers)
    ]
    foreach(x -> checked_interval(x; guaranteed = true), radii)
    normalized = [
        IA.setdecoration(
            IA.isthinzero(r) ? zero(I) : IA.interval(T, -1, 1),
            IA.decoration(x)
        ) for (x, r) in zip(physical, radii)
    ]
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
function model_scalar(::Type{I}, value::Real) where {I <: Interval}
    value isa TaylorPolynomial && throw(
        ArgumentError(
            "Supply an explicit absolute remainder when constructing a model " *
                "from a polynomial"
        )
    )
    return checked_interval(DA.coefficient_convert(I, value); guaranteed = true)
end
function TaylorModel(value::Real, reference::TaylorModel{I}) where {I <: Interval}
    model_context(reference)
    c = model_scalar(I, value)
    p = TaylorPolynomial{I}(DA.model_valid(reference), deepcopy(c))
    return model_result(p, zero(I), reference._coordinates, reference._order)
end

# A constant has no variable dependence, but its range still carries every
# coordinate's decoration, including unused fixed axes.
function constant_bound(c::I, coordinates) where {I <: Interval}
    result = zero(I)
    for axis in coordinates.normalized
        result *= IA.pown(asinterval(I, axis), 0)
    end
    return result + checked_interval(asinterval(I, c))
end

function model_scalar_sum(
        op, a::TaylorModel{I}, value::Real, reverse = false
    ) where {I <: Interval}
    model_context(a)
    c = model_scalar(I, value)
    z = zero(I)
    sign = reverse && op === (-) ? -one(I) : one(I)
    p = DA.scaled(a._polynomial, sign; filter = false)
    # Match polynomial addition's endpoint precision and own BigFloat zeros.
    map!(coefficient -> coefficient + z, p.coeffs, p.coeffs)
    p.coeffs[1] = op === (-) && !reverse ? p.coeffs[1] - c : p.coeffs[1] + c
    DA.finish!(p, p.len)
    r = reverse ? op(z, a._remainder) : op(a._remainder, z)
    enclosure = if a._enclosure === nothing
        nothing
    else
        range, scalar = model_bound(a), constant_bound(c, a._coordinates)
        reverse ? op(scalar, range) : op(range, scalar)
    end
    return model_result(p, r, a._coordinates, a._order, enclosure)
end
for op in (:+, :-)
    @eval begin
        Base.$op(a::TaylorModel, b::Real) = model_scalar_sum($op, a, b)
        Base.$op(a::Real, b::TaylorModel) = model_scalar_sum($op, b, a, true)
    end
end
function TaylorModel(p::TaylorPolynomial, r::Real, reference::TaylorModel{I}) where {I <: Interval}
    model_context(reference)
    DA.compatible(p, reference._polynomial)
    DA.degree(p) <= reference._order || throw(ArgumentError("Polynomial exceeds the retained model order; no implicit truncation is allowed"))
    coeffs = [deepcopy(checked_interval(asinterval(I, p.coeffs[k]); guaranteed = true)) for k in 1:p.len]
    owned = TaylorPolynomial{I}(coeffs, p.len, p.algebra)
    return model_result(owned, deepcopy(checked_interval(DA.coefficient_convert(I, r); guaranteed = true)), reference._coordinates, reference._order)
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
        result[i] = if IA.isthinzero(c.radius[i])
            IA.setdecoration(zero(I), min(IA.decoration(x), IA.decoration(c.box[i])))
        else
            IA.intersect_interval(
                (x - c.center[i]) / c.radius[i], c.normalized[i];
                dec = :auto
            )
        end
    end
    return result
end
function intersect_range(a, bound)
    return a._enclosure === nothing ? bound :
        IA.intersect_interval(bound, a._enclosure; dec = :auto)
end
polynomial_bound(a::TaylorModel, powers, I) = polynomial_bound(a._polynomial, powers, I)
polynomial_bound(a::DA.CompiledTaylorModel, powers, I) =
    monomial_bound(a._coefficients, a._exponents, powers, I)

# Live models and owned snapshots share range arithmetic; only coefficient
# access differs. Public queries copy the result at their ownership boundary.
model_bound(a, powers, ::Type{I}) where {I <: Interval} =
    intersect_range(a, polynomial_bound(a, powers, I) + a._remainder)
function model_bound(
        a::Union{TaylorModel{I}, DA.CompiledTaylorModel{I}},
        normalized = a._coordinates.normalized
    ) where {I <: Interval}
    powers = power_cache(normalized, DA.degree(a), I)
    return model_bound(a, powers, I)
end
function DA.enclose(a::Union{TaylorModel, DA.CompiledTaylorModel})
    validated_rounding()
    a isa TaylorModel && DA.model_valid(a)
    return deepcopy(model_bound(a))
end
function DA.enclose(a::Union{TaylorModel, DA.CompiledTaylorModel}, physical_subbox)
    normalized = normalized_box(a, physical_subbox)
    return deepcopy(model_bound(a, normalized))
end
function DA.evaluate(a::Union{TaylorModel, DA.CompiledTaylorModel}, physical_point)
    normalized = normalized_box(a, physical_point; point = true)
    return deepcopy(model_bound(a, normalized))
end

Base.zero(a::TaylorModel) = TaylorModel(0, a)
Base.one(a::TaylorModel) = TaylorModel(1, a)
function Base.:-(a::TaylorModel)
    model_context(a)
    enclosure = a._enclosure === nothing ? nothing : -a._enclosure
    return model_result(
        -a._polynomial, -a._remainder, a._coordinates, a._order, enclosure
    )
end
for op in (:+, :-)
    @eval function Base.$op(a::TaylorModel{I}, b::TaylorModel{I}) where {I <: Interval}
        compatible_models(a, b)
        enclosure = a._enclosure === nothing && b._enclosure === nothing ? nothing :
            $op(model_bound(a), model_bound(b))
        return model_result(
            $op(a._polynomial, b._polynomial), $op(a._remainder, b._remainder),
            a._coordinates, a._order, enclosure
        )
    end
end

function Base.:*(a::TaylorModel{I}, b::TaylorModel{I}) where {I <: Interval}
    ctx = compatible_models(a, b)
    p, q = a._polynomial, b._polynomial
    powers = power_cache(a._coordinates.normalized, DA.degree(p) + DA.degree(q), I)
    tail = discarded_bound(p, q, a._order, powers, I)
    rp, rq = polynomial_bound(p, powers, I), polynomial_bound(q, powers, I)
    r = tail + rp * b._remainder + rq * a._remainder + a._remainder * b._remainder
    enclosure = if a._enclosure === nothing && b._enclosure === nothing
        nothing
    else
        ra = intersect_range(a, rp + a._remainder)
        a === b ? IA.pown(ra, 2) : ra * intersect_range(b, rq + b._remainder)
    end
    retained = DA.allocate(ctx, I, ctx.basis.ends[a._order + 1])
    DA.multiply!(retained, p, q, a._order)
    return model_result(retained, r, a._coordinates, a._order, enclosure)
end
for op in (:*, :/)
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

function inverse_factorial(::Type{I}, k) where {I <: Interval}
    T = IA.numtype(I)
    factor = one(I)
    for j in 1:k
        factor /= IA.interval(T, j)
    end
    return factor
end

# f^(k)(X)/k!, valid on the hull from the expansion point to the input range.
function model_derivative(fn, x::I, k::Int) where {I <: Interval}
    factor = inverse_factorial(I, k)
    fn === exp && return exp(x) * factor
    s, c = sincos(x)
    values = fn === sin ? (s, c, -s, -c) : (c, -s, -c, s)
    return values[mod(k, 4) + 1] * factor
end

function model_horner(a, h, coefficients)
    result = TaylorModel(coefficients[end], a)
    for k in (length(coefficients) - 1):-1:1
        result = result * h + coefficients[k]
    end
    return result
end

# For |u| <= r <= 1 the absolute binomial tail is bounded by
# sqrt(c)*binomial(2n,n)/4^n * r^(n+1). The central factor is also
# (2n-1)*abs(sqrt(c)*binomial(1/2,n)), available from the last coefficient.
function sqrt_remainder(x::I, center::I, coefficients) where {I <: Interval}
    T = IA.numtype(I)
    n = length(coefficients) - 1
    # The chosen center is inside x and at least sup(x)/2, proving |u| <= 1.
    u = IA.intersect_interval((x - center) / center, IA.interval(T, -1, 1); dec = :auto)
    radius = IA.interval(T, IA.mag(u))
    power = IA.pown(radius, n + 1)
    central = IA.interval(T, 2n - 1) * abs(coefficients[end])
    upper = IA.sup(central * power)
    if IA.sup(radius) < 1
        omitted = coefficients[end] * IA.interval(T, 1 // 2 - n) / IA.interval(T, n + 1)
        upper = min(upper, IA.sup(abs(omitted) * power / (one(I) - radius)))
    end
    # Even-order derivatives of sqrt are negative throughout the positive axis.
    return IA.interval(T, -upper, isodd(n) ? zero(T) : upper)
end

# Compose in u=(a-c)/c, avoiding large or tiny powers of c in coefficients.
function relative_model(fn, a::TaylorModel{I}, x::I, center::I) where {I <: Interval}
    T = IA.numtype(I)
    n = a._order
    coefficients = Vector{I}(undef, n + 1)
    coefficients[1] = fn(center)
    for k in 1:n
        coefficients[k + 1] = if fn === inv
            -coefficients[k]
        elseif fn === log
            (isodd(k) ? one(I) : -one(I)) / IA.interval(T, k)
        else
            coefficients[k] * IA.interval(T, 1 // 2 - (k - 1)) / IA.interval(T, k)
        end
    end
    # The displacement owns its polynomial. Divide its coefficients directly:
    # an intermediate reciprocal can overflow at a subnormal center even when
    # every normalized coefficient (and the final sqrt/log result) is finite.
    h = a - center
    p = h._polynomial
    for k in 1:p.len
        p.coeffs[k] /= center
    end
    h = model_result(
        DA.finish!(p, p.len; filter = false), h._remainder / center, a._coordinates, n,
        h._enclosure === nothing ? nothing : h._enclosure / center
    )
    result = model_horner(a, h, coefficients)
    tail = if fn === inv
        IA.pown(-((x - center) / center), n + 1) / x
    elseif fn === log
        # Integral remainder: (-1)^n*u^(n+1)*int_0^1 s^n/(1+s*u) ds.
        # Since x/c contains 1, it bounds every denominator along the segment.
        factor = isodd(n) ? -one(I) : one(I)
        u = (x - center) / center
        bound = factor * IA.pown(u, n + 1) /
            (IA.interval(T, n + 1) * (x / center))
        if IA.isbounded(bound) && IA.decoration(bound) >= IA.def
            bound
        else
            # A rounded x/c may include zero for an extremely wide positive x.
            # Evaluate a separate finite residual instead of using that bound.
            approximation = coefficients[end]
            for k in n:-1:2
                approximation = approximation * u + coefficients[k]
            end
            residual = log(x) - coefficients[1] - approximation * u
            isodd(n) ? IA.intersect_interval(
                    residual, IA.interval(T, IA.inf(residual), 0); dec = :auto
                ) : residual
        end
    else
        sqrt_remainder(x, center, coefficients)
    end
    return model_result(
        result._polynomial, result._remainder + tail, a._coordinates, n, fn(x)
    )
end

function model_center(fn, x::I) where {I <: Interval}
    T = IA.numtype(I)
    midpoint = IA.mid(x)
    if fn === sqrt
        # Outward half-upper also handles rounding when the midpoint is near
        # sup(x)/2. It stays inside x and makes the binomial tail valid at r=1.
        midpoint = max(midpoint, IA.sup(IA.interval(T, IA.sup(x)) / 2))
    end
    center = IA.interval(T, midpoint)
    IA.issubset_interval(center, x) || throw(ArgumentError("Could not construct an expansion point inside the input enclosure"))
    return center
end

function elementary_model(fn, a::TaylorModel{I}) where {I <: Interval}
    model_context(a)
    x = checked_interval(DA.enclose(a); guaranteed = true)
    fn === sqrt && IA.isthinzero(x) && return zero(a)
    if fn === inv && IA.in_interval(0, x)
        throw(DomainError(x, "The reciprocal input enclosure must exclude zero"))
    elseif fn in (log, sqrt) && IA.inf(x) <= 0
        throw(DomainError(x, "This expansion requires a strictly positive whole enclosure"))
    end
    center = model_center(fn, x)
    fn in (inv, log, sqrt) && return relative_model(fn, a, x, center)
    h = a - center
    n = a._order
    coefficients = if fn === exp
        values = Vector{I}(undef, n + 1)
        values[1] = exp(center)
        for k in 1:n
            values[k + 1] = values[k] / DA.degree_factor(values[k], k)
        end
        values
    else
        s, c = sincos(center)
        values = fn === sin ? (s, c, -s, -c) : (c, -s, -c, s)
        DA.cycle_coefficients(values, n)
    end
    result = model_horner(a, h, coefficients)
    # Taylor's theorem: the (n+1)th normalized derivative is bounded on x;
    # h^(n+1) here is an interval power of the whole input displacement.
    tail = model_derivative(fn, x, n + 1) * IA.pown(x - center, n + 1)
    return model_result(
        result._polynomial, result._remainder + tail, a._coordinates, n, fn(x)
    )
end
for fn in (:inv, :exp, :log, :sin, :cos, :sqrt)
    @eval Base.$fn(a::TaylorModel) = elementary_model($fn, a)
end


function Base.sincos(a::TaylorModel{I}) where {I <: Interval}
    model_context(a)
    x = checked_interval(DA.enclose(a); guaranteed = true)
    center = model_center(sin, x)
    h = a - center
    n = a._order
    s0, c0 = sincos(center)
    s, c = DA.trig_coefficients(s0, c0, n)
    sine, cosine = model_horner(a, h, s), model_horner(a, h, c)
    sx, cx = sincos(x)
    sine_cycle, cosine_cycle = (sx, cx, -sx, -cx), (cx, -sx, -cx, sx)
    index = mod(n + 1, 4) + 1
    power = IA.pown(x - center, n + 1)
    factor = inverse_factorial(I, n + 1)
    return (
        model_result(
            sine._polynomial, sine._remainder + sine_cycle[index] * factor * power,
            a._coordinates, n, sx
        ),
        model_result(
            cosine._polynomial, cosine._remainder + cosine_cycle[index] * factor * power,
            a._coordinates, n, cx
        ),
    )
end

function compile_models(models::AbstractVector{TaylorModel{I}}) where {I <: Interval}
    validated_rounding()
    reference = first(models)
    b = DA.model_valid(reference).basis
    foreach(m -> DA.model_valid(m), models)
    # Each patch owns one coordinate family and exponent matrix. Outputs share
    # these read-only snapshots, but own their coefficients and remainders.
    coordinates = deepcopy(reference._coordinates)
    exponents = Matrix{Int}(b.exponents[:, 1:maximum(m -> m._polynomial.len, models)])
    return [
        DA.CompiledTaylorModel(
            deepcopy(m._polynomial.coeffs[1:m._polynomial.len]), exponents,
            deepcopy(m._remainder), deepcopy(m._enclosure), coordinates,
            m._order, DA.degree(m)
        ) for m in models
    ]
end
DA.compile(a::TaylorModel{I}) where {I <: Interval} = only(compile_models([a]))

function snapshot_bounds(
        models::AbstractVector{<:DA.CompiledTaylorModel{I}}, box
    ) where {I <: Interval}
    normalized = normalized_box(first(models), box)
    powers = power_cache(normalized, maximum(DA.degree, models), I)
    return [deepcopy(model_bound(m, powers, I)) for m in models]
end

# Polynomial functions enclose coefficient arithmetic only. They do not acquire
# a Taylor-model remainder. For interval coefficients sqrt avoids interval order
# comparisons and uses an exact enclosed half rather than a floating exponent.
function Base.sqrt(a::TaylorPolynomial{I}) where {I <: Interval}
    x = DA.constant_term(a)
    DA.coefficient_iszero(a) && return zero(a)
    IA.inf(checked_interval(x)) > 0 || throw(
        DA.TaylorError(
            "Polynomial square-root expansion requires a strictly positive constant interval"
        )
    )
    c0 = sqrt(x)
    return isempty(a.algebra.basis.products) ? DA.power_series(a, IA.interval(IA.numtype(I), 1 // 2), c0) : DA.square_root(a, c0)
end
