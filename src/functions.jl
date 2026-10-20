# Compose univariate Taylor coefficients with the nonconstant part of a TaylorPolynomial.
# Horner stages need only successively increasing total degrees, and reuse buffers.
function series(a::TaylorPolynomial, coefficients::AbstractVector{T}; filter::Bool = true) where {T}
    ctx = valid(a)
    R = promote_type(coefficient_type(a), T)
    a.len == 1 && return finish!(TaylorPolynomial{R}([coefficient_convert(R, coefficients[1])], 1, ctx), 1; filter)
    n = min(working_order(ctx), length(coefficients) - 1)
    h = TaylorPolynomial{R}(a); h.coeffs[1] = zero(R)
    result, scratch = allocate(ctx, R), allocate(ctx, R)
    result.coeffs[1] = coefficients[n + 1]
    for k in (n - 1):-1:0
        multiply!(scratch, h, result, n - k; filter = false)
        scratch.coeffs[1] += coefficients[k + 1]
        finish!(scratch, scratch.len; filter = false)
        result, scratch = scratch, result
    end
    return finish!(result, result.len; filter)
end

# A pair shares the nonconstant argument and rotates three Horner buffers.
function series_pair(
        a::TaylorPolynomial, first_coeffs::AbstractVector{T},
        second_coeffs::AbstractVector{T}
    ) where {T}
    ctx = valid(a)
    R = promote_type(coefficient_type(a), T)
    if a.len == 1
        return TaylorPolynomial{R}(ctx, first_coeffs[1]),
            TaylorPolynomial{R}(ctx, second_coeffs[1])
    end
    n = min(working_order(ctx), length(first_coeffs) - 1, length(second_coeffs) - 1)
    h = TaylorPolynomial{R}(a)
    h.coeffs[1] = zero(R)
    ncoeffs = ctx.basis.ends[n + 1]
    s, c, scratch = allocate(ctx, R, ncoeffs), allocate(ctx, R, ncoeffs),
        allocate(ctx, R, ncoeffs)
    s.coeffs[1], c.coeffs[1] = first_coeffs[n + 1], second_coeffs[n + 1]
    for k in (n - 1):-1:0
        multiply!(scratch, h, s, n - k; filter = false)
        scratch.coeffs[1] += first_coeffs[k + 1]
        finish!(scratch, scratch.len; filter = false)
        s, scratch = scratch, s
        multiply!(scratch, h, c, n - k; filter = false)
        scratch.coeffs[1] += second_coeffs[k + 1]
        finish!(scratch, scratch.len; filter = false)
        c, scratch = scratch, c
    end
    return finish!(s, s.len), finish!(c, c.len)
end

# Normalization is an intermediate analytic step: a small coefficient can grow
# in the subsequent composition. Only the completed public result is filtered.
function normalized_argument(a::TaylorPolynomial, x::Real)
    R = promote_type(coefficient_type(a), typeof(x))
    factor = inv(coefficient_convert(R, x))
    return scale!(allocate_undef(valid(a), typeof(factor), a.len), a, factor; filter = false)
end
function domain_call(f, x)
    return try
        value = f(x)
        value isa Real && isnan(value) &&
            throw(TaylorError("The real expansion is undefined at $x"))
        value
    catch err
        err isa DomainError || rethrow()
        throw(TaylorError(sprint(showerror, err)))
    end
end

# Euler's homogeneous derivative E(x^alpha)=|alpha|x^alpha turns
# E(exp(a))=exp(a)E(a), a E(log(a))=E(a), and a E(a^p)=p a^p E(a)
# into triangular coefficient recurrences. One traversal replaces repeated
# polynomial products, with no cutoff applied to intermediate Taylor factors.
function recurrence(a::TaylorPolynomial, c0, p, ::Val{MODE}; filter::Bool = true) where {MODE}
    ctx = valid(a); basis = ctx.basis
    T = promote_type(coefficient_type(a), typeof(c0), typeof(p))
    a.len == 1 && return finish!(TaylorPolynomial{T}([coefficient_convert(T, c0)], 1, ctx), 1; filter)
    n = basis.ends[working_order(ctx) + 1]
    out = allocate_undef(ctx, T, n); out.coeffs[1] = c0
    a0 = a.coeffs[1]
    @inbounds for k in 2:n
        degree = basis.degrees[k]
        value = zero(T)
        for t in (basis.splits[k] + 1):basis.splits[k + 1]
            i, j = Int(basis.left[t]), Int(basis.right[t])
            i > a.len && break
            ai = a.coeffs[i]
            coefficient_iszero(ai) && continue
            weight = if MODE === :exp
                degree_factor(ai, basis.degrees[i])
            elseif MODE === :log
                degree_factor(ai, -basis.degrees[j])
            else
                p * degree_factor(ai, basis.degrees[i]) - degree_factor(ai, basis.degrees[j])
            end
            value += weight * ai * out.coeffs[j]
        end
        if MODE === :log
            value += k <= a.len ? degree_factor(a.coeffs[k], degree) * a.coeffs[k] : zero(T)
        end
        out.coeffs[k] = if MODE === :exp
            value / degree_factor(value, degree)
        else
            value / (degree_factor(a0, degree) * a0)
        end
    end
    return finish!(out, n; filter)
end

function trig_recurrence(
        a::TaylorPolynomial, hyperbolic::Bool, ::Val{OUTPUT} = Val(:both)
    ) where {OUTPUT}
    x = constant_term(a); ctx = a.algebra; basis = ctx.basis
    s0, c0 = hyperbolic ? (sinh(x), cosh(x)) : sincos(x)
    T = typeof(s0)
    if a.len == 1
        OUTPUT === :both && return (
            TaylorPolynomial{T}(ctx, s0), TaylorPolynomial{T}(ctx, c0),
        )
        return TaylorPolynomial{T}(ctx, OUTPUT === :sin ? s0 : c0)
    end
    order = working_order(ctx)
    n = basis.ends[order + 1]
    # A single result needs its partner only through the preceding degree.
    shared = OUTPUT === :both ? n : basis.ends[order]
    ns, nc = OUTPUT === :cos ? shared : n, OUTPUT === :sin ? shared : n
    s, c = allocate_undef(ctx, T, ns), allocate_undef(ctx, T, nc)
    s.coeffs[1], c.coeffs[1] = s0, c0
    @inbounds for k in 2:shared
        sv, cv = zero(T), zero(T)
        for t in (basis.splits[k] + 1):basis.splits[k + 1]
            i, j = Int(basis.left[t]), Int(basis.right[t])
            i > a.len && break
            ai = a.coeffs[i]
            coefficient_iszero(ai) && continue
            weight = degree_factor(ai, basis.degrees[i]) * ai
            sv += weight * c.coeffs[j]
            cv += weight * s.coeffs[j]
        end
        s.coeffs[k] = sv / degree_factor(sv, basis.degrees[k])
        c.coeffs[k] = (hyperbolic ? cv : -cv) / degree_factor(cv, basis.degrees[k])
    end
    OUTPUT === :both && return finish!(s, n), finish!(c, n)
    result, partner = OUTPUT === :sin ? (s, c) : (c, s)
    @inbounds for k in (shared + 1):n
        value = zero(T)
        for t in (basis.splits[k] + 1):basis.splits[k + 1]
            i, j = Int(basis.left[t]), Int(basis.right[t])
            i > a.len && break
            ai = a.coeffs[i]
            coefficient_iszero(ai) && continue
            value += degree_factor(ai, basis.degrees[i]) * ai * partner.coeffs[j]
        end
        OUTPUT === :cos && !hyperbolic && (value = -value)
        result.coeffs[k] = value / degree_factor(value, basis.degrees[k])
    end
    return finish!(result, n)
end

function exponential(a::TaylorPolynomial; filter::Bool = true)
    x = constant_term(a); c0 = exp(x)
    !isempty(a.algebra.basis.products) &&
        return recurrence(a, c0, zero(c0), Val(:exp); filter)
    c = Vector{typeof(c0)}(undef, working_order(a.algebra) + 1); c[1] = c0
    @inbounds for k in 1:(length(c) - 1)
        c[k + 1] = c[k] / degree_factor(c[k], k)
    end
    return series(a, c; filter)
end
Base.exp(a::TaylorPolynomial) = exponential(a)
function logarithm(a::TaylorPolynomial; filter::Bool = true)
    x = constant_term(a); c0 = domain_call(log, x)
    !invertible_constant(x) && throw(TaylorError("Logarithm requires a constant excluding zero"))
    !isempty(a.algebra.basis.products) &&
        return recurrence(a, c0, zero(c0), Val(:log); filter)
    c = Vector{typeof(c0)}(undef, working_order(a.algebra) + 1); c[1] = c0
    factor = one(x)
    @inbounds for k in 1:(length(c) - 1)
        c[k + 1] = factor / degree_factor(factor, k)
        factor = -factor
    end
    return series(normalized_argument(a, x), c; filter)
end
Base.log(a::TaylorPolynomial) = logarithm(a)
Base.log2(a::TaylorPolynomial) = log(2, a)
Base.log10(a::TaylorPolynomial) = log(10, a)
function Base.log(b::Real, a::TaylorPolynomial)
    R = promote_type(typeof(b), coefficient_type(a))
    return logarithm(a; filter = false) / log(coefficient_convert(R, b))
end
Base.log(::Irrational{:ℯ}, a::TaylorPolynomial) = log(a)

function power_series(a::TaylorPolynomial, p::Real, c0; filter::Bool = true)
    x = constant_term(a)
    coefficient_iszero(x) && throw(TaylorError("Noninteger powers require a nonzero expansion point"))
    !isempty(a.algebra.basis.products) && return recurrence(a, c0, p, Val(:power); filter)
    T = promote_type(coefficient_type(a), typeof(c0), typeof(p))
    x = coefficient_convert(T, x)
    c = Vector{T}(undef, working_order(a.algebra) + 1); c[1] = c0
    @inbounds for k in 1:(length(c) - 1)
        c[k + 1] = c[k] * (p - degree_factor(p, k - 1)) / degree_factor(c[k], k)
    end
    return series(normalized_argument(a, x), c; filter)
end
function Base.:^(a::TaylorPolynomial, p::AbstractFloat)
    R = promote_type(coefficient_type(a), typeof(p))
    isinteger(p) && typemin(Int32) < p <= typemax(Int32) && return convert(TaylorPolynomial{R}, a)^Int(p)
    coefficient_iszero(a) && p > 0 &&
        return TaylorPolynomial{R}(valid(a), zero(R))
    return power_series(a, p, domain_call(x -> x^p, constant_term(a)))
end
Base.:^(a::TaylorPolynomial, p::Rational) = a^convert(float(promote_type(coefficient_type(a), typeof(p))), p)
# A zero base is locally constant for positive exponents. At a varying zero
# exponent the two-sided real function has no analytic expansion.
function zero_base_power(base, exponent, identically_zero)
    exponent > 0 && return zero(float(base))
    identically_zero && return one(float(base))
    throw(TaylorError("A zero base needs a positive exponent or an identically zero exponent"))
end
function Base.:^(a::Real, b::TaylorPolynomial)
    R = promote_type(typeof(a), coefficient_type(b))
    base = coefficient_convert(R, a)
    if coefficient_iszero(base)
        c0 = zero_base_power(base, constant_term(b), iszero(b))
        return finish!(scalar_argument(b, c0), 1)
    end
    if b.len == 1
        c0 = domain_call(x -> float(base)^x, constant_term(b))
        return finish!(scalar_argument(b, c0), 1)
    end
    return exp(scaled(b, domain_call(log, base); filter = false))
end
Base.:^(a::TaylorPolynomial, b::TaylorPolynomial) =
    exp(product(logarithm(a; filter = false), b; filter = false))
Base.:^(::Irrational{:ℯ}, b::TaylorPolynomial) = exp(b)
@inline coefficient_at(a::TaylorPolynomial, k::Int) = k <= a.len ? a.coeffs[k] : zero(coefficient_type(a))
@inline coefficient_at(a::Real, k::Int) = k == 1 ? a : zero(a)

# b*y=a is triangular by total degree. Division needs one recurrence, avoiding
# both the homogeneous-degree weights of a generic power and a second product.
function quotient(a::Real, b::TaylorPolynomial; filter::Bool = true)
    ctx = a isa TaylorPolynomial ? compatible(a, b) : valid(b)
    b0 = constant_term(b)
    !invertible_constant(b0) && throw(TaylorError("Division requires a denominator constant excluding zero"))
    T = promote_type(
        a isa TaylorPolynomial ? coefficient_type(a) : typeof(a), coefficient_type(b)
    )
    T = negative_coefficient_type(typeof(one(T) / one(T)))
    if b.len == 1
        if a isa TaylorPolynomial
            return scaled(a, inv(coefficient_convert(T, b0)); filter)
        end
        numerator = coefficient_operand(T, a)
        result = numerator / coefficient_convert(T, b0)
        return finish!(TaylorPolynomial{typeof(result)}([result], 1, ctx), 1; filter)
    end
    basis = ctx.basis
    if isempty(basis.products)
        reciprocal = power_series(
            b, -one(T), inv(coefficient_convert(T, b0)); filter = false
        )
        return a isa TaylorPolynomial ? product(a, reciprocal; filter) :
            scaled(reciprocal, a; filter)
    end
    n = basis.ends[working_order(ctx) + 1]
    result = allocate_undef(ctx, T, n)
    result.coeffs[1] = coefficient_operand(T, constant_term(a)) / b0
    @inbounds for k in 2:n
        value = coefficient_convert(T, coefficient_at(a, k))
        @simd for t in (basis.splits[k] + 1):split_end(basis, k, b.len)
            i, j = Int(basis.left[t]), Int(basis.right[t])
            bi = b.coeffs[i]
            value -= coefficient_iszero(bi) ? zero(T) : bi * result.coeffs[j]
        end
        result.coeffs[k] = value / b0
    end
    return finish!(result, n; filter)
end
Base.inv(a::TaylorPolynomial) = quotient(one(coefficient_type(a)), a)

# y^2=a: use each unordered pair of nonconstant monomials only once.
function square_root(a::TaylorPolynomial, c0; filter::Bool = true)
    ctx = valid(a); basis = ctx.basis
    a.len == 1 && return finish!(TaylorPolynomial{typeof(c0)}([c0], 1, ctx), 1; filter)
    n = basis.ends[working_order(ctx) + 1]
    result = allocate_undef(ctx, typeof(c0), n)
    result.coeffs[1] = c0
    denominator = c0 + c0
    @inbounds for k in 2:n
        value = coefficient_at(a, k)
        # Mirrored pairs are ordered oppositely; the first half has i <= j
        # and excludes the terminal (k, 1) pair.
        stop = basis.splits[k] + (basis.splits[k + 1] - basis.splits[k]) ÷ 2
        @simd for t in (basis.splits[k] + 1):stop
            i, j = Int(basis.left[t]), Int(basis.right[t])
            term = result.coeffs[i] * result.coeffs[j]
            value -= i == j ? term : term + term
        end
        result.coeffs[k] = value / denominator
    end
    return finish!(result, n; filter)
end
function root_series(a::TaylorPolynomial, p::Integer; filter::Bool = true)
    p != 0 || throw(ArgumentError("Zeroth root is undefined"))
    x = constant_term(a)
    coefficient_iszero(x) && !(p > 0 && coefficient_iszero(a)) &&
        throw(TaylorError("Root is not analytic at a zero constant"))
    iseven(p) && x < 0 && throw(TaylorError("Even root of a negative constant"))
    exponent = one(float(x)) / p
    c0 = domain_call(
        x -> p == 2 ? sqrt(x) :
            p == 3 ? cbrt(x) : copysign(abs(x)^exponent, x), x
    )
    coefficient_iszero(a) && return finish!(scalar_argument(a, c0), 1; filter)
    p == 2 && !isempty(a.algebra.basis.products) && return square_root(a, c0; filter)
    return power_series(a, exponent, c0; filter)
end
# The identity root preserves the stored coefficient type and degrees. Other
# roots use the scalar root's working type, including at the zero polynomial.
nthroot(a::TaylorPolynomial, p::Integer = 2) = p == 1 ? copy(a) : root_series(a, p)
Base.sqrt(a::TaylorPolynomial) = root_series(a, 2)
Base.cbrt(a::TaylorPolynomial) = root_series(a, 3)
function nthroot(x::Real, p::Integer = 2)
    p != 0 || throw(ArgumentError("Zeroth root is undefined"))
    p == 1 && return x
    exponent = one(float(x)) / p
    return isodd(p) ? copysign(abs(x)^exponent, x) : x^exponent
end

for fn in (:sin, :cos, :sinh, :cosh)
    @eval function Base.$fn(a::TaylorPolynomial)
        ctx = valid(a)
        if a.len == 1
            c0 = $fn(a.coeffs[1])
            return TaylorPolynomial{typeof(c0)}(ctx, c0)
        end
        if !isempty(ctx.basis.products)
            return trig_recurrence(
                a, $(fn in (:sinh, :cosh)),
                $(fn in (:sin, :sinh) ? :(Val(:sin)) : :(Val(:cos)))
            )
        end
        x = constant_term(a)
        s, c = $(fn in (:sin, :cos) ? :(sincos(x)) : :((sinh(x), cosh(x))))
        values = $(fn == :sin ? :((s, c, -s, -c)) : fn == :cos ? :((c, -s, -c, s)) : fn == :sinh ? :((s, c, s, c)) : :((c, s, c, s)))
        return series(a, cycle_coefficients(values, working_order(a.algebra)))
    end
end
function cycle_coefficients(values, n)
    coefficients = Vector{typeof(values[1])}(undef, n + 1)
    factor = one(values[1])
    @inbounds for k in 0:n
        coefficients[k + 1] = values[mod(k, 4) + 1] * factor
        factor /= degree_factor(factor, k + 1)
    end
    return coefficients
end
function trig_coefficients(s0, c0, n)
    T = typeof(s0)
    s, c = Vector{T}(undef, n + 1), Vector{T}(undef, n + 1)
    sine, cosine = (s0, c0, -s0, -c0), (c0, -s0, -c0, s0)
    factor = one(s0)
    @inbounds for k in 0:n
        index = mod(k, 4) + 1
        s[k + 1], c[k + 1] = sine[index] * factor, cosine[index] * factor
        factor /= degree_factor(factor, k + 1)
    end
    return s, c
end
function trig_series_pair(a::TaylorPolynomial)
    s0, c0 = sincos(constant_term(a))
    T = typeof(s0)
    if a.len == 1
        return TaylorPolynomial{T}(a.algebra, s0), TaylorPolynomial{T}(a.algebra, c0)
    end
    s, c = trig_coefficients(s0, c0, working_order(a.algebra))
    return series_pair(a, s, c)
end
Base.sincos(a::TaylorPolynomial) = isempty(valid(a).basis.products) ?
    trig_series_pair(a) : trig_recurrence(a, false)
for (fn, sign) in ((:tan, 1), (:tanh, -1))
    @eval function Base.$fn(a::TaylorPolynomial)
        x = constant_term(a); c0 = $fn(x)
        c = zeros(typeof(c0), working_order(a.algebra) + 1); c[1] = c0
        @inbounds for k in 0:(length(c) - 2)
            value = k == 0 ? one(c0) : zero(c0)
            for j in 0:k
                value += $sign * c[j + 1] * c[k - j + 1]
            end
            c[k + 2] = value / (k + 1)
        end
        return series(a, c)
    end
end

# If y=(q0+q1*t+q2*t^2)^p, q*y'=p*q'*y gives this linear recurrence.
function quadratic_power(q0, q1, q2, p, n)
    coefficient_iszero(q0) && throw(TaylorError("Singular derivative at the expansion point"))
    c0 = domain_call(x -> x^p, q0)
    c = zeros(typeof(c0), n + 1); c[1] = c0
    @inbounds for k in 0:(n - 1)
        value = (p - k) * q1 * c[k + 1]
        k > 0 && (value += (2p - k + 1) * q2 * c[k])
        c[k + 2] = value / ((k + 1) * q0)
    end
    return c
end
function integral_series(a::TaylorPolynomial, c0, derivative; filter::Bool = true)
    T = promote_type(typeof(c0), eltype(derivative))
    c = Vector{T}(undef, length(derivative) + 1); c[1] = c0
    @inbounds for k in eachindex(derivative)
        c[k + 1] = derivative[k] / k
    end
    return series(a, c; filter)
end
for fn in (:asin, :acos, :atan, :asinh, :acosh, :atanh)
    @eval function inverse_series(
            a::TaylorPolynomial, ::Val{$(QuoteNode(fn))}; filter::Bool = true
        )
        x = constant_term(a); u = one(float(x)); c0 = domain_call($fn, x)
        a.len == 1 && return series(a, [c0]; filter)
        x = coefficient_convert(typeof(c0), x)
        # Scale floating centers before squaring. Compose in the normalized
        # coordinate so tiny derivatives can still multiply large input slopes.
        if $(fn in (:asinh, :acosh)) && x isa AbstractFloat && abs(x) > 2
            radius = $(
                fn === :asinh ? :(hypot(x, u)) :
                    :(x * sqrt((u - inv(x)) * (u + inv(x))))
            )
            derivative = quadratic_power(
                u, 2(x / radius), u, -u / 2,
                working_order(a.algebra) - 1
            )
            return integral_series(normalized_argument(a, radius), c0, derivative; filter)
        end
        q0, q1, q2, p = $(
            fn in (:asin, :acos) ? :(((u - x) * (u + x), -2x, -u, -u / 2)) :
                fn == :atan ? :((u + x * x, 2x, u, -u)) :
                fn == :asinh ? :((u + x * x, 2x, u, -u / 2)) :
                fn == :acosh ? :(((x - u) * (x + u), 2x, u, -u / 2)) :
                :(((u - x) * (u + x), -2x, -u, -u))
        )
        derivative = quadratic_power(q0, q1, q2, p, working_order(a.algebra) - 1)
        $(fn == :acos) && (derivative .*= -one(eltype(derivative)))
        return integral_series(a, c0, derivative; filter)
    end
    @eval Base.$fn(a::TaylorPolynomial) = inverse_series(a, Val($(QuoteNode(fn))))
end
function Base.atan(y::TaylorPolynomial, x::TaylorPolynomial)
    compatible(y, x)
    x0, y0 = constant_term(x), constant_term(y)
    coefficient_iszero(x0) && coefficient_iszero(y0) && throw(TaylorError("atan is not analytic at the origin"))
    use_x = abs(x0) >= abs(y0)
    ratio = use_x ? quotient(y, x; filter = false) : quotient(x, y; filter = false)
    p = inverse_series(ratio, Val(:atan); filter = false)
    !use_x && scale!(p, p, -one(coefficient_type(p)); filter = false)
    p.coeffs[1] = atan(y0, x0)
    return finish!(p, p.len)
end
# Scalar promotion is internal to these operations and must not filter the input.
function scalar_argument(a::TaylorPolynomial, b::Real)
    R = promote_type(coefficient_type(a), typeof(b))
    return TaylorPolynomial{R}([coefficient_convert(R, b)], 1, valid(a))
end
Base.atan(a::TaylorPolynomial, b::Real) = atan(a, scalar_argument(a, b))
Base.atan(a::Real, b::TaylorPolynomial) = atan(scalar_argument(b, a), b)
function Base.hypot(a::TaylorPolynomial, b::TaylorPolynomial)
    compatible(a, b)
    a0, b0 = constant_term(a), constant_term(b)
    return abs(b0) > abs(a0) ? dominant_hypot(b, a, b0, a0) :
        dominant_hypot(a, b, a0, b0)
end
# Keep each coordinate's concrete type when selecting the dominant argument.
function dominant_hypot(a::TaylorPolynomial, b::TaylorPolynomial, a0, b0)
    if coefficient_iszero(a0)
        iszero(a) && iszero(b) && return scalar_argument(a, hypot(a0, b0))
        throw(TaylorError("hypot is not analytic at the origin"))
    end
    iszero(b) && return scaled(a, sign(a0) * one(hypot(a0, b0)))
    # The dominant coordinate stays outside the root: dividing both inputs by
    # a tiny center can overflow their squared slopes before rescaling them.
    ratio = quotient(b, a; filter = false)
    squares = product(ratio, ratio; filter = false)
    squares.coeffs[1] += one(coefficient_type(squares))
    factor = root_series(squares, 2; filter = false)
    return scaled(product(a, factor; filter = false), sign(a0))
end
Base.hypot(a::TaylorPolynomial, b::Real) = hypot(a, scalar_argument(a, b))
Base.hypot(a::Real, b::TaylorPolynomial) = hypot(scalar_argument(b, a), b)
function Base.mod(a::TaylorPolynomial, p::Real)
    out = TaylorPolynomial{promote_type(coefficient_type(a), typeof(p))}(a)
    out.coeffs[1] = mod(constant_term(a), p)
    return finish!(out, out.len)
end

for fn in (:erf, :erfc)
    @eval function SpecialFunctions.$fn(a::TaylorPolynomial)
        x = constant_term(a); c0 = SpecialFunctions.$fn(x)
        n = working_order(a.algebra)
        derivative = zeros(typeof(c0), n)
        derivative[1] = $(fn == :erf ? 2 : -2) * exp(-x * x) / sqrt(convert(typeof(c0), π))
        @inbounds for k in 1:(n - 1)
            derivative[k + 1] = (-2x * derivative[k] - (k > 1 ? 2derivative[k - 1] : zero(c0))) / k
        end
        return integral_series(a, c0, derivative)
    end
end
function derivative_series(f, a::TaylorPolynomial, c0; filter::Bool = true)
    a.len == 1 && return series(a, [c0]; filter)
    x = constant_term(a)
    c = Vector{typeof(c0)}(undef, working_order(a.algebra) + 1); c[1] = c0
    factor = one(c0)
    for k in 1:(length(c) - 1)
        factor /= k
        c[k + 1] = f(k, x) * factor
    end
    return series(a, c; filter)
end
function gamma_argument(a::TaylorPolynomial)
    x = constant_term(a)
    x <= 0 && coefficient_isinteger(x) &&
        throw(TaylorError("Gamma-family expansion has a pole at $x"))
    return x
end
gamma_series(a::TaylorPolynomial, c0; filter::Bool = true) =
    derivative_series((k, x) -> scalar_psi(k - 1, x), a, c0; filter)
function SpecialFunctions.loggamma(a::TaylorPolynomial)
    return gamma_series(a, domain_call(SpecialFunctions.loggamma, gamma_argument(a)))
end
function SpecialFunctions.gamma(a::TaylorPolynomial)
    value, sign = domain_call(SpecialFunctions.logabsgamma, gamma_argument(a))
    return scaled(exponential(gamma_series(a, value; filter = false); filter = false), sign)
end
function polygamma_series(a::TaylorPolynomial, n::Integer)
    n >= 0 || throw(ArgumentError("Invalid polygamma order"))
    x = gamma_argument(a); c0 = domain_call(x -> scalar_psi(n, x), x)
    return derivative_series((k, x) -> scalar_psi(n + k, x), a, c0)
end
SpecialFunctions.polygamma(n::Integer, a::TaylorPolynomial) = polygamma_series(a, n)
SpecialFunctions.digamma(a::TaylorPolynomial) = polygamma_series(a, 0)

# Neighboring-order derivative identities keep the regular J/I solutions stable
# near zero. Each stage forms Taylor coefficients directly, avoiding large
# binomial/factorial intermediates and reusing the scalar seeds.
function regular_bessel_series(fn::F, n, a, c0, s) where {F}
    x = constant_term(a)
    coefficient_iszero(x) && fn === SpecialFunctions.besselix &&
        throw(TaylorError("Scaled Bessel I is not analytic at zero"))
    cutoff = working_order(a.algebra)
    values = [
        k == 0 ? c0 : scalar_bessel(fn, n + k, x) for k in (-cutoff):cutoff
    ]
    c = Vector{typeof(c0)}(undef, cutoff + 1)
    c[1] = c0
    subtract = fn === SpecialFunctions.besselj
    for k in 1:cutoff
        for j in 1:(2(cutoff - k) + 1)
            neighbor = subtract ? (values[j] - values[j + 2]) / 2 :
                values[j] / 2 + values[j + 2] / 2
            values[j] = (neighbor - s * values[j + 1]) / k
        end
        c[k + 1] = values[cutoff - k + 1]
    end
    return series(a, c)
end

# DLMF 10.32.9 gives e^x*K_n(x) as a positive integral. Its k-th derivative
# inserts (1-cosh(t))^k, avoiding cancellation against the exponential factor.
# Normalize the polynomial argument by x and leave the scalar seed outside the
# composition so tiny intermediate derivatives can still multiply large slopes.
function large_besselk_series(n, a, c0::T) where {T <: AbstractFloat}
    x = coefficient_convert(T, constant_term(a))
    count = working_order(a.algebra) + 1
    nu = T(n)
    root = sqrt(x)
    peak = asinh(abs(nu) / x)
    width = sqrt(max(hypot(x, nu), one(T)))
    logtwo = log(T(2))
    logcosh(t) = abs(t) + log1p(exp(-2abs(t))) - logtwo
    decay(t) = 2(root * sinh(t / 2))^2
    shift = -decay(peak) + logcosh(nu * peak)
    function integrand!(values, u)
        t = peak + u / width
        z = decay(t)
        if !isfinite(z)
            fill!(values, zero(T))
            return nothing
        end
        values[1] = exp(-z + logcosh(nu * t) - shift) / width
        for k in 1:(count - 1)
            values[k + 1] = -z * values[k] / k
        end
        return nothing
    end
    quadrature_order = max(15, cld(precision(T), 8))
    limits = (-peak * width, zero(T), T(Inf))
    approximate, _ = QuadGK.quadgk!(
        integrand!, zeros(T, count), limits...;
        order = quadrature_order, rtol = sqrt(eps(T))
    )
    # Relative component weights also resolve orders whose derivatives are much
    # smaller than the value or neighboring derivatives.
    scales = map(c -> iszero(c) ? one(T) : abs(c), approximate)
    relative_norm(v) = maximum(i -> abs(v[i]) / scales[i], eachindex(v))
    coefficients, error = QuadGK.quadgk!(
        integrand!, zeros(T, count), limits...;
        order = quadrature_order, rtol = 8eps(T), norm = relative_norm
    )
    error <= 64eps(T) * relative_norm(coefficients) ||
        throw(TaylorError("Scaled Bessel K derivative quadrature did not converge"))
    coefficients ./= first(coefficients)
    result = series(normalized_argument(a, x), coefficients; filter = false)
    return scaled(result, c0)
end

# Bessel's second-order differential equation supplies the singular Y/K solutions.
# With y=exp(s*x)*z, scaled I/K remain scaled throughout (no overflow-prone unscaling).
function bessel_series(fn::F, n::Integer, a::TaylorPolynomial, q::Int, s) where {F}
    typemin(Int32) < n <= typemax(Int32) || throw(ArgumentError("Bessel order out of range"))
    x = constant_term(a)
    coefficient_iszero(x) &&
        fn in (
        SpecialFunctions.bessely, SpecialFunctions.besselk, SpecialFunctions.besselkx,
    ) &&
        throw(TaylorError("Bessel expansion is singular at zero"))
    c0 = domain_call(x -> scalar_bessel(fn, n, x), x)
    a.len == 1 && return TaylorPolynomial{typeof(float(c0))}(valid(a), c0)
    if fn in (SpecialFunctions.besselj, SpecialFunctions.besseli, SpecialFunctions.besselix)
        return regular_bessel_series(fn, n, a, c0, s)
    end
    fn === SpecialFunctions.besselkx && x isa AbstractFloat && isfinite(x) && x > 64 &&
        return large_besselk_series(n, a, c0)
    cutoff = working_order(a.algebra)
    coefficient_iszero(x) && throw(TaylorError("Bessel expansion is singular at zero"))
    c = zeros(typeof(c0), cutoff + 1); c[1] = c0
    lower, upper = scalar_bessel(fn, n - 1, x), scalar_bessel(fn, n + 1, x)
    neighbor = fn in (SpecialFunctions.besselj, SpecialFunctions.bessely) ? (lower - upper) / 2 :
        fn in (SpecialFunctions.besseli, SpecialFunctions.besselix) ? (lower + upper) / 2 : -(lower + upper) / 2
    c[2] = neighbor - s * c0
    r = s * s + q
    @inbounds for k in 0:(cutoff - 2)
        value = -(2x * k + x + 2s * x * x) * (k + 1) * c[k + 2] - (k * k + 4s * x * k + r * x * x + s * x - n * n) * c[k + 1]
        k >= 1 && (value -= (2s * (k - 1) + 2r * x + s) * c[k])
        k >= 2 && (value -= r * c[k - 1])
        c[k + 3] = value / (x * x * (k + 2) * (k + 1))
    end
    return series(a, c)
end
for (fn, q) in ((:besselj, 1), (:bessely, 1), (:besseli, -1), (:besselk, -1), (:besselix, -1), (:besselkx, -1))
    @eval SpecialFunctions.$fn(n::Integer, a::TaylorPolynomial) = bessel_series(
        SpecialFunctions.$fn, n, a, $q,
        $(fn == :besselix ? :(sign(constant_term(a))) : fn == :besselkx ? :(-one(constant_term(a))) : :(zero(constant_term(a))))
    )
end
