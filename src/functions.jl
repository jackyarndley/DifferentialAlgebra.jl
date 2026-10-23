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
        second_coeffs::AbstractVector{T}; filter::Bool = true
    ) where {T}
    ctx = valid(a)
    R = promote_type(coefficient_type(a), T)
    if a.len == 1
        return TaylorPolynomial{R}(ctx, first_coeffs[1]; filter),
            TaylorPolynomial{R}(ctx, second_coeffs[1]; filter)
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
    return finish!(s, s.len; filter), finish!(c, c.len; filter)
end

# Normalization is an intermediate analytic step: a small coefficient can grow
# in the subsequent composition. Only the completed public result is filtered.
function normalized_argument(a::TaylorPolynomial, x::Real)
    return divided(a, x; filter = false)
end
function domain_call(f, x)
    return try
        value = f(x)
        value isa Real && isnan(value) && !isnan(x) &&
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
function recurrence(a::TaylorPolynomial, c0, p, mode; filter::Bool = true)
    return isempty(valid(a).basis.products) ?
        recurrence(a, c0, p, mode, Val(false); filter) :
        recurrence(a, c0, p, mode, Val(true); filter)
end

# This formula is shared by stored pairs and on-demand complementary monomials.
function recurrence_term(ai::T, yj::T, a0::T, p::T, di, dj, degree, ::Val{MODE}) where {T, MODE}
    denominator = coefficient_convert(T, degree_factor(ai, degree))
    if MODE === :exp
        return coefficient_product_quotient(
            ai, yj, denominator, coefficient_convert(T, degree_factor(ai, di))
        )
    end
    MODE === :log && dj == 0 && return zero(T)
    left = coefficient_convert(T, degree_factor(ai, di)) / denominator
    right = coefficient_convert(T, degree_factor(ai, dj)) / denominator
    # p*di can overflow even when p*(di/degree) and the result are finite.
    weight = MODE === :log ? -right : p * left - right
    return coefficient_product_quotient(ai, yj, a0, weight)
end

# With coefficients in this conservative range, every weighted product and
# sum stays normal. Ordinary problems can accumulate one numerator per output;
# extreme inputs use the range-preserving contribution formula above.
function recurrence_range(a, c0, p, ::Type{T}, n) where {T}
    T <: Union{Float32, Float64, BigFloat} || return nothing
    upper = sqrt(sqrt(floatmax(T) / (T(a.algebra.basis.order) * T(min(a.len, n)))))
    upper > one(T) || return nothing
    lower = inv(upper)
    normal(x) = iszero(x) || lower <= abs(x) <= upper
    normal(c0) && normal(p) && all(normal, @view(a.coeffs[1:a.len])) || return nothing
    return (lower, upper)
end

function recurrence(
        a::TaylorPolynomial, c0, p, mode::Val{MODE}, table::Val{TABLE};
        filter::Bool = true
    ) where {MODE, TABLE}
    ctx = valid(a); basis = ctx.basis
    T = promote_type(coefficient_type(a), typeof(c0), typeof(p))
    a.len == 1 && return finish!(TaylorPolynomial{T}([coefficient_convert(T, c0)], 1, ctx), 1; filter)
    n = basis.ends[working_order(ctx) + 1]
    out = allocate_undef(ctx, T, n); out.coeffs[1] = c0
    a0, exponent = coefficient_convert(T, a.coeffs[1]), coefficient_convert(T, p)
    powers = TABLE ? Int[] : zeros(Int, basis.variables)
    range = recurrence_range(a, c0, p, T, n)
    direct = range !== nothing
    @inbounds for k in 2:n
        out.coeffs[k] = direct ?
            recurrence_coefficient(a, out, a0, exponent, k, powers, mode, table, Val(:fast)) :
            recurrence_coefficient(a, out, a0, exponent, k, powers, mode, table, Val(:checked))
        if direct && !coefficient_iszero(out.coeffs[k])
            direct = range[1] <= abs(out.coeffs[k]) <= range[2]
        end
    end
    return finish!(out, n; filter)
end

@inline function recurrence_coefficient(
        a, out::TaylorPolynomial{T}, a0, p, k, powers,
        mode::Val{MODE}, table::Val{TABLE}, ::Val{POLICY}
    ) where {T, MODE, TABLE, POLICY}
    basis = a.algebra.basis
    degree = basis.degrees[k]
    # Preserve cancellation in the finite numerator before dividing. Only a
    # range failure retries with separately normalized contributions.
    value = zero(T)
    pairs = TABLE ? ((basis.splits[k] + 1):split_end(basis, k, a.len)) :
        monomial_divisors(basis, k, a.len, powers)
    @inbounds for pair in pairs
        i = TABLE ? Int(basis.left[pair]) : pair[1]
        ai = coefficient_convert(T, a.coeffs[i])
        coefficient_iszero(ai) && continue
        j = TABLE ? Int(basis.right[pair]) : pair[2]
        if POLICY !== :normalized
            di, dj = degree_factor(ai, basis.degrees[i]), degree_factor(ai, basis.degrees[j])
            weight = MODE === :exp ? di : MODE === :log ? -dj : p * di - dj
            term = weight * ai * out.coeffs[j]
            if POLICY === :checked
                primal = coefficient_primal(term)
                if primal isa AbstractFloat &&
                        (!isfinite(primal) || abs(primal) < floatmin(typeof(primal))) &&
                        !coefficient_iszero(weight) && !coefficient_iszero(out.coeffs[j])
                    return recurrence_coefficient(a, out, a0, p, k, powers, mode, table, Val(:normalized))
                end
            end
            value += term
            if POLICY === :checked && !isfinite(value)
                return recurrence_coefficient(a, out, a0, p, k, powers, mode, table, Val(:normalized))
            end
        else
            value += recurrence_term(
                ai, out.coeffs[j], a0, p,
                basis.degrees[i], basis.degrees[j], degree, mode
            )
        end
    end
    if POLICY !== :normalized
        factor = degree_factor(a0, degree)
        MODE === :log && (value += factor * coefficient_convert(T, coefficient_at(a, k)))
        denominator = MODE === :exp ? factor : factor * a0
        if POLICY === :checked && (!isfinite(value) || !isfinite(denominator))
            return recurrence_coefficient(a, out, a0, p, k, powers, mode, table, Val(:normalized))
        end
        return value / denominator
    end
    return MODE === :log ? coefficient_convert(T, coefficient_at(a, k)) / a0 + value : value
end

function trig_recurrence(
        a::TaylorPolynomial, hyperbolic::Bool, ::Val{OUTPUT} = Val(:both);
        filter::Bool = true
    ) where {OUTPUT}
    x = constant_term(a); ctx = a.algebra; basis = ctx.basis
    if isempty(basis.products)
        s, c = trig_series_pair(a, hyperbolic; filter)
        return OUTPUT === :both ? (s, c) : OUTPUT === :sin ? s : c
    end
    s0, c0 = hyperbolic ? (sinh(x), cosh(x)) : sincos(x)
    T = typeof(s0)
    if a.len == 1
        OUTPUT === :both && return (
            TaylorPolynomial{T}(ctx, s0; filter), TaylorPolynomial{T}(ctx, c0; filter),
        )
        return TaylorPolynomial{T}(ctx, OUTPUT === :sin ? s0 : c0; filter)
    end
    order = working_order(ctx)
    n = basis.ends[order + 1]
    # A single result needs its partner only through the preceding degree.
    shared = OUTPUT === :both ? n : basis.ends[order]
    ns, nc = OUTPUT === :cos ? shared : n, OUTPUT === :sin ? shared : n
    s, c = allocate_undef(ctx, T, ns), allocate_undef(ctx, T, nc)
    s.coeffs[1], c.coeffs[1] = s0, c0
    range = recurrence_range(a, s0, c0, T, n)
    direct = range !== nothing
    @inbounds for k in 2:shared
        if !direct
            s.coeffs[k] = recurrence_coefficient(a, c, one(T), zero(T), k, (), Val(:exp), Val(true), Val(:checked))
            cv = recurrence_coefficient(a, s, one(T), zero(T), k, (), Val(:exp), Val(true), Val(:checked))
            c.coeffs[k] = hyperbolic ? cv : -cv
            continue
        end
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
        direct = all(v -> iszero(v) || range[1] <= abs(v) <= range[2], (s.coeffs[k], c.coeffs[k]))
    end
    OUTPUT === :both && return finish!(s, n; filter), finish!(c, n; filter)
    result, partner = OUTPUT === :sin ? (s, c) : (c, s)
    @inbounds for k in (shared + 1):n
        value = direct ?
            recurrence_coefficient(a, partner, one(T), zero(T), k, (), Val(:exp), Val(true), Val(:fast)) :
            recurrence_coefficient(a, partner, one(T), zero(T), k, (), Val(:exp), Val(true), Val(:checked))
        OUTPUT === :cos && !hyperbolic && (value = -value)
        result.coeffs[k] = value
    end
    return finish!(result, n; filter)
end

function exponential(a::TaylorPolynomial; filter::Bool = true)
    c0 = exp(constant_term(a))
    return recurrence(a, c0, zero(c0), Val(:exp); filter)
end
Base.exp(a::TaylorPolynomial) = exponential(a)
function logarithm(a::TaylorPolynomial; filter::Bool = true)
    x = constant_term(a); c0 = domain_call(log, x)
    !invertible_constant(x) && throw(TaylorError("Logarithm requires a constant excluding zero"))
    return recurrence(a, c0, zero(c0), Val(:log); filter)
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
    invertible_constant(constant_term(a)) ||
        throw(TaylorError("Noninteger powers require an expansion point excluding zero"))
    return recurrence(a, c0, p, Val(:power); filter)
end
function Base.:^(a::TaylorPolynomial, p::AbstractFloat)
    R = promote_type(coefficient_type(a), typeof(p))
    isinteger(p) && typemin(Int32) < p <= typemax(Int32) && return convert(TaylorPolynomial{R}, a)^Int(p)
    coefficient_iszero(a) && p > 0 &&
        return TaylorPolynomial{R}(valid(a), zero(R))
    return power_series(a, p, domain_call(x -> x^p, constant_term(a)))
end
function Base.:^(a::TaylorPolynomial, p::Rational)
    R = float(promote_type(coefficient_type(a), typeof(p)))
    isinteger(p) && return convert(TaylorPolynomial{R}, a)^numerator(p)
    coefficient_iszero(a) && p > 0 && return TaylorPolynomial{R}(valid(a), zero(R))
    exponent = coefficient_convert(R, p)
    center = coefficient_convert(R, constant_term(a))
    return power_series(a, exponent, domain_call(x -> x^exponent, center))
end
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
    invertible_constant(base) ||
        throw(TaylorError("A varying exponent requires a nonzero base primal"))
    return exp(scaled(b, domain_call(log, base); filter = false))
end
Base.:^(a::TaylorPolynomial, b::TaylorPolynomial) =
    exp(product(logarithm(a; filter = false), b; filter = false))
Base.:^(::Irrational{:ℯ}, b::TaylorPolynomial) = exp(b)

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
            return divided(a, b0; filter)
        end
        numerator = coefficient_operand(T, a)
        result = numerator / coefficient_convert(T, b0)
        return finish!(TaylorPolynomial{typeof(result)}([result], 1, ctx), 1; filter)
    end
    return isempty(ctx.basis.products) ?
        quotient_recurrence(a, b, T, Val(false); filter) :
        quotient_recurrence(a, b, T, Val(true); filter)
end
function quotient_recurrence(a, b, ::Type{T}, ::Val{TABLE}; filter::Bool = true) where {T, TABLE}
    ctx = b.algebra; basis = ctx.basis
    b0 = coefficient_convert(T, constant_term(b))
    n = basis.ends[working_order(ctx) + 1]
    result = allocate_undef(ctx, T, n)
    result.coeffs[1] = coefficient_convert(T, constant_term(a)) / b0
    powers = TABLE ? Int[] : zeros(Int, basis.variables)
    @inbounds for k in 2:n
        # Subtract before division so a-a can cancel even when a/b0 overflows.
        # Retry in normalized coordinates only if a product loses its range.
        for normalized in (false, true)
            value = coefficient_convert(T, coefficient_at(a, k))
            normalized && (value /= b0)
            rescale = false
            pairs = TABLE ? ((basis.splits[k] + 1):split_end(basis, k, b.len)) :
                monomial_divisors(basis, k, b.len, powers)
            for pair in pairs
                i = TABLE ? Int(basis.left[pair]) : pair[1]
                bi = coefficient_convert(T, b.coeffs[i])
                coefficient_iszero(bi) && continue
                j = TABLE ? Int(basis.right[pair]) : pair[2]
                term = normalized ? coefficient_product_quotient(bi, result.coeffs[j], b0, one(T)) :
                    bi * result.coeffs[j]
                primal = coefficient_primal(term)
                if !normalized && primal isa AbstractFloat &&
                        (!isfinite(primal) || abs(primal) < floatmin(typeof(primal))) &&
                        !coefficient_iszero(result.coeffs[j])
                    rescale = true
                    break
                end
                value -= term
                primal = coefficient_primal(value)
                if !normalized && primal isa AbstractFloat && !isfinite(primal)
                    rescale = true
                    break
                end
            end
            rescale && continue
            result.coeffs[k] = normalized ? value : value / b0
            break
        end
    end
    return finish!(result, n; filter)
end
Base.inv(a::TaylorPolynomial) = quotient(one(coefficient_type(a)), a)

# y^2=a: use each unordered pair of nonconstant monomials only once.
function square_root(a::TaylorPolynomial, c0; filter::Bool = true)
    ctx = valid(a); basis = ctx.basis
    a.len == 1 && return finish!(TaylorPolynomial{typeof(c0)}([c0], 1, ctx), 1; filter)
    isempty(basis.products) && return tablefree_square_root(a, c0; filter)
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

# The same triangular square-root recurrence, ranking complementary exponents
# directly when a multiplication table is unavailable. Used for tiny hypot centers,
# where normalizing before the root would overflow representable derivatives.
function tablefree_square_root(a::TaylorPolynomial, c0; filter::Bool = true)
    ctx = valid(a); basis = ctx.basis
    n = basis.ends[working_order(ctx) + 1]
    result = allocate_undef(ctx, typeof(c0), n)
    result.coeffs[1] = c0
    powers = zeros(Int, basis.variables)
    @inbounds for k in 2:n
        value = coefficient_at(a, k)
        for (i, j) in monomial_divisors(basis, k, k - 1, powers)
            i > j && continue
            term = result.coeffs[i] * result.coeffs[j]
            value -= i == j ? term : term + term
        end
        result.coeffs[k] = value / (c0 + c0)
    end
    return finish!(result, n; filter)
end
function root_series(a::TaylorPolynomial, p::Integer; filter::Bool = true)
    p != 0 || throw(ArgumentError("Zeroth root is undefined"))
    x = constant_term(a)
    !invertible_constant(x) && !(p > 0 && coefficient_iszero(a)) &&
        throw(TaylorError("Root is not analytic at a zero constant"))
    iseven(p) && x < 0 && throw(TaylorError("Even root of a negative constant"))
    exponent = one(float(x)) / p
    c0 = domain_call(
        x -> p == 2 ? sqrt(x) :
            p == 3 ? cbrt(x) : copysign(abs(x)^exponent, x), x
    )
    coefficient_iszero(a) && return finish!(scalar_argument(a, zero(c0)), 1; filter)
    p == 2 && return square_root(a, c0; filter)
    return power_series(a, exponent, c0; filter)
end
"""
    DifferentialAlgebra.nthroot(p::TaylorPolynomial, n = 2)

Expand the real `n`th root; negative `n` returns the reciprocal root and zero
is invalid. Even roots require a nonnegative center; odd roots also allow negative
centers. Except for an identically zero polynomial with positive `n`, the center
must exclude zero for an analytic expansion. The result follows the working order
and scalar root's coefficient type. `n = 1` instead returns a copy preserving all
stored degrees and the original coefficient type.
"""
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
    sine, cosine = (s0, c0, -s0, -c0), (c0, -s0, -c0, s0)
    return cycle_coefficients(sine, n), cycle_coefficients(cosine, n)
end
function trig_series_pair(a::TaylorPolynomial, hyperbolic::Bool = false; filter::Bool = true)
    x = constant_term(a)
    s0, c0 = hyperbolic ? (sinh(x), cosh(x)) : sincos(x)
    T = typeof(s0)
    if a.len == 1
        return TaylorPolynomial{T}(a.algebra, s0; filter), TaylorPolynomial{T}(a.algebra, c0; filter)
    end
    n = working_order(a.algebra)
    s, c = hyperbolic ? (cycle_coefficients((s0, c0, s0, c0), n), cycle_coefficients((c0, s0, c0, s0), n)) :
        trig_coefficients(s0, c0, n)
    return series_pair(a, s, c; filter)
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
    invertible_constant(q0) || throw(TaylorError("Singular derivative at the expansion point"))
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
        a.len == 1 && coefficient_isequal(x, coefficient_primal(x)) &&
            return series(a, [c0]; filter)
        x = coefficient_convert(typeof(c0), x)
        if $(fn === :atan) && x isa AbstractFloat && abs(x) > 2
            reciprocal = quotient(one(x), a; filter = false)
            result = inverse_series(reciprocal, Val(:atan); filter = false)
            scale!(result, result, -one(coefficient_type(result)); filter = false)
            result.coeffs[1] = c0
            return finish!(result, result.len; filter)
        end
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
    !invertible_constant(x0) && !invertible_constant(y0) && throw(TaylorError("atan is not analytic at the origin"))
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
    if !invertible_constant(a0)
        iszero(a) && iszero(b) && return scalar_argument(a, zero(hypot(a0, b0)))
        throw(TaylorError("hypot is not analytic at the origin"))
    end
    iszero(b) && return scaled(a, sign(a0) * one(hypot(a0, b0)))
    scale = abs(a0)
    large = scale isa AbstractFloat && scale > sqrt(floatmax(typeof(scale))) &&
        isfinite(hypot(a0, b0))
    large && (scale = sqrt(hypot(a0, b0)))
    tiny = scale isa AbstractFloat && scale < sqrt(floatmin(typeof(scale))) &&
        max(coefficient_norm(a, 0), coefficient_norm(b, 0)) >
        scale * sqrt(floatmax(typeof(scale)))
    left, right = tiny ? (a, b) :
        (scaled(a, inv(scale); filter = false), scaled(b, inv(scale); filter = false))
    left, right = product(left, left; filter = false), product(right, right; filter = false)
    R = promote_type(coefficient_type(left), coefficient_type(right))
    squares = weighted_sum!(
        allocate_undef(a.algebra, R, max(left.len, right.len)),
        left, one(R), right, one(R); filter = false
    )
    # The scalar hypot preserves a tiny constant whose square underflows.
    tiny && return square_root(squares, hypot(a0, b0))
    # A second normalization would underflow small, representable derivatives
    # at huge centers. Balanced scaling keeps both squares and derivatives finite.
    root = large ? square_root(squares, hypot(a0, b0) / scale; filter = false) :
        root_series(squares, 2; filter = false)
    return scaled(root, scale)
end
Base.hypot(a::TaylorPolynomial, b::Real) = hypot(a, scalar_argument(a, b))
Base.hypot(a::Real, b::TaylorPolynomial) = hypot(scalar_argument(b, a), b)
function Base.mod(a::TaylorPolynomial, p::Real)
    out = TaylorPolynomial{promote_type(coefficient_type(a), typeof(p))}(a)
    out.coeffs[1] = mod(constant_term(a), p)
    limit = out.algebra.basis.ends[working_order(out.algebra) + 1]
    return finish!(out, min(out.len, limit))
end
function Base.mod(a::TaylorPolynomial, p::TaylorPolynomial)
    compatible(a, p)
    p.len == 1 || throw(ArgumentError("A polynomial modulus must be constant"))
    return mod(a, constant_term(p))
end
