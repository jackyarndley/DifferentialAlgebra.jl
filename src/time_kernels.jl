# Shared coefficient recurrences for eager TimeSeries arithmetic and recorded ODEs.
# Coefficients use zero-based mathematical degrees: degree k lives at index k + 1.
# Callers fill earlier degrees before requesting k; exact and polynomial arithmetic
# retain separate kernels where their accumulation and storage rules differ.

# Polynomial coefficients inside a time operation are intermediate values.
# Preserve them until the complete series is returned; a small uncertainty
# coefficient can be amplified by a later time coefficient or recorded node.
time_coefficient(f, args...) = f(args...)
function time_coefficient(f::Union{typeof(+), typeof(-), typeof(*), typeof(/)}, a::Real, b::Real)
    T = promote_type(typeof(a), typeof(b))
    return f(coefficient_convert(T, a), coefficient_convert(T, b))
end
time_muladd(a, b, c) = a isa TaylorPolynomial || b isa TaylorPolynomial || c isa TaylorPolynomial ?
    time_coefficient(+, time_coefficient(*, a, b), c) : coefficient_muladd(a, b, c)
time_identity(x::Real) = one(x)
time_identity(x::TaylorPolynomial{T}) where {T} = TaylorPolynomial{T}(valid(x), one(T); filter = false)
function time_coefficient_sum(f, a::TaylorPolynomial, b::TaylorPolynomial)
    ctx = compatible(a, b)
    T = promote_type(coefficient_type(a), coefficient_type(b))
    f === (-) && (T = negative_coefficient_type(T))
    n = min(max(a.len, b.len), ctx.basis.ends[working_order(ctx) + 1])
    return weighted_sum!(allocate_undef(ctx, T, n), a, one(T), b, f === (+) ? one(T) : -one(T); filter = false)
end
time_coefficient(::typeof(+), a::TaylorPolynomial, b::TaylorPolynomial) = time_coefficient_sum(+, a, b)
time_coefficient(::typeof(-), a::TaylorPolynomial, b::TaylorPolynomial) = time_coefficient_sum(-, a, b)
time_coefficient(::typeof(+), a::TaylorPolynomial, b::Real) = scalar_sum(a, b, Val(:add); filter = false)
time_coefficient(::typeof(+), a::Real, b::TaylorPolynomial) = time_coefficient(+, b, a)
time_coefficient(::typeof(-), a::TaylorPolynomial, b::Real) = scalar_sum(a, b, Val(:subtract); filter = false)
time_coefficient(::typeof(-), a::Real, b::TaylorPolynomial) = scalar_sum(b, a, Val(:subtract_from); filter = false)
time_coefficient(::typeof(-), a::TaylorPolynomial) = scaled(a, -one(negative_coefficient_type(coefficient_type(a))); filter = false)
time_coefficient(::typeof(*), a::TaylorPolynomial, b::TaylorPolynomial) = product(a, b; filter = false)
time_coefficient(::typeof(*), a::TaylorPolynomial, b::Real) = scaled(a, b; filter = false)
time_coefficient(::typeof(*), a::Real, b::TaylorPolynomial) = scaled(b, a; filter = false)
time_coefficient(::typeof(/), a::TaylorPolynomial, b::TaylorPolynomial) = quotient(a, b; filter = false)
time_coefficient(::typeof(/), a::Real, b::TaylorPolynomial) = quotient(a, b; filter = false)
time_coefficient(::typeof(/), a::TaylorPolynomial, b::Real) = divided(a, b; filter = false)
time_coefficient(::typeof(exp), a::TaylorPolynomial) = exponential(a; filter = false)
time_coefficient(::typeof(log), a::TaylorPolynomial) = logarithm(a; filter = false)
time_coefficient(::typeof(sqrt), a::TaylorPolynomial) = root_series(a, 2; filter = false)
time_coefficient(::typeof(cbrt), a::TaylorPolynomial) = root_series(a, 3; filter = false)
function time_coefficient(::typeof(^), a::TaylorPolynomial, p::Integer)
    p >= 0 && return integer_power(a, p; filter = false)
    reciprocal = quotient(one(coefficient_type(a)), a; filter = false)
    return product(integer_power(reciprocal, -(p + 1); filter = false), reciprocal; filter = false)
end
function time_coefficient(::typeof(^), a::TaylorPolynomial, p::Real)
    if p isa Union{AbstractFloat, Rational} && isinteger(p)
        return time_coefficient(^, a, BigInt(p))
    end
    R = promote_type(coefficient_type(a), typeof(p))
    p isa Rational && (R = float(R))
    exponent = coefficient_convert(R, p)
    center = coefficient_convert(R, constant_term(a))
    iszero(a) && p > 0 && return TaylorPolynomial{R}(valid(a), zero(R); filter = false)
    return power_series(a, exponent, domain_call(x -> x^exponent, center); filter = false)
end
function time_coefficient(::typeof(^), a::Real, b::TaylorPolynomial)
    R = promote_type(typeof(a), coefficient_type(b))
    base = coefficient_convert(R, a)
    if coefficient_iszero(base) || b.len == 1
        value = coefficient_iszero(base) ? zero_base_power(base, constant_term(b), iszero(b)) :
            domain_call(x -> float(base)^x, constant_term(b))
        return TaylorPolynomial{typeof(value)}(valid(b), value; filter = false)
    end
    invertible_constant(base) || throw(TaylorError("A varying exponent requires a nonzero base primal"))
    return exponential(scaled(b, domain_call(log, base); filter = false); filter = false)
end
time_coefficient(::typeof(^), a::TaylorPolynomial, b::TaylorPolynomial) =
    exponential(product(logarithm(a; filter = false), b; filter = false); filter = false)
time_coefficient(::typeof(sincos), a::TaylorPolynomial) = trig_recurrence(a, false; filter = false)
for (f, hyperbolic, output) in ((:sin, false, :sin), (:cos, false, :cos), (:sinh, true, :sin), (:cosh, true, :cos))
    @eval function time_coefficient(::typeof($f), a::TaylorPolynomial)
        return trig_recurrence(a, $hyperbolic, Val($(QuoteNode(output))); filter = false)
    end
end

function time_convolution(a, b, k, prototype)
    value = zero(prototype)
    @inbounds for j in 0:k
        aj, bj = a[j + 1], b[k - j + 1]
        (iszero(aj) || iszero(bj)) && continue
        value = muladd(aj, bj, value)
    end
    return value
end
function time_convolution(a::Vector{T}, b::Vector{T}, k, prototype::T) where {T <: Union{Float32, Float64}}
    value = zero(T)
    @inbounds @simd for j in 0:k
        value += a[j + 1] * b[k - j + 1]
    end
    return value
end
function time_convolution(a::Vector{<:TaylorPolynomial}, b::Vector{<:TaylorPolynomial}, k, prototype::P) where {P <: TaylorPolynomial}
    value, scratch = zero(prototype), zero(prototype)
    @inbounds for j in 0:k
        (iszero(a[j + 1]) || iszero(b[k - j + 1])) && continue
        multiply!(scratch, a[j + 1], b[k - j + 1], working_order(value.algebra); filter = false)
        weighted_sum!(value, value, one(coefficient_type(value)), scratch, one(coefficient_type(value)); filter = false)
    end
    return finish!(value, value.len; filter = false)
end

@inline time_weight(j, k, p, unit, ::Val{:exp}) = degree_factor(unit, j) * unit
@inline time_weight(j, k, p, unit, ::Val{:log}) = degree_factor(unit, j - k) * unit
@inline time_weight(j, k, p, unit, ::Val{:power}) =
    coefficient_convert(typeof(unit), p) * (degree_factor(unit, j) * unit) - degree_factor(unit, k - j) * unit
@inline time_weight(j, k, p, unit, ::Val{:quotient}) = -one(unit)

# Normalize scalar exponent weights once, keeping exact coefficient arithmetic
# and the wider working type used by BigInt rationals. Polynomial and dual
# coefficients only supply the scalar type; their values remain untouched.
time_exponent(p, value) = p
time_exponent(p::Rational, value::Real) =
    coefficient_convert(promote_type(typeof(p), typeof(value)), p)
time_exponent(p::Rational, value::TaylorPolynomial) =
    time_exponent(p, constant_term(value))

function time_weighted_sum(a, c, k, p, mode, denominator = time_identity(c[1]))
    return time_divided_sum(a, c, k, p, mode, denominator, zero(c[1]))
end
function time_divided_sum(a, c, k, p, mode, denominator, numerator)
    T = promote_type(eltype(a), eltype(c), typeof(p), typeof(denominator))
    a0, exponent = coefficient_convert(T, denominator), coefficient_convert(T, p)
    for normalized in (false, true)
        value = coefficient_convert(T, numerator)
        normalized && (value /= a0)
        rescale = false
        @inbounds for j in 1:k
            weight = coefficient_convert(T, time_weight(j, k, exponent, one(T) / degree_factor(one(T), k), mode))
            iszero(weight) && continue
            aj, cj = coefficient_convert(T, a[j + 1]), coefficient_convert(T, c[k - j + 1])
            (coefficient_iszero(aj) || coefficient_iszero(cj)) && continue
            term = normalized ? coefficient_product_quotient(aj, cj, a0, weight) : weight * aj * cj
            primal = coefficient_primal(term)
            if !normalized && primal isa AbstractFloat &&
                    (!isfinite(term) || abs(primal) < floatmin(typeof(primal)))
                rescale = true
                break
            end
            value += term
            if !normalized && !isfinite(value)
                rescale = true
                break
            end
        end
        rescale && continue
        return normalized ? value : value / a0
    end
    return
end
function time_weighted_sum(a::Vector{<:TaylorPolynomial}, c::Vector{P}, k, p, mode, denominator = time_identity(c[1]); numerator = nothing) where {P <: TaylorPolynomial}
    time_normal_products(a, c, k, p, denominator) ||
        return time_scaled_sum(a, c, k, p, mode, denominator, numerator)
    value, scratch = numerator === nothing ? zero(c[1]) : P(numerator), zero(c[1])
    unit = one(coefficient_type(value))
    @inbounds for j in 1:k
        weight = time_weight(j, k, p, unit / degree_factor(unit, k), mode)
        (iszero(weight) || iszero(a[j + 1]) || iszero(c[k - j + 1])) && continue
        multiply!(scratch, a[j + 1], c[k - j + 1], working_order(value.algebra); filter = false)
        weighted_sum!(value, value, unit, scratch, weight; filter = false)
    end
    return mode isa Val{:exp} ? finish!(value, value.len; filter = false) : quotient(value, denominator; filter = false)
end
function time_quotient_coefficient(a, b, c, k, denominator)
    return time_divided_sum(b, c, k, 0, Val(:quotient), denominator, a[k + 1])
end
function time_quotient_coefficient(a::Vector{<:TaylorPolynomial}, b::Vector{<:TaylorPolynomial}, c::Vector{P}, k, denominator) where {P <: TaylorPolynomial}
    time_normal_products(b, c, k, 0, denominator) ||
        return time_scaled_sum(b, c, k, 0, Val(:quotient), denominator, a[k + 1])
    value, scratch = P(a[k + 1]), zero(c[1])
    unit = one(coefficient_type(value))
    @inbounds for j in 1:k
        (iszero(b[j + 1]) || iszero(c[k - j + 1])) && continue
        multiply!(scratch, b[j + 1], c[k - j + 1], working_order(value.algebra); filter = false)
        weighted_sum!(value, value, unit, scratch, -unit; filter = false)
    end
    return quotient(value, denominator; filter = false)
end

# Keep the buffered polynomial products for ordinary scales. Extreme floating
# coefficients need normalization inside each scalar product, before it can
# overflow or underflow; dual coefficients use that path for their partials too.
function time_normal_products(a::Vector{<:TaylorPolynomial}, c::Vector{P}, k, p, denominator) where {P <: TaylorPolynomial}
    T = coefficient_type(c[1])
    coefficient_primal(zero(T)) isa AbstractFloat || return true
    T <: AbstractFloat || return false
    n = c[1].algebra.basis.ends[working_order(c[1].algebra) + 1]
    upper = sqrt(sqrt(floatmax(T) / (T(k) * T(n))))
    lower = inv(upper)
    normal(x) = iszero(x) || lower <= abs(x) <= upper
    normal(p) && all(normal, @view(denominator.coeffs[1:denominator.len])) || return false
    return all(1:k) do j
        all(normal, @view(a[j + 1].coeffs[1:a[j + 1].len])) &&
            all(normal, @view(c[k - j + 1].coeffs[1:c[k - j + 1].len]))
    end
end

function time_scaled_sum(a::Vector{<:TaylorPolynomial}, c::Vector{P}, k, p, mode, denominator, numerator = nothing; filter::Bool = false) where {P <: TaylorPolynomial}
    ctx = compatible(c[1], denominator)
    for j in 1:k
        compatible(a[j + 1], c[1])
        compatible(c[k - j + 1], c[1])
    end
    numerator === nothing || compatible(numerator, c[1])
    T = coefficient_type(c[1]); basis = ctx.basis
    n = basis.ends[working_order(ctx) + 1]
    result = allocate_undef(ctx, T, n)
    powers = zeros(Int, basis.variables)
    center, unit = coefficient_convert(T, constant_term(denominator)), one(T)
    @inbounds for r in 1:n
        for normalized in (false, true)
            value = numerator === nothing ? zero(T) : coefficient_convert(T, coefficient_at(numerator, r))
            normalized && (value /= center)
            rescale = false
            for j in 1:k
                weight = coefficient_convert(T, time_weight(j, k, p, unit / degree_factor(unit, k), mode))
                iszero(weight) && continue
                aj, cj = a[j + 1], c[k - j + 1]
                for i in 1:min(aj.len, r)
                    ai = coefficient_convert(T, aj.coeffs[i])
                    coefficient_iszero(ai) && continue
                    q = i == 1 ? r : complementary_index(basis, r, i, powers)
                    1 <= q <= cj.len || continue
                    coefficient_iszero(cj.coeffs[q]) && continue
                    term = normalized ? coefficient_product_quotient(ai, cj.coeffs[q], center, weight) :
                        weight * ai * cj.coeffs[q]
                    primal = coefficient_primal(term)
                    if !normalized && primal isa AbstractFloat &&
                            (!isfinite(term) || abs(primal) < floatmin(typeof(primal)))
                        rescale = true
                        break
                    end
                    value += term
                    if !normalized && !isfinite(value)
                        rescale = true
                        break
                    end
                end
                rescale && break
            end
            rescale && continue
            for i in 2:min(denominator.len, r)
                di = coefficient_convert(T, denominator.coeffs[i])
                coefficient_iszero(di) && continue
                q = complementary_index(basis, r, i, powers)
                q == 0 && continue
                coefficient_iszero(result.coeffs[q]) && continue
                term = normalized ? coefficient_product_quotient(di, result.coeffs[q], center, unit) :
                    di * result.coeffs[q]
                primal = coefficient_primal(term)
                if !normalized && primal isa AbstractFloat &&
                        (!isfinite(term) || abs(primal) < floatmin(typeof(primal)))
                    rescale = true
                    break
                end
                value -= term
                if !normalized && !isfinite(value)
                    rescale = true
                    break
                end
            end
            rescale && continue
            result.coeffs[r] = normalized ? value : value / center
            break
        end
    end
    return finish!(result, n; filter)
end

# A single coefficient of exp(a), log(a), or a^p. The caller checks analyticity
# and supplies a[0] for log/power; exp does not need a denominator.
function time_unary_coefficient(a, c, k, p, denominator, mode::Val{MODE}) where {MODE}
    MODE === :log && return time_divided_sum(a, c, k, p, mode, denominator, a[k + 1])
    return MODE === :exp ? time_weighted_sum(a, c, k, p, mode) :
        time_weighted_sum(a, c, k, p, mode, denominator)
end
function time_unary_coefficient(a::Vector{<:TaylorPolynomial}, c::Vector{P}, k, p, denominator, mode::Val{:log}) where {P <: TaylorPolynomial}
    return time_weighted_sum(a, c, k, p, mode, denominator; numerator = a[k + 1])
end
