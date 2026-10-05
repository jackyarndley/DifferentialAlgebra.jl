function Base.copy(a::TaylorPolynomial{T}) where {T <: Real}
    ctx = valid(a)
    return TaylorPolynomial{T}(a.coeffs[1:a.len], a.len, ctx)
end
function TaylorPolynomial{T}(a::TaylorPolynomial) where {T <: Real}
    ctx = valid(a)
    return TaylorPolynomial{T}([coefficient_convert(T, a.coeffs[i]) for i in 1:a.len], a.len, ctx)
end
Base.convert(::Type{TaylorPolynomial}, a::TaylorPolynomial) = a
Base.convert(::Type{TaylorPolynomial}, x::Real) = TaylorPolynomial(x)
Base.convert(::Type{TaylorPolynomial{T}}, a::TaylorPolynomial{T}) where {T <: Real} = a
Base.convert(::Type{TaylorPolynomial{T}}, a::TaylorPolynomial) where {T <: Real} = TaylorPolynomial{T}(a)
Base.convert(::Type{TaylorPolynomial{T}}, x::Real) where {T <: Real} = TaylorPolynomial{T}(x)
Base.promote_rule(::Type{TaylorPolynomial{T}}, ::Type{TaylorPolynomial{S}}) where {T <: Real, S <: Real} = TaylorPolynomial{promote_type(T, S)}
Base.promote_rule(::Type{TaylorPolynomial{T}}, ::Type{S}) where {T <: Real, S <: Real} = TaylorPolynomial{promote_type(T, S)}
Base.promote_rule(::Type{TaylorPolynomial}, ::Type{S}) where {S <: Real} = TaylorPolynomial
# BigFloat otherwise promotes every Real to BigFloat, conflicting with TaylorPolynomial's
# rule and recursively asking Julia to promote TaylorPolynomial{BigFloat} with BigFloat.
Base.promote_rule(::Type{BigFloat}, ::Type{TaylorPolynomial{T}}) where {T <: Real} = TaylorPolynomial{promote_type(BigFloat, T)}
Base.promote_rule(::Type{BigFloat}, ::Type{TaylorPolynomial}) = TaylorPolynomial
Base.promote_rule(::Type{TaylorPolynomial{T}}, ::Type{TaylorPolynomial}) where {T <: Real} = TaylorPolynomial
Base.zero(::Type{TaylorPolynomial}) = TaylorPolynomial()
Base.one(::Type{TaylorPolynomial}) = TaylorPolynomial(1.0)
Base.zero(::Type{TaylorPolynomial{T}}) where {T <: Real} = TaylorPolynomial{T}(zero(T))
Base.one(::Type{TaylorPolynomial{T}}) where {T <: Real} = TaylorPolynomial{T}(one(T))
Base.zero(a::TaylorPolynomial{T}) where {T <: Real} = (valid(a); TaylorPolynomial{T}(zero(T)))
Base.one(a::TaylorPolynomial{T}) where {T <: Real} = (valid(a); TaylorPolynomial{T}(one(T)))
Base.float(a::TaylorPolynomial) = a
Base.float(::Type{TaylorPolynomial}) = TaylorPolynomial{Float64}
Base.float(::Type{TaylorPolynomial{T}}) where {T <: Real} = TaylorPolynomial{float(T)}
Base.eps(::Type{TaylorPolynomial}) = eps(Float64)
Base.eps(::Type{TaylorPolynomial{T}}) where {T <: Real} = eps(T)
Base.eps(a::TaylorPolynomial) = eps(constant_term(a))
Base.Float64(a::TaylorPolynomial) = Float64(constant_term(a))
Base.Float32(a::TaylorPolynomial) = Float32(constant_term(a))
Base.BigFloat(a::TaylorPolynomial) = BigFloat(constant_term(a))
Base.conj(a::TaylorPolynomial) = a
Base.real(a::TaylorPolynomial) = a
Base.abs2(a::TaylorPolynomial) = a * a
function Base.:+(a::TaylorPolynomial)
    ctx = valid(a)
    return a.len <= ctx.basis.ends[ctx.cutoff + 1] ? a : trim(a, 0, ctx.cutoff)
end
Base.:-(a::TaylorPolynomial) = -one(coefficient_type(a)) * a
Base.signbit(a::TaylorPolynomial) = signbit(constant_term(a))
Base.sign(a::TaylorPolynomial) = sign(constant_term(a))
Base.abs(a::TaylorPolynomial) = signbit(a) ? -a : copy(a)
Base.iszero(a::TaylorPolynomial) = (valid(a); a.len == 1 && coefficient_iszero(a.coeffs[1]))
Base.isnan(a::TaylorPolynomial) = (valid(a); any(isnan, @view a.coeffs[1:a.len]))
Base.isinf(a::TaylorPolynomial) = (valid(a); any(isinf, @view a.coeffs[1:a.len]))
Base.isfinite(a::TaylorPolynomial) = !isnan(a) && !isinf(a)
Base.hash(a::TaylorPolynomial, h::UInt) = hash(constant_term(a), h)
for op in (:(==), :<, :<=, :>, :>=, :isless, :isequal)
    @eval begin
        Base.$op(a::TaylorPolynomial, b::TaylorPolynomial) = Base.$op(constant_term(a), constant_term(b))
        Base.$op(a::TaylorPolynomial, b::Real) = Base.$op(constant_term(a), b)
        Base.$op(a::Real, b::TaylorPolynomial) = Base.$op(a, constant_term(b))
    end
end
for op in (:isless, :isequal, :(==))
    for S in (:AbstractFloat, :AbstractIrrational)
        @eval begin
            Base.$op(a::TaylorPolynomial, b::$S) = Base.$op(constant_term(a), b)
            Base.$op(a::$S, b::TaylorPolynomial) = Base.$op(a, constant_term(b))
        end
    end
end

"Set out = alpha*a + beta*b. Aliasing either input is supported."
@inline function weighted_sum!(out::TaylorPolynomial, a::TaylorPolynomial, alpha::Real, b::TaylorPolynomial, beta::Real; order::Int = a.algebra.cutoff)
    ctx = compatible(a, b)
    valid(out) === ctx || throw(ArgumentError("Different output context"))
    n = min(max(a.len, b.len), ctx.basis.ends[order + 1])
    reserve!(out, n)
    R = coefficient_type(out)
    alpha, beta = coefficient_operand(R, alpha), coefficient_operand(R, beta)
    @inbounds @simd for i in 1:n
        ac = i <= a.len ? a.coeffs[i] : zero(coefficient_type(a))
        bc = i <= b.len ? b.coeffs[i] : zero(coefficient_type(b))
        out.coeffs[i] = alpha * coefficient_operand(R, ac) + beta * coefficient_operand(R, bc)
    end
    return finish!(out, n)
end

function Base.muladd(alpha::Real, a::TaylorPolynomial, b::TaylorPolynomial)
    ctx = compatible(a, b)
    T = promote_type(typeof(alpha), coefficient_type(a), coefficient_type(b))
    n = min(max(a.len, b.len), ctx.basis.ends[ctx.cutoff + 1])
    return weighted_sum!(allocate_undef(ctx, T, n), a, alpha, b, one(T))
end
Base.muladd(a::TaylorPolynomial, alpha::Real, b::TaylorPolynomial) = muladd(alpha, a, b)
Base.muladd(a::TaylorPolynomial, b::TaylorPolynomial, c::TaylorPolynomial) = a * b + c
function add!(out::TaylorPolynomial, a::TaylorPolynomial, b::TaylorPolynomial)
    return weighted_sum!(out, a, one(coefficient_type(a)), b, one(coefficient_type(b)))
end
for (op, sign) in ((:+, 1), (:-, -1))
    @eval function Base.$op(a::TaylorPolynomial{T}, b::TaylorPolynomial{S}) where {T <: Real, S <: Real}
        ctx = compatible(a, b)
        out = allocate_undef(ctx, promote_type(T, S), min(max(a.len, b.len), ctx.basis.ends[ctx.cutoff + 1]))
        return weighted_sum!(out, a, one(T), b, degree_factor(one(S), $sign) * one(S))
    end
end
function scale!(out::TaylorPolynomial, a::TaylorPolynomial, factor::Real; order::Int = a.algebra.cutoff)
    ctx = compatible(out, a)
    n = min(a.len, ctx.basis.ends[order + 1])
    reserve!(out, n)
    R = coefficient_type(out)
    factor = coefficient_operand(R, factor)
    @inbounds @simd for i in 1:n
        out.coeffs[i] = factor * coefficient_operand(R, a.coeffs[i])
    end
    return finish!(out, n)
end
function Base.:*(a::TaylorPolynomial{T}, b::Real) where {T <: Real}
    ctx = valid(a)
    out = allocate_undef(ctx, promote_type(T, typeof(b)), min(a.len, ctx.basis.ends[ctx.cutoff + 1]))
    return scale!(out, a, b)
end
Base.:*(a::Real, b::TaylorPolynomial) = b * a
function Base.:+(a::TaylorPolynomial{T}, b::Real) where {T <: Real}
    ctx = valid(a)
    R = promote_type(T, typeof(b))
    n = min(a.len, ctx.basis.ends[ctx.cutoff + 1])
    out = allocate_undef(ctx, R, n)
    @inbounds for i in 1:n
        out.coeffs[i] = coefficient_convert(R, a.coeffs[i])
    end
    out.coeffs[1] += coefficient_convert(R, b)
    return finish!(out, n)
end
Base.:+(a::Real, b::TaylorPolynomial) = b + a
Base.:-(a::TaylorPolynomial, b::Real) = a + (-b)
Base.:-(a::Real, b::TaylorPolynomial) = (-b) + a
function Base.:/(a::TaylorPolynomial, b::Real)
    coefficient_iszero(b) && throw(TaylorError("Division by zero"))
    # Invert in the coefficient precision. inv(3) would promote Float32 to
    # Float64 and inject a rounded Float64 reciprocal into BigFloat coefficients.
    return a * inv(coefficient_convert(promote_type(coefficient_type(a), typeof(b)), b))
end
Base.:/(a::Real, b::TaylorPolynomial) = quotient(a, b)
Base.:/(a::TaylorPolynomial, b::TaylorPolynomial) = quotient(a, b)

# For fixed i all product indices differ, so the inner loop has no overlapping writes.
function convolve!(out, a, b, basis, cutoff, ::Val{TABLE}, ::Val{SQUARE}) where {TABLE, SQUARE}
    n = min(a.len, basis.ends[cutoff + 1])
    @inbounds for i in 1:n
        ac = a.coeffs[i]
        coefficient_iszero(ac) && continue
        stop = min(b.len, basis.ends[cutoff - basis.degrees[i] + 1])
        first = SQUARE ? i : 1
        offset = basis.offsets[i]
        for j in first:stop
            bc = b.coeffs[j]
            coefficient_iszero(bc) && continue
            k = TABLE ? Int(basis.products[offset + j]) : product_rank(basis, i, j)
            if SQUARE && j != i
                value = coefficient_operand(eltype(out), ac) * coefficient_operand(eltype(out), bc)
                out[k] += value + value
            else
                out[k] = coefficient_muladd(ac, bc, out[k])
            end
        end
    end
    return out
end
function multiply!(out::TaylorPolynomial, a::TaylorPolynomial, b::TaylorPolynomial, cutoff::Int)
    ctx = compatible(a, b)
    valid(out) === ctx || throw(ArgumentError("Different output context"))
    # Scalar products need neither lookup tables nor a full coefficient buffer.
    a.len == 1 && return scale!(out, b, a.coeffs[1]; order = cutoff)
    b.len == 1 && return scale!(out, a, b.coeffs[1]; order = cutoff)
    basis = ctx.basis
    n = basis.ends[min(cutoff, basis.degrees[a.len] + basis.degrees[b.len]) + 1]
    if out === a || out === b || Base.mightalias(out.coeffs, a.coeffs) || Base.mightalias(out.coeffs, b.coeffs)
        tmp = allocate_undef(ctx, coefficient_type(out), n)
        multiply!(tmp, a, b, cutoff)
        out.coeffs, tmp.coeffs = tmp.coeffs, out.coeffs
        out.len = tmp.len
        return out
    end
    reserve!(out, n)
    @inbounds for i in 1:n
        out.coeffs[i] = zero(coefficient_type(out))
    end
    # Keep value parameters literal: Val(runtime_bool) causes dynamic dispatch
    # and boxing on Julia 1.10/1.13 rather than selecting a specialized kernel.
    if isempty(basis.products)
        if a === b
            convolve!(out.coeffs, a, b, basis, cutoff, Val(false), Val(true))
        else
            convolve!(out.coeffs, a, b, basis, cutoff, Val(false), Val(false))
        end
    elseif a === b
        convolve!(out.coeffs, a, b, basis, cutoff, Val(true), Val(true))
    else
        convolve!(out.coeffs, a, b, basis, cutoff, Val(true), Val(false))
    end
    return finish!(out, n)
end
"Multiply into reusable polynomial storage; scalar arithmetic may still allocate."
LinearAlgebra.mul!(out::TaylorPolynomial, a::TaylorPolynomial, b::TaylorPolynomial) = multiply!(out, a, b, a.algebra.cutoff)
function Base.:*(a::TaylorPolynomial{T}, b::TaylorPolynomial{S}) where {T <: Real, S <: Real}
    ctx = compatible(a, b)
    n = ctx.basis.ends[min(ctx.cutoff, ctx.basis.degrees[a.len] + ctx.basis.degrees[b.len]) + 1]
    return multiply!(allocate_undef(ctx, promote_type(T, S), n), a, b, ctx.cutoff)
end
function Base.:^(a::TaylorPolynomial, n::Integer)
    # Preserve the explicit range check in the existing API, avoiding abs(typemin).
    typemin(Int32) < n <= typemax(Int32) || throw(ArgumentError("Exponent out of range"))
    valid(a)
    n == 0 && return one(a)
    n == 1 && return copy(+a)
    n == 2 && return a * a
    n < 0 && return inv(a)^(-n)
    result, power = one(a), a
    while n > 0
        isodd(n) && (result = result * power)
        n >>= 1
        n > 0 && (power = power * power)
    end
    return result
end
"""
    coefficient_norm(p::TaylorPolynomial, type = 0)

Compute a norm over all polynomial coefficients.
Type zero returns the maximum absolute coefficient; type one returns their
absolute sum; higher integer types use the corresponding power norm.
Unlike scalar polynomial comparisons, this includes nonconstant coefficients.
"""
function coefficient_norm(a::TaylorPolynomial, p::Real = 0)
    valid(a)
    isinteger(p) && p >= 0 || throw(ArgumentError("Use a nonnegative integer coefficient norm"))
    c = @view a.coeffs[1:a.len]
    p == 0 && return maximum(abs, c)
    p == 1 && return sum(abs, c)
    return sum(x -> abs(x)^p, c)^(one(float(constant_term(a))) / p)
end
