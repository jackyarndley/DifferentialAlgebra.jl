function Base.copy(a::TaylorPolynomial{T}) where {T <: Real}
    ctx = valid(a)
    return TaylorPolynomial{T}(a.coeffs[1:a.len], a.len, ctx)
end
function TaylorPolynomial{T}(a::TaylorPolynomial) where {T <: Real}
    ctx = valid(a)
    out = TaylorPolynomial{T}(
        [coefficient_convert(T, a.coeffs[i]) for i in 1:a.len], a.len, ctx
    )
    return finish!(out, a.len; filter = false)
end
Base.convert(::Type{TaylorPolynomial}, a::TaylorPolynomial) = copy(a)
Base.convert(::Type{TaylorPolynomial}, x::Real) = TaylorPolynomial(x)
Base.convert(::Type{TaylorPolynomial{T}}, a::TaylorPolynomial{T}) where {T <: Real} = copy(a)
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
Base.zero(a::TaylorPolynomial{T}) where {T <: Real} =
    TaylorPolynomial{T}(valid(a), zero(T))
Base.one(a::TaylorPolynomial{T}) where {T <: Real} =
    TaylorPolynomial{T}(valid(a), one(T))
Base.float(a::TaylorPolynomial) = convert(float(typeof(a)), a)
Base.float(::Type{TaylorPolynomial}) = TaylorPolynomial{Float64}
Base.float(::Type{TaylorPolynomial{T}}) where {T <: Real} = TaylorPolynomial{float(T)}
Base.eps(::Type{TaylorPolynomial}) = eps(Float64)
Base.eps(::Type{TaylorPolynomial{T}}) where {T <: Real} = eps(T)
Base.eps(a::TaylorPolynomial) = eps(constant_term(a))
function scalar_value(::Type{T}, a::TaylorPolynomial) where {T}
    valid(a)
    a.len == 1 || throw(InexactError(:convert, T, a))
    return T(a.coeffs[1])
end
Base.Float64(a::TaylorPolynomial) = scalar_value(Float64, a)
Base.Float32(a::TaylorPolynomial) = scalar_value(Float32, a)
Base.BigFloat(a::TaylorPolynomial) = scalar_value(BigFloat, a)
Base.Bool(a::TaylorPolynomial) = scalar_value(Bool, a)
Base.conj(a::TaylorPolynomial) = copy(a)
Base.real(a::TaylorPolynomial) = copy(a)
Base.abs2(a::TaylorPolynomial) = a * a
function Base.:+(a::TaylorPolynomial)
    ctx = valid(a)
    return a.len <= ctx.basis.ends[working_order(ctx) + 1] ?
        copy(a) : trim(a, 0, working_order(ctx))
end
negative_coefficient_type(::Type{T}) where {T} = T
negative_coefficient_type(::Type{Bool}) = Int
negative_coefficient_type(::Type{T}) where {T <: Unsigned} =
    widen(typeof(signed(zero(T))))
negative_coefficient_type(::Type{Rational{T}}) where {T <: Unsigned} =
    Rational{negative_coefficient_type(T)}
Base.:-(a::TaylorPolynomial) =
    scaled(a, -one(negative_coefficient_type(coefficient_type(a))))
Base.signbit(a::TaylorPolynomial) = signbit(constant_term(a))
Base.sign(a::TaylorPolynomial) = sign(constant_term(a))
Base.abs(a::TaylorPolynomial) = a * (signbit(a) ? -one(coefficient_type(a)) : one(coefficient_type(a)))
Base.iszero(a::TaylorPolynomial) = (valid(a); a.len == 1 && coefficient_iszero(a.coeffs[1]))
Base.isone(a::TaylorPolynomial) = (valid(a); a.len == 1 && coefficient_isone(a.coeffs[1]))
Base.isinteger(a::TaylorPolynomial) =
    (valid(a); a.len == 1 && coefficient_isinteger(a.coeffs[1]))
Base.isnan(a::TaylorPolynomial) = (valid(a); any(isnan, @view a.coeffs[1:a.len]))
Base.isinf(a::TaylorPolynomial) = (valid(a); any(isinf, @view a.coeffs[1:a.len]))
Base.isfinite(a::TaylorPolynomial) = !isnan(a) && !isinf(a)
# Structural identity remains available after reinitialization, so dictionary
# keys do not become unhashable when numerical use of their algebra expires.
function Base.hash(a::TaylorPolynomial, h::UInt)
    a.len == 1 && return hash(a.coeffs[1], h)
    return hash((a.algebra.basis.variables, @view(a.coeffs[1:a.len])), h)
end
function Base.isequal(a::TaylorPolynomial, b::TaylorPolynomial)
    a.len == b.len || return false
    a.len > 1 && a.algebra.basis.variables != b.algebra.basis.variables && return false
    return all(i -> coefficient_isequal(a.coeffs[i], b.coeffs[i]), 1:a.len)
end
Base.isequal(a::TaylorPolynomial, b::Real) =
    a.len == 1 && coefficient_isequal(a.coeffs[1], b)
Base.isequal(a::Real, b::TaylorPolynomial) = isequal(b, a)
for S in (:AbstractFloat, :AbstractIrrational)
    @eval begin
        Base.isequal(a::TaylorPolynomial, b::$S) =
            a.len == 1 && coefficient_isequal(a.coeffs[1], b)
        Base.isequal(a::$S, b::TaylorPolynomial) = isequal(b, a)
    end
end
# Base defines > and >= by reversing < and <=.
for op in (:(==), :<, :<=, :isless)
    @eval begin
        Base.$op(a::TaylorPolynomial, b::TaylorPolynomial) = Base.$op(constant_term(a), constant_term(b))
        Base.$op(a::TaylorPolynomial, b::Real) = Base.$op(constant_term(a), b)
        Base.$op(a::Real, b::TaylorPolynomial) = Base.$op(a, constant_term(b))
    end
end
for op in (:isless, :(==))
    for S in (:AbstractFloat, :AbstractIrrational)
        @eval begin
            Base.$op(a::TaylorPolynomial, b::$S) = Base.$op(constant_term(a), b)
            Base.$op(a::$S, b::TaylorPolynomial) = Base.$op(a, constant_term(b))
        end
    end
end

"Set out = alpha*a + beta*b. Aliasing either input is supported."
@inline function weighted_sum!(
        out::TaylorPolynomial, a::TaylorPolynomial, alpha::Real,
        b::TaylorPolynomial, beta::Real;
        order::Int = working_order(a.algebra), filter::Bool = true
    )
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
    return finish!(out, n; filter)
end

function Base.muladd(alpha::Real, a::TaylorPolynomial, b::TaylorPolynomial)
    ctx = compatible(a, b)
    T = promote_type(typeof(alpha), coefficient_type(a), coefficient_type(b))
    order = working_order(ctx)
    n = min(max(a.len, b.len), ctx.basis.ends[order + 1])
    return weighted_sum!(allocate_undef(ctx, T, n), a, alpha, b, one(T); order)
end
Base.muladd(a::TaylorPolynomial, alpha::Real, b::TaylorPolynomial) = muladd(alpha, a, b)
Base.muladd(a::TaylorPolynomial, b::TaylorPolynomial, c::TaylorPolynomial) =
    product(a, b; filter = false) + c
"""
    DifferentialAlgebra.add!(out, p, q)

Internal addition buffer kernel. Store `p + q` in `out`, applying its
coefficient type and the working order
and tolerance. All three polynomials must share an algebra. Either input may
alias `out`; operand precision is retained until assignment to its coefficients.
"""
function add!(out::TaylorPolynomial, a::TaylorPolynomial, b::TaylorPolynomial)
    return weighted_sum!(out, a, one(coefficient_type(a)), b, one(coefficient_type(b)))
end
for (op, sign) in ((:+, 1), (:-, -1))
    @eval function Base.$op(a::TaylorPolynomial{T}, b::TaylorPolynomial{S}) where {T <: Real, S <: Real}
        ctx = compatible(a, b)
        R = promote_type(T, S)
        $sign < 0 && (R = negative_coefficient_type(R))
        order = working_order(ctx)
        n = min(max(a.len, b.len), ctx.basis.ends[order + 1])
        out = allocate_undef(ctx, R, n)
        return weighted_sum!(
            out, a, one(R), b, degree_factor(one(R), $sign) * one(R); order
        )
    end
end
function scale!(
        out::TaylorPolynomial, a::TaylorPolynomial, factor::Real;
        order::Int = working_order(a.algebra), filter::Bool = true
    )
    ctx = compatible(out, a)
    n = min(a.len, ctx.basis.ends[order + 1])
    reserve!(out, n)
    R = coefficient_type(out)
    factor = coefficient_operand(R, factor)
    @inbounds @simd for i in 1:n
        out.coeffs[i] = factor * coefficient_operand(R, a.coeffs[i])
    end
    return finish!(out, n; filter)
end
function scaled(a::TaylorPolynomial{T}, b::Real; filter::Bool = true) where {T <: Real}
    ctx = valid(a)
    order = working_order(ctx)
    out = allocate_undef(
        ctx, promote_type(T, typeof(b)),
        min(a.len, ctx.basis.ends[order + 1])
    )
    return scale!(out, a, b; order, filter)
end
Base.:*(a::TaylorPolynomial, b::Real) = scaled(a, b)
Base.:*(a::Real, b::TaylorPolynomial) = b * a
function scalar_sum(
        a::TaylorPolynomial{T}, b::Real, ::Val{OP}
    ) where {T <: Real, OP}
    ctx = valid(a)
    R = promote_type(T, typeof(b))
    OP !== :add && (R = negative_coefficient_type(R))
    n = min(a.len, ctx.basis.ends[working_order(ctx) + 1])
    out = allocate_undef(ctx, R, n)
    @inbounds for i in 1:n
        c = coefficient_convert(R, a.coeffs[i])
        out.coeffs[i] = OP === :subtract_from ? -c : c
    end
    c = coefficient_convert(R, b)
    out.coeffs[1] = OP === :subtract ? out.coeffs[1] - c : out.coeffs[1] + c
    return finish!(out, n)
end
Base.:+(a::TaylorPolynomial, b::Real) = scalar_sum(a, b, Val(:add))
Base.:+(a::Real, b::TaylorPolynomial) = b + a
Base.:-(a::TaylorPolynomial, b::Real) = scalar_sum(a, b, Val(:subtract))
Base.:-(a::Real, b::TaylorPolynomial) = scalar_sum(b, a, Val(:subtract_from))
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
function multiply!(out::TaylorPolynomial, a::TaylorPolynomial, b::TaylorPolynomial, cutoff::Int; filter::Bool = true)
    ctx = compatible(a, b)
    valid(out) === ctx || throw(ArgumentError("Different output context"))
    # Scalar products need neither lookup tables nor a full coefficient buffer.
    a.len == 1 && return scale!(out, b, a.coeffs[1]; order = cutoff, filter)
    b.len == 1 && return scale!(out, a, b.coeffs[1]; order = cutoff, filter)
    basis = ctx.basis
    n = basis.ends[min(cutoff, basis.degrees[a.len] + basis.degrees[b.len]) + 1]
    if out === a || out === b || Base.mightalias(out.coeffs, a.coeffs) || Base.mightalias(out.coeffs, b.coeffs)
        tmp = allocate_undef(ctx, coefficient_type(out), n)
        multiply!(tmp, a, b, cutoff; filter)
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
    return finish!(out, n; filter)
end
"Multiply into reusable polynomial storage; scalar arithmetic may still allocate."
LinearAlgebra.mul!(out::TaylorPolynomial, a::TaylorPolynomial, b::TaylorPolynomial) =
    multiply!(out, a, b, working_order(a.algebra))
function product(
        a::TaylorPolynomial{T}, b::TaylorPolynomial{S}; filter::Bool = true
    ) where {T <: Real, S <: Real}
    ctx = compatible(a, b)
    order = min(working_order(ctx), ctx.basis.degrees[a.len] + ctx.basis.degrees[b.len])
    n = ctx.basis.ends[order + 1]
    return multiply!(
        allocate_undef(ctx, promote_type(T, S), n), a, b, working_order(ctx); filter
    )
end
Base.:*(a::TaylorPolynomial, b::TaylorPolynomial) = product(a, b)
# Base's negative literal-power fallback filters inv(a) before raising it again.
Base.literal_pow(::typeof(^), a::TaylorPolynomial, ::Val{N}) where {N} = a^N
# Nonnegative powers keep one concrete accumulator type. Local filtering also
# lets negative powers retain reciprocal coefficients until the final product.
function integer_power(a::TaylorPolynomial, n::Integer; filter::Bool = true)
    ctx = valid(a)
    if a.len == 1
        c0 = a.coeffs[1]^n
        return finish!(TaylorPolynomial{typeof(c0)}([c0], 1, ctx), 1; filter)
    end
    n == 0 && return finish!(
        TaylorPolynomial{coefficient_type(a)}([one(coefficient_type(a))], 1, ctx), 1;
        filter
    )
    n > working_order(ctx) && coefficient_iszero(a.coeffs[1]) && return zero(a)
    n == 1 && return scaled(a, one(coefficient_type(a)); filter)
    n == 2 && return product(a, a; filter)
    while iseven(n)
        a = product(a, a; filter = false)
        n >>= 1
    end
    result = scaled(a, one(coefficient_type(a)); filter = false)
    n >>= 1
    while n > 0
        a = product(a, a; filter = false)
        isodd(n) && (result = product(result, a; filter = false))
        n >>= 1
    end
    return finish!(result, result.len; filter)
end
function Base.:^(a::TaylorPolynomial, n::Integer)
    n >= 0 && return integer_power(a, n)
    n == -1 && return quotient(one(coefficient_type(a)), a)
    reciprocal = quotient(one(coefficient_type(a)), a; filter = false)
    # -(n+1) is representable even for typemin of a signed integer type.
    return product(integer_power(reciprocal, -(n + 1); filter = false), reciprocal)
end

"""
    coefficient_norm(p::TaylorPolynomial, type = 0)

Compute a norm over all polynomial coefficients.
Type zero returns the maximum absolute coefficient; type one returns their
absolute sum; higher integer types use the corresponding power norm.
Unlike scalar polynomial comparisons, this includes nonconstant coefficients.
The result uses the floating norm type of the coefficients, as in `LinearAlgebra.norm`.
Interval coefficients produce an outward-rounded enclosure of the possible norms.
They require `IntervalArithmetic.configure(rounding=:correct)`.
"""
function coefficient_norm(a::TaylorPolynomial, p::Real = 0)
    valid(a)
    isinteger(p) && p >= 0 || throw(ArgumentError("Use a nonnegative integer coefficient norm"))
    c = @view a.coeffs[1:a.len]
    return LinearAlgebra.norm(c, p == 0 ? Inf : p)
end
