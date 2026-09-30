function Base.copy(a::DA{T}) where {T<:Real}
    ctx = valid(a)
    DA{T}(a.coeffs[1:a.len],a.len,ctx,true)
end
function DA{T}(a::DA) where {T<:Real}
    ctx = valid(a)
    DA{T}(T.(a.coeffs[1:a.len]),a.len,ctx,true)
end
Base.convert(::Type{DA},a::DA) = a
Base.convert(::Type{DA},x::Real) = DA(x)
Base.convert(::Type{DA{T}},a::DA{T}) where {T<:Real} = a
Base.convert(::Type{DA{T}},a::DA) where {T<:Real} = DA{T}(a)
Base.convert(::Type{DA{T}},x::Real) where {T<:Real} = DA{T}(x)
Base.promote_rule(::Type{DA{T}},::Type{DA{S}}) where {T<:Real,S<:Real} = DA{promote_type(T,S)}
Base.promote_rule(::Type{DA{T}},::Type{S}) where {T<:Real,S<:Real} = DA{promote_type(T,S)}
Base.promote_rule(::Type{DA},::Type{S}) where {S<:Real} = DA
# BigFloat otherwise promotes every Real to BigFloat, conflicting with DA's
# rule and recursively asking Julia to promote DA{BigFloat} with BigFloat.
Base.promote_rule(::Type{BigFloat},::Type{DA{T}}) where {T<:Real} = DA{promote_type(BigFloat,T)}
Base.promote_rule(::Type{BigFloat},::Type{DA}) = DA
Base.promote_rule(::Type{DA{T}},::Type{DA}) where {T<:Real} = DA
Base.zero(::Type{DA}) = DA()
Base.one(::Type{DA}) = DA(1.0)
Base.zero(::Type{DA{T}}) where {T<:Real} = DA{T}(zero(T))
Base.one(::Type{DA{T}}) where {T<:Real} = DA{T}(one(T))
Base.zero(a::DA{T}) where {T<:Real} = (valid(a); DA{T}(zero(T)))
Base.one(a::DA{T}) where {T<:Real} = (valid(a); DA{T}(one(T)))
Base.float(a::DA) = a
Base.float(::Type{DA}) = DA{Float64}
Base.float(::Type{DA{T}}) where {T<:Real} = DA{float(T)}
Base.eps(::Type{DA}) = eps(Float64)
Base.eps(::Type{DA{T}}) where {T<:Real} = eps(T)
Base.eps(a::DA) = eps(cons(a))
Base.Float64(a::DA) = Float64(cons(a))
Base.Float32(a::DA) = Float32(cons(a))
Base.BigFloat(a::DA) = BigFloat(cons(a))
Base.conj(a::DA) = a
Base.real(a::DA) = a
Base.abs2(a::DA) = a*a
Base.:+(a::DA) = a
Base.:-(a::DA) = -one(coefftype(a))*a
Base.signbit(a::DA) = signbit(cons(a))
Base.sign(a::DA) = sign(cons(a))
Base.abs(a::DA) = signbit(a) ? -a : copy(a)
Base.iszero(a::DA) = (valid(a); a.len == 1 && iszero(a.coeffs[1]))
Base.isnan(a::DA) = (valid(a); any(isnan,@view a.coeffs[1:a.len]))
Base.isinf(a::DA) = (valid(a); any(isinf,@view a.coeffs[1:a.len]))
Base.isfinite(a::DA) = !isnan(a) && !isinf(a)
Base.hash(a::DA,h::UInt) = hash(cons(a),h)
for op in (:(==), :<, :<=, :>, :>=, :isless, :isequal)
    @eval begin
        Base.$op(a::DA,b::DA) = Base.$op(cons(a),cons(b))
        Base.$op(a::DA,b::Real) = Base.$op(cons(a),b)
        Base.$op(a::Real,b::DA) = Base.$op(a,cons(b))
    end
end
for op in (:isless,:isequal,:(==))
    for S in (:AbstractFloat,:AbstractIrrational)
        @eval begin
            Base.$op(a::DA,b::$S) = Base.$op(cons(a),b)
            Base.$op(a::$S,b::DA) = Base.$op(a,cons(b))
        end
    end
end

"Set out = alpha*a + beta*b. Aliasing either input is supported."
@inline function weighted_sum!(out::DA,a::DA,alpha::Real,b::DA,beta::Real; order::Int=a.context.cutoff)
    ctx = compatible(a,b)
    valid(out) === ctx || throw(ArgumentError("Different output context"))
    n = min(max(a.len,b.len),ctx.basis.ends[order+1])
    reserve!(out,n)
    @inbounds for i in 1:n
        ac = i <= a.len ? a.coeffs[i] : zero(coefftype(a))
        bc = i <= b.len ? b.coeffs[i] : zero(coefftype(b))
        out.coeffs[i] = alpha*ac + beta*bc
    end
    finish!(out,n)
end

function Base.muladd(alpha::Real,a::DA,b::DA)
    ctx = compatible(a,b)
    T = promote_type(typeof(alpha),coefftype(a),coefftype(b))
    n = min(max(a.len,b.len),ctx.basis.ends[ctx.cutoff+1])
    weighted_sum!(allocate_undef(ctx,T,n),a,alpha,b,one(T))
end
Base.muladd(a::DA,alpha::Real,b::DA) = muladd(alpha,a,b)
Base.muladd(a::DA,b::DA,c::DA) = a*b+c
function add!(out::DA,a::DA,b::DA)
    weighted_sum!(out,a,one(coefftype(a)),b,one(coefftype(b)))
end
for (op,sign) in ((:+,1),(:-,-1))
    @eval function Base.$op(a::DA{T},b::DA{S}) where {T<:Real,S<:Real}
        ctx = compatible(a,b)
        out = allocate_undef(ctx,promote_type(T,S),min(max(a.len,b.len),ctx.basis.ends[ctx.cutoff+1]))
        weighted_sum!(out,a,one(T),b,$sign*one(S))
    end
end
function scale!(out::DA,a::DA,factor::Real; order::Int=a.context.cutoff)
    ctx = compatible(out,a)
    n = min(a.len,ctx.basis.ends[order+1])
    reserve!(out,n)
    @inbounds @simd for i in 1:n
        out.coeffs[i] = factor*a.coeffs[i]
    end
    finish!(out,n)
end
function Base.:*(a::DA{T},b::Real) where {T<:Real}
    ctx = valid(a)
    out = allocate_undef(ctx,promote_type(T,typeof(b)),min(a.len,ctx.basis.ends[ctx.cutoff+1]))
    scale!(out,a,b)
end
Base.:*(a::Real,b::DA) = b*a
function Base.:+(a::DA{T},b::Real) where {T<:Real}
    ctx = valid(a)
    out = DA{promote_type(T,typeof(b))}(a)
    out.coeffs[1] += b
    finish!(out)
end
Base.:+(a::Real,b::DA) = b+a
Base.:-(a::DA,b::Real) = a+(-b)
Base.:-(a::Real,b::DA) = (-b)+a
function Base.:/(a::DA,b::Real)
    iszero(b) && throw(DAError("Division by zero"))
    # Invert in the coefficient precision. inv(3) would promote Float32 to
    # Float64 and inject a rounded Float64 reciprocal into BigFloat coefficients.
    a*inv(convert(promote_type(coefftype(a),typeof(b)),b))
end
Base.:/(a::Real,b::DA) = quotient(a,b)
Base.:/(a::DA,b::DA) = quotient(a,b)

# For fixed i all product indices differ, so the inner loop has no overlapping writes.
function convolve!(out,a,b,basis,cutoff,::Val{TABLE},::Val{SQUARE}) where {TABLE,SQUARE}
    n = min(a.len,basis.ends[cutoff+1])
    @inbounds for i in 1:n
        ac = a.coeffs[i]
        iszero(ac) && continue
        stop = min(b.len,basis.ends[cutoff-basis.degrees[i]+1])
        first = SQUARE ? i : 1
        offset = basis.offsets[i]
        for j in first:stop
            bc = b.coeffs[j]
            iszero(bc) && continue
            k = TABLE ? Int(basis.products[offset+j]) : product_rank(basis,i,j)
            if SQUARE && j != i
                value = ac*bc
                out[k] += value+value
            else
                out[k] = muladd(ac,bc,out[k])
            end
        end
    end
    out
end
function multiply!(out::DA,a::DA,b::DA,cutoff::Int)
    ctx = compatible(a,b)
    valid(out) === ctx || throw(ArgumentError("Different output context"))
    # Scalar products need neither lookup tables nor a full coefficient buffer.
    a.len == 1 && return scale!(out,b,a.coeffs[1]; order=cutoff)
    b.len == 1 && return scale!(out,a,b.coeffs[1]; order=cutoff)
    basis = ctx.basis
    n = basis.ends[min(cutoff,basis.degrees[a.len]+basis.degrees[b.len])+1]
    if out === a || out === b || Base.mightalias(out.coeffs,a.coeffs) || Base.mightalias(out.coeffs,b.coeffs)
        tmp = allocate_undef(ctx,coefftype(out),n)
        multiply!(tmp,a,b,cutoff)
        out.coeffs,tmp.coeffs = tmp.coeffs,out.coeffs
        out.len = tmp.len
        return out
    end
    reserve!(out,n)
    @inbounds for i in 1:n
        out.coeffs[i] = zero(coefftype(out))
    end
    # Keep value parameters literal: Val(runtime_bool) causes dynamic dispatch
    # and boxing on Julia 1.10/1.13 rather than selecting a specialized kernel.
    if isempty(basis.products)
        if a === b
            convolve!(out.coeffs,a,b,basis,cutoff,Val(false),Val(true))
        else
            convolve!(out.coeffs,a,b,basis,cutoff,Val(false),Val(false))
        end
    elseif a === b
        convolve!(out.coeffs,a,b,basis,cutoff,Val(true),Val(true))
    else
        convolve!(out.coeffs,a,b,basis,cutoff,Val(true),Val(false))
    end
    finish!(out,n)
end
"Multiply into a reusable polynomial; disjoint Float32/Float64 buffers allocate nothing after sizing."
LinearAlgebra.mul!(out::DA,a::DA,b::DA) = multiply!(out,a,b,a.context.cutoff)
function Base.:*(a::DA{T},b::DA{S}) where {T<:Real,S<:Real}
    ctx = compatible(a,b)
    n = ctx.basis.ends[min(ctx.cutoff,ctx.basis.degrees[a.len]+ctx.basis.degrees[b.len])+1]
    multiply!(allocate_undef(ctx,promote_type(T,S),n),a,b,ctx.cutoff)
end
function Base.:^(a::DA,n::Integer)
    # Preserve the explicit range check in the existing API, avoiding abs(typemin).
    typemin(Int32) < n <= typemax(Int32) || throw(ArgumentError("Exponent out of range"))
    valid(a)
    n == 0 && return one(a)
    n == 1 && return copy(a)
    n == 2 && return a*a
    n < 0 && return inv(a)^(-n)
    result, power = one(a), a
    while n > 0
        isodd(n) && (result = result*power)
        n >>= 1
        n > 0 && (power = power*power)
    end
    result
end
powi(a::DA,n::Integer) = a^n
sqr(a::DA) = a*a
sqr(a::Real) = a*a
function norm(a::DA,p::Real=0)
    valid(a)
    isinteger(p) && p >= 0 || throw(ArgumentError("Use a nonnegative integer coefficient norm"))
    c = @view a.coeffs[1:a.len]
    p == 0 && return maximum(abs,c)
    p == 1 && return sum(abs,c)
    sum(x -> abs(x)^p,c)^(one(float(cons(a)))/p)
end
