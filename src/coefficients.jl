struct Monomial{T <: Real}
    coefficient::T
    exponents::Vector{UInt32}
end
Monomial{T}() where {T <: Real} = Monomial{T}(zero(T), zeros(UInt32, getMaxVariables()))
Monomial() = Monomial{Float64}()
getCoefficient(m::Monomial) = m.coefficient
getExponents(m::Monomial) = copy(m.exponents)
order(m::Monomial) = sum(m.exponents)
function exponent_vector(jj::AbstractVector{<:Integer}, basis::Basis = ready().basis)
    all(j -> 0 <= j <= typemax(UInt32), jj) || throw(ArgumentError("Invalid monomial exponents"))
    out = zeros(Int, basis.variables)
    for (i, j) in enumerate(jj)
        i > length(out) && break
        out[i] = j
    end
    return out
end
function getCoefficient(a::DA{T}, jj::AbstractVector{<:Integer}) where {T}
    b = valid(a).basis
    exponents = exponent_vector(jj, b)
    degree = sum(exponents)
    degree > b.order && return zero(T)
    i = rank(b, exponents, degree)
    return i <= a.len ? a.coeffs[i] : zero(T)
end
function coefficient(a::DA, jj::AbstractVector{<:Integer})
    b = valid(a).basis
    length(jj) == b.variables || throw(DimensionMismatch("One exponent per variable is required"))
    all(>=(0), jj) && sum(big, jj) <= b.order || throw(ArgumentError("Invalid exponents"))
    return getCoefficient(a, jj)
end
function setCoefficient!(a::DA{T}, jj::AbstractVector{<:Integer}, value::Real) where {T}
    ctx = valid(a)
    exponents = exponent_vector(jj, ctx.basis)
    degree = sum(exponents)
    degree <= ctx.basis.order || throw(ArgumentError("Monomial exceeds initialized order"))
    i = rank(ctx.basis, exponents, degree)
    reserve!(a, i)
    @inbounds for j in (a.len + 1):i
        a.coeffs[j] = zero(T)
    end
    c = convert(T, value)
    a.coeffs[i] = keep(c, ctx) ? c : zero(T)
    return finish!(a, max(i, a.len))
end
function getMonomials(a::DA)
    b = valid(a).basis
    return [Monomial(a.coeffs[i], UInt32.(b.exponents[:, i])) for i in 1:a.len if !iszero(a.coeffs[i])]
end
function getMonomial(a::DA, pos::Integer)
    pos >= 1 || throw(BoundsError(a, pos))
    b = valid(a).basis
    found = 0
    for i in 1:a.len
        if !iszero(a.coeffs[i])
            found += 1
            found == pos && return Monomial(a.coeffs[i], UInt32.(b.exponents[:, i]))
        end
    end
    throw(BoundsError(a, pos))
end
function trim(a::DA, low::Integer, high::Integer = getMaxOrder())
    ctx = valid(a)
    low >= 0 && high >= 0 || throw(ArgumentError("Invalid order range"))
    low > high && return zero(a)
    n = min(a.len, ctx.basis.ends[min(high, ctx.basis.order) + 1])
    out = allocate(ctx, coefftype(a), n)
    @inbounds for i in 1:n
        ctx.basis.degrees[i] >= low && (out.coeffs[i] = a.coeffs[i])
    end
    return finish!(out, n)
end
function multiplyMonomials(a::DA, b::DA)
    ctx = compatible(a, b)
    n = min(a.len, b.len, ctx.basis.ends[ctx.cutoff + 1])
    out = allocate(ctx, promote_type(coefftype(a), coefftype(b)), n)
    @inbounds for i in 1:n
        out.coeffs[i] = a.coeffs[i] * b.coeffs[i]
    end
    return finish!(out, n)
end
function calculus(a::DA, v::Integer, operation::Symbol, p::Integer = 1)
    ctx = valid(a); b = ctx.basis
    1 <= v <= b.variables && 0 <= p <= b.order || throw(ArgumentError("Invalid variable power"))
    T = operation === :integrate ? typeof(one(coefftype(a)) / 1) : coefftype(a)
    out = allocate(ctx, T)
    exponents = zeros(Int, b.variables)
    for i in 1:a.len
        c = a.coeffs[i]
        iszero(c) && continue
        exponents .= @view b.exponents[:, i]
        degree = b.degrees[i]
        if operation === :integrate
            degree >= ctx.cutoff && continue
            exponents[v] += 1
            out.coeffs[rank(b, exponents, degree + 1)] = c / exponents[v]
        elseif operation === :deriv
            exponents[v] == 0 && continue
            power = exponents[v]; exponents[v] -= 1
            degree - 1 <= ctx.cutoff && (out.coeffs[rank(b, exponents, degree - 1)] = power * c)
        else
            exponents[v] < p && throw(DAError("Polynomial is not divisible by this variable power"))
            exponents[v] -= p
            degree - p <= ctx.cutoff && (out.coeffs[rank(b, exponents, degree - p)] = c)
        end
    end
    return finish!(out)
end
deriv(a::DA, i::Integer) = calculus(a, i, :deriv)
integrate(a::DA, i::Integer) = calculus(a, i, :integrate)
divide(a::DA, i::Integer, p::Integer = 1) = calculus(a, i, :divide, p)
for fn in (:deriv, :integrate)
    @eval function $fn(a::DA, counts::AbstractVector{<:Integer})
        all(>=(0), counts) || throw(ArgumentError("Counts must be nonnegative"))
        out = copy(a)
        for (i, n) in enumerate(counts)
            i > a.context.basis.variables && break
            for _ in 1:n
                out = $fn(out, i)
            end
        end
        return out
    end
end
const integ = integrate
function orderNorm(a::DA, v::Integer = 0, p::Integer = 0)
    b = valid(a).basis
    0 <= v <= b.variables && p >= 0 || throw(ArgumentError("Invalid norm arguments"))
    out = zeros(typeof(abs(cons(a))), b.order + 1)
    @inbounds for i in 1:a.len
        degree = v == 0 ? b.degrees[i] : b.exponents[v, i]
        c = abs(a.coeffs[i])
        out[degree + 1] = p == 0 ? max(out[degree + 1], c) : out[degree + 1] + c^p
    end
    p > 1 && (out .= out .^ (one(eltype(out)) / p))
    return out
end
function estimNorm(a::DA, v::Integer = 0, p::Integer = 0, n::Integer = getMaxOrder(); errors::Bool = false)
    n >= 0 || throw(ArgumentError("Invalid estimated order"))
    norms = orderNorm(a, v, p)
    length(norms) >= 3 || throw(DAError("At least second order is required to estimate norms"))
    indices = findall(i -> keep(norms[i + 1], a.context), 1:(length(norms) - 1))
    if length(indices) < 2
        @warn "Insufficient nonzero orders to estimate coefficient norms"
        estimates = zeros(eltype(norms), n + 1)
        return errors ? (estimates, zeros(eltype(norms), min(n, length(norms) - 1) + 1)) : estimates
    end
    x = indices; y = log.(norms[indices .+ 1])
    count = length(x); sx = sum(x); sy = sum(y)
    slope = (count * sum(x .* y) - sx * sy) / (count * sum(abs2, x) - sx * sx)
    intercept = (sy - slope * sx) / count
    estimates = [exp(intercept + slope * i) for i in 0:n]
    errors || return estimates
    count = min(length(norms), n + 1)
    return estimates, max.(zero(eltype(norms)), norms[1:count] - estimates[1:count])
end
struct Interval{T <: Real}
    m_lb::T
    m_ub::T
end
function bound(a::DA)
    b = valid(a).basis
    lo = hi = cons(a)
    for i in 2:a.len
        c = a.coeffs[i]
        if all(iseven, @view b.exponents[:, i])
            lo += min(c, zero(c)); hi += max(c, zero(c))
        else
            lo -= abs(c); hi += abs(c)
        end
    end
    return Interval(lo, hi)
end
function random(filling::Real = -1.0, ::Type{T} = Float64) where {T <: AbstractFloat}
    ctx = ready(); a = allocate(ctx, T)
    density = min(abs(filling), 1)
    for i in eachindex(a.coeffs)
        if Random.rand() < density
            a.coeffs[i] = 2Random.rand(T) - one(T)
            filling < 0 && (a.coeffs[i] = ldexp(a.coeffs[i], -ctx.basis.degrees[i]))
        end
    end
    return finish!(a)
end
# Map calculus preserves the shape of ordinary arrays, including views.
trim(a::AbstractArray{<:DA}, low::Integer, high::Integer = getMaxOrder()) = trim.(a, low, high)
deriv(a::AbstractArray{<:DA}, i::Integer) = deriv.(a, i)
integrate(a::AbstractArray{<:DA}, i::Integer) = integrate.(a, i)
