# Shared coefficient recurrences for eager TimeSeries arithmetic and recorded ODEs.
# Coefficients use zero-based mathematical degrees: degree k lives at index k + 1.
# Callers fill earlier degrees before requesting k; exact and polynomial arithmetic
# retain separate kernels where their accumulation and storage rules differ.

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
function time_convolution(a::Vector{P}, b::Vector{P}, k, prototype::P) where {P <: TaylorPolynomial}
    value, scratch = zero(prototype), zero(prototype)
    @inbounds for j in 0:k
        (iszero(a[j + 1]) || iszero(b[k - j + 1])) && continue
        LinearAlgebra.mul!(scratch, a[j + 1], b[k - j + 1])
        add!(value, value, scratch)
    end
    return value
end

@inline time_weight(j, k, p, ::Val{:exp}) = j
@inline time_weight(j, k, p, ::Val{:log}) = j - k
@inline time_weight(j, k, p, ::Val{:power}) = p * j - (k - j)

# Normalize scalar exponent weights once, keeping exact coefficient arithmetic
# and the wider working type used by BigInt rationals. Polynomial and dual
# coefficients only supply the scalar type; their values remain untouched.
time_exponent(p, value) = p
time_exponent(p::Rational, value::AbstractFloat) =
    convert(promote_type(typeof(p), typeof(value)), p)
time_exponent(p::Rational, value::TaylorPolynomial) =
    time_exponent(p, constant_term(value))

function time_weighted_sum(a, c, k, p, mode)
    value = zero(c[1])
    @inbounds for j in 1:k
        weight = time_weight(j, k, p, mode)
        iszero(weight) && continue
        value += weight * a[j + 1] * c[k - j + 1]
    end
    return value
end
function time_weighted_sum(a::Vector{T}, c::Vector{T}, k, p, mode) where {T <: Union{Float32, Float64}}
    value = zero(T)
    @inbounds @simd for j in 1:k
        value += time_weight(j, k, p, mode) * a[j + 1] * c[k - j + 1]
    end
    return value
end
function time_weighted_sum(a::Vector{P}, c::Vector{P}, k, p, mode) where {P <: TaylorPolynomial}
    value, scratch = zero(c[1]), zero(c[1])
    unit = one(coefficient_type(value))
    @inbounds for j in 1:k
        weight = time_weight(j, k, p, mode)
        (iszero(weight) || iszero(a[j + 1]) || iszero(c[k - j + 1])) && continue
        LinearAlgebra.mul!(scratch, a[j + 1], c[k - j + 1])
        weighted_sum!(value, value, unit, scratch, weight)
    end
    return value
end
time_normalize(value, reciprocal, k) = (value / k) * reciprocal
function time_normalize(value::TaylorPolynomial, reciprocal::TaylorPolynomial, k)
    scale!(value, value, inv(convert(coefficient_type(value), k)))
    return value * reciprocal
end
function time_quotient_coefficient(a, b, c, k, reciprocal)
    value = a[k + 1]
    @inbounds for j in 1:k
        value -= b[j + 1] * c[k - j + 1]
    end
    return time_normalize(value, reciprocal, 1)
end
function time_quotient_coefficient(a::Vector{P}, b::Vector{P}, c::Vector{P}, k, reciprocal) where {P <: TaylorPolynomial}
    value, scratch = copy(a[k + 1]), zero(c[1])
    unit = one(coefficient_type(value))
    @inbounds for j in 1:k
        (iszero(b[j + 1]) || iszero(c[k - j + 1])) && continue
        LinearAlgebra.mul!(scratch, b[j + 1], c[k - j + 1])
        weighted_sum!(value, value, unit, scratch, -unit)
    end
    return time_normalize(value, reciprocal, 1)
end

# A single coefficient of exp(a), log(a), or a^p. The caller checks analyticity
# and supplies 1/a[0] for log/power; exp never reads the reciprocal argument.
function time_unary_coefficient(a, c, k, p, reciprocal, mode::Val{MODE}) where {MODE}
    value = time_weighted_sum(a, c, k, p, mode)
    MODE === :log && (value += k * a[k + 1])
    return MODE === :exp ? value / k : time_normalize(value, reciprocal, k)
end
