# Compose univariate Taylor coefficients with the nonconstant part of a TaylorPolynomial.
# Horner stages need only successively increasing total degrees, and reuse buffers.
function series(a::TaylorPolynomial, coefficients::AbstractVector{T}) where {T}
    ctx = valid(a)
    R = promote_type(coefficient_type(a), T)
    a.len == 1 && return TaylorPolynomial{R}(coefficients[1])
    n = min(ctx.cutoff, length(coefficients) - 1)
    h = TaylorPolynomial{R}(a); h.coeffs[1] = zero(R)
    result, scratch = allocate(ctx, R), allocate(ctx, R)
    result.coeffs[1] = coefficients[n + 1]
    for k in (n - 1):-1:0
        multiply!(scratch, h, result, n - k)
        scratch.coeffs[1] += coefficients[k + 1]
        finish!(scratch, scratch.len)
        result, scratch = scratch, result
    end
    return result
end
function domain_call(f, x)
    return try
        f(x)
    catch err
        err isa DomainError || rethrow()
        throw(TaylorError(sprint(showerror, err)))
    end
end

# Euler's homogeneous derivative E(x^alpha)=|alpha|x^alpha turns
# E(exp(a))=exp(a)E(a), a E(log(a))=E(a), and a E(a^p)=p a^p E(a)
# into triangular coefficient recurrences. One traversal replaces repeated
# polynomial products, with no cutoff applied to intermediate Taylor factors.
function recurrence(a::TaylorPolynomial, c0, p, ::Val{MODE}) where {MODE}
    ctx = valid(a); basis = ctx.basis
    T = promote_type(coefficient_type(a), typeof(c0), typeof(p))
    a.len == 1 && return TaylorPolynomial{T}(c0)
    n = basis.ends[ctx.cutoff + 1]
    out = allocate_undef(ctx, T, n); out.coeffs[1] = c0
    a0 = a.coeffs[1]
    @inbounds for k in 2:n
        degree = basis.degrees[k]
        value = zero(T)
        for t in (basis.splits[k] + 1):basis.splits[k + 1]
            i, j = Int(basis.left[t]), Int(basis.right[t])
            i > a.len && break
            ai = a.coeffs[i]
            iszero(ai) && continue
            weight = MODE === :exp ? basis.degrees[i] : MODE === :log ? -basis.degrees[j] : p * basis.degrees[i] - basis.degrees[j]
            value += weight * ai * out.coeffs[j]
        end
        if MODE === :log
            value += k <= a.len ? degree * a.coeffs[k] : zero(T)
        end
        out.coeffs[k] = MODE === :exp ? value / degree : value / (degree * a0)
    end
    return finish!(out, n)
end

function trig_recurrence(a::TaylorPolynomial, hyperbolic::Bool)
    x = constant_term(a); ctx = a.algebra; basis = ctx.basis
    s0, c0 = hyperbolic ? (sinh(x), cosh(x)) : sincos(x)
    T = typeof(s0)
    a.len == 1 && return TaylorPolynomial{T}(s0), TaylorPolynomial{T}(c0)
    n = basis.ends[ctx.cutoff + 1]
    s, c = allocate_undef(ctx, T, n), allocate_undef(ctx, T, n)
    s.coeffs[1], c.coeffs[1] = s0, c0
    @inbounds for k in 2:n
        sv, cv = zero(T), zero(T)
        for t in (basis.splits[k] + 1):basis.splits[k + 1]
            i, j = Int(basis.left[t]), Int(basis.right[t])
            i > a.len && break
            ai = a.coeffs[i]
            iszero(ai) && continue
            weight = basis.degrees[i] * ai
            sv += weight * c.coeffs[j]
            cv += weight * s.coeffs[j]
        end
        s.coeffs[k] = sv / basis.degrees[k]
        c.coeffs[k] = (hyperbolic ? cv : -cv) / basis.degrees[k]
    end
    return finish!(s, n), finish!(c, n)
end

function Base.exp(a::TaylorPolynomial)
    x = constant_term(a); c0 = exp(x)
    !isempty(a.algebra.basis.products) && return recurrence(a, c0, zero(c0), Val(:exp))
    c = Vector{typeof(c0)}(undef, a.algebra.cutoff + 1); c[1] = c0
    @inbounds for k in 1:(length(c) - 1)
        c[k + 1] = c[k] / k
    end
    return series(a, c)
end
function Base.log(a::TaylorPolynomial)
    x = constant_term(a); c0 = domain_call(log, x)
    iszero(x) && throw(TaylorError("Logarithm requires a nonzero constant"))
    !isempty(a.algebra.basis.products) && return recurrence(a, c0, zero(c0), Val(:log))
    c = Vector{typeof(c0)}(undef, a.algebra.cutoff + 1); c[1] = c0
    factor = one(x)
    @inbounds for k in 1:(length(c) - 1)
        c[k + 1] = factor / k
        factor = -factor
    end
    return series(a / x, c)
end
Base.log2(a::TaylorPolynomial) = log(a) / log(convert(typeof(float(constant_term(a))), 2))
Base.log10(a::TaylorPolynomial) = log(a) / log(convert(typeof(float(constant_term(a))), 10))
Base.log(b::Real, a::TaylorPolynomial) = log(a) / log(b)
Base.log(::Irrational{:ℯ}, a::TaylorPolynomial) = log(a)

function power_series(a::TaylorPolynomial, p::Real, c0)
    x = constant_term(a)
    iszero(x) && throw(TaylorError("Noninteger powers require a nonzero expansion point"))
    !isempty(a.algebra.basis.products) && return recurrence(a, c0, p, Val(:power))
    T = promote_type(typeof(c0), typeof(p))
    c = Vector{T}(undef, a.algebra.cutoff + 1); c[1] = c0
    @inbounds for k in 1:(length(c) - 1)
        c[k + 1] = c[k] * (p - (k - 1)) / k
    end
    return series(a / x, c)
end
function Base.:^(a::TaylorPolynomial, p::AbstractFloat)
    R = promote_type(coefficient_type(a), typeof(p))
    isinteger(p) && typemin(Int32) < p <= typemax(Int32) && return convert(TaylorPolynomial{R}, a)^Int(p)
    iszero(a) && p > 0 && return zero(TaylorPolynomial{R})
    return power_series(a, p, domain_call(x -> x^p, constant_term(a)))
end
Base.:^(a::TaylorPolynomial, p::Rational) = a^convert(float(promote_type(coefficient_type(a), typeof(p))), p)
Base.:^(a::Real, b::TaylorPolynomial) = exp(log(a) * b)
Base.:^(a::TaylorPolynomial, b::TaylorPolynomial) = exp(log(a) * b)
Base.:^(::Irrational{:ℯ}, b::TaylorPolynomial) = exp(b)
@inline coefficient_at(a::TaylorPolynomial, k::Int) = k <= a.len ? a.coeffs[k] : zero(coefficient_type(a))
@inline coefficient_at(a::Real, k::Int) = k == 1 ? a : zero(a)

# b*y=a is triangular by total degree. Division needs one recurrence, avoiding
# both the homogeneous-degree weights of a generic power and a second product.
function quotient(a::Real, b::TaylorPolynomial)
    ctx = a isa TaylorPolynomial ? compatible(a, b) : valid(b)
    b0 = constant_term(b)
    iszero(b0) && throw(TaylorError("Division requires a nonzero denominator constant"))
    if b.len == 1
        result = a / b0
        return result isa TaylorPolynomial ? result : TaylorPolynomial{typeof(result)}(result)
    end
    basis = ctx.basis
    if isempty(basis.products)
        return a * power_series(b, -one(b0), inv(b0))
    end
    T = promote_type(a isa TaylorPolynomial ? coefficient_type(a) : typeof(a), coefficient_type(b))
    T = typeof(one(T) / one(T))
    n = basis.ends[ctx.cutoff + 1]
    result = allocate_undef(ctx, T, n)
    result.coeffs[1] = constant_term(a) / b0
    @inbounds for k in 2:n
        value = convert(T, coefficient_at(a, k))
        for t in (basis.splits[k] + 1):basis.splits[k + 1]
            i, j = Int(basis.left[t]), Int(basis.right[t])
            i > b.len && break
            bi = b.coeffs[i]
            iszero(bi) && continue
            value = muladd(-bi, result.coeffs[j], value)
        end
        result.coeffs[k] = value / b0
    end
    return finish!(result, n)
end
Base.inv(a::TaylorPolynomial) = quotient(one(coefficient_type(a)), a)

# y^2=a: use each unordered pair of nonconstant monomials only once.
function square_root(a::TaylorPolynomial, c0)
    ctx = valid(a); basis = ctx.basis
    a.len == 1 && return TaylorPolynomial(c0)
    n = basis.ends[ctx.cutoff + 1]
    result = allocate_undef(ctx, typeof(c0), n)
    result.coeffs[1] = c0
    denominator = c0 + c0
    @inbounds for k in 2:n
        value = coefficient_at(a, k)
        for t in (basis.splits[k] + 1):basis.splits[k + 1]
            i, j = Int(basis.left[t]), Int(basis.right[t])
            i > j && break
            term = result.coeffs[i] * result.coeffs[j]
            value -= i == j ? term : term + term
        end
        result.coeffs[k] = value / denominator
    end
    return finish!(result, n)
end
function nthroot(a::TaylorPolynomial, p::Integer = 2)
    p != 0 || throw(DomainError(p, "Zeroth root is undefined"))
    p == 1 && return copy(a)
    p > 0 && iszero(a) && return zero(a)
    x = constant_term(a)
    iszero(x) && throw(TaylorError("Root is not analytic at a zero constant"))
    iseven(p) && x < 0 && throw(TaylorError("Even root of a negative constant"))
    exponent = one(float(x)) / p
    c0 = p == 2 ? sqrt(x) : p == 3 ? cbrt(x) : copysign(abs(x)^exponent, x)
    p == 2 && !isempty(a.algebra.basis.products) && return square_root(a, c0)
    return power_series(a, exponent, c0)
end
Base.sqrt(a::TaylorPolynomial) = nthroot(a, 2)
Base.cbrt(a::TaylorPolynomial) = nthroot(a, 3)
function nthroot(x::Real, p::Integer = 2)
    p != 0 || throw(DomainError(p, "Zeroth root is undefined"))
    p == 1 && return x
    exponent = one(float(x)) / p
    return isodd(p) ? copysign(abs(x)^exponent, x) : x^exponent
end

for fn in (:sin, :cos, :sinh, :cosh)
    @eval function Base.$fn(a::TaylorPolynomial)
        if !isempty(valid(a).basis.products)
            pair = trig_recurrence(a, $(fn in (:sinh, :cosh)))
            return pair[$(fn in (:sin, :sinh) ? 1 : 2)]
        end
        x = constant_term(a)
        s, c = $(fn in (:sin, :cos) ? :(sincos(x)) : :((sinh(x), cosh(x))))
        values = $(fn == :sin ? :((s, c, -s, -c)) : fn == :cos ? :((c, -s, -c, s)) : fn == :sinh ? :((s, c, s, c)) : :((c, s, c, s)))
        coeffs = Vector{typeof(s)}(undef, a.algebra.cutoff + 1)
        factor = one(s)
        @inbounds for k in 0:(length(coeffs) - 1)
            coeffs[k + 1] = values[mod(k, 4) + 1] * factor
            factor /= k + 1
        end
        return series(a, coeffs)
    end
end
Base.sincos(a::TaylorPolynomial) = isempty(valid(a).basis.products) ? (sin(a), cos(a)) : trig_recurrence(a, false)
for (fn, sign) in ((:tan, 1), (:tanh, -1))
    @eval function Base.$fn(a::TaylorPolynomial)
        x = constant_term(a); c0 = $fn(x)
        c = zeros(typeof(c0), a.algebra.cutoff + 1); c[1] = c0
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
    iszero(q0) && throw(TaylorError("Singular derivative at the expansion point"))
    c0 = domain_call(x -> x^p, q0)
    c = zeros(typeof(c0), n + 1); c[1] = c0
    @inbounds for k in 0:(n - 1)
        value = (p - k) * q1 * c[k + 1]
        k > 0 && (value += (2p - k + 1) * q2 * c[k])
        c[k + 2] = value / ((k + 1) * q0)
    end
    return c
end
function integral_series(a::TaylorPolynomial, c0, derivative)
    T = promote_type(typeof(c0), eltype(derivative))
    c = Vector{T}(undef, length(derivative) + 1); c[1] = c0
    @inbounds for k in eachindex(derivative)
        c[k + 1] = derivative[k] / k
    end
    return series(a, c)
end
for fn in (:asin, :acos, :atan, :asinh, :acosh, :atanh)
    @eval function Base.$fn(a::TaylorPolynomial)
        x = constant_term(a); u = one(float(x)); c0 = domain_call($fn, x)
        a.len == 1 && return TaylorPolynomial(c0)
        q0, q1, q2, p = $(
            fn in (:asin, :acos) ? :((u - x * x, -2x, -u, -u / 2)) :
                fn == :atan ? :((u + x * x, 2x, u, -u)) : fn == :asinh ? :((u + x * x, 2x, u, -u / 2)) :
                fn == :acosh ? :((x * x - u, 2x, u, -u / 2)) : :((u - x * x, -2x, -u, -u))
        )
        derivative = quadratic_power(q0, q1, q2, p, a.algebra.cutoff - 1)
        $(fn == :acos) && (derivative .*= -one(eltype(derivative)))
        return integral_series(a, c0, derivative)
    end
end
function Base.atan(y::TaylorPolynomial, x::TaylorPolynomial)
    compatible(y, x)
    x0, y0 = constant_term(x), constant_term(y)
    iszero(x0) && iszero(y0) && throw(TaylorError("atan is not analytic at the origin"))
    p = abs(x0) >= abs(y0) ? atan(y / x) : -atan(x / y)
    return p + (atan(y0, x0) - constant_term(p))
end
Base.atan(a::TaylorPolynomial, b::Real) = atan(promote(a, b)...)
Base.atan(a::Real, b::TaylorPolynomial) = atan(promote(a, b)...)
function Base.hypot(a::TaylorPolynomial, b::TaylorPolynomial)
    compatible(a, b)
    scale = max(abs(constant_term(a)), abs(constant_term(b)))
    iszero(scale) && return sqrt(a * a + b * b)
    # Normalize before squaring to avoid overflow/underflow at finite centers.
    return scale * sqrt((a / scale)^2 + (b / scale)^2)
end
Base.hypot(a::TaylorPolynomial, b::Real) = hypot(promote(a, b)...)
Base.hypot(a::Real, b::TaylorPolynomial) = hypot(promote(a, b)...)
for fn in (:round, :trunc)
    @eval function Base.$fn(a::TaylorPolynomial)
        out = copy(a)
        out.coeffs[1] = $fn(constant_term(a))
        return finish!(out, out.len)
    end
end
function Base.mod(a::TaylorPolynomial, p::Real)
    out = TaylorPolynomial{promote_type(coefficient_type(a), typeof(p))}(a)
    out.coeffs[1] = mod(constant_term(a), p)
    return finish!(out, out.len)
end

for fn in (:erf, :erfc)
    @eval const $fn = SpecialFunctions.$fn
    @eval function SpecialFunctions.$fn(a::TaylorPolynomial)
        x = constant_term(a); c0 = SpecialFunctions.$fn(x)
        n = a.algebra.cutoff
        derivative = zeros(typeof(c0), n)
        derivative[1] = $(fn == :erf ? 2 : -2) * exp(-x * x) / sqrt(convert(typeof(c0), π))
        @inbounds for k in 1:(n - 1)
            derivative[k + 1] = (-2x * derivative[k] - (k > 1 ? 2derivative[k - 1] : zero(c0))) / k
        end
        return integral_series(a, c0, derivative)
    end
end
const loggamma = SpecialFunctions.loggamma
const gamma = SpecialFunctions.gamma
function derivative_series(f, a::TaylorPolynomial, c0)
    a.len == 1 && return TaylorPolynomial(c0)
    x = constant_term(a)
    c = Vector{typeof(c0)}(undef, a.algebra.cutoff + 1); c[1] = c0
    factor = one(c0)
    for k in 1:(length(c) - 1)
        factor /= k
        c[k + 1] = f(k, x) * factor
    end
    return series(a, c)
end
gamma_series(a::TaylorPolynomial, c0) = derivative_series((k, x) -> scalar_psi(k - 1, x), a, c0)
SpecialFunctions.loggamma(a::TaylorPolynomial) = gamma_series(a, domain_call(SpecialFunctions.loggamma, constant_term(a)))
function SpecialFunctions.gamma(a::TaylorPolynomial)
    value, sign = domain_call(SpecialFunctions.logabsgamma, constant_term(a))
    return sign * exp(gamma_series(a, value))
end
function polygamma_series(a::TaylorPolynomial, n::Integer)
    n >= 0 || throw(ArgumentError("Invalid polygamma order"))
    x = constant_term(a); c0 = scalar_psi(n, x)
    return derivative_series((k, x) -> scalar_psi(n + k, x), a, c0)
end
SpecialFunctions.polygamma(n::Integer, a::TaylorPolynomial) = polygamma_series(a, n)
SpecialFunctions.digamma(a::TaylorPolynomial) = polygamma_series(a, 0)

# Bessel's second-order differential equation supplies all higher derivatives.
# With y=exp(s*x)*z, scaled I/K remain scaled throughout (no overflow-prone unscaling).
function bessel_series(fn, n::Integer, a::TaylorPolynomial, q::Int, s)
    typemin(Int32) < n <= typemax(Int32) || throw(ArgumentError("Bessel order out of range"))
    x = constant_term(a); c0 = scalar_bessel(fn, n, x)
    a.len == 1 && return TaylorPolynomial(c0)
    N = a.algebra.cutoff
    if iszero(x)
        fn in (SpecialFunctions.besselj, SpecialFunctions.besseli) || throw(TaylorError("Bessel expansion is singular at zero"))
        c = zeros(typeof(c0), N + 1)
        order = abs(n)
        if order <= N
            sign = n < 0 && fn === SpecialFunctions.besselj && isodd(order) ? -1 : 1
            value = convert(typeof(c0), sign) / convert(typeof(c0), big(2)^order * factorial(big(order)))
            for k in order:2:N
                c[k + 1] = value
                m = (k - order) ÷ 2 + 1
                value *= -q / (convert(typeof(c0), 4) * m * (m + order))
            end
        end
        return series(a, c)
    end
    c = zeros(typeof(c0), N + 1); c[1] = c0
    lower, upper = scalar_bessel(fn, n - 1, x), scalar_bessel(fn, n + 1, x)
    neighbor = fn in (SpecialFunctions.besselj, SpecialFunctions.bessely) ? (lower - upper) / 2 :
        fn in (SpecialFunctions.besseli, SpecialFunctions.besselix) ? (lower + upper) / 2 : -(lower + upper) / 2
    c[2] = neighbor - s * c0
    r = s * s + q
    @inbounds for k in 0:(N - 2)
        value = -(2x * k + x + 2s * x * x) * (k + 1) * c[k + 2] - (k * k + 4s * x * k + r * x * x + s * x - n * n) * c[k + 1]
        k >= 1 && (value -= (2s * (k - 1) + 2r * x + s) * c[k])
        k >= 2 && (value -= r * c[k - 1])
        c[k + 3] = value / (x * x * (k + 2) * (k + 1))
    end
    return series(a, c)
end
for (fn, q) in ((:besselj, 1), (:bessely, 1), (:besseli, -1), (:besselk, -1), (:besselix, -1), (:besselkx, -1))
    @eval const $fn = SpecialFunctions.$fn
    @eval SpecialFunctions.$fn(n::Integer, a::TaylorPolynomial) = bessel_series(
        SpecialFunctions.$fn, n, a, $q,
        $(fn == :besselix ? :(sign(constant_term(a))) : fn == :besselkx ? :(-one(constant_term(a))) : :(zero(constant_term(a))))
    )
end
