# Polynomial SpecialFunctions methods share the elementary composition kernels.
# Scalar seeds, including precision-preserving fallbacks, live in special_scalars.jl.
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
    primal = coefficient_primal(x)
    primal <= 0 && isinteger(primal) &&
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
    return scaled(
        exponential(gamma_series(a, value; filter = false); filter = false),
        isnan(value) ? one(sign) : sign
    )
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
    limits = (-peak * width, zero(T), T(Inf))
    coefficients = bessel_moments(integrand!, T, count, limits)
    result = series(normalized_argument(a, x), coefficients; filter = false)
    return scaled(result, c0)
end

function bessel_moments(integrand!, ::Type{T}, count, limits; norm = nothing) where {T}
    quadrature_order = max(15, cld(precision(T), 8))
    if norm === nothing
        approximate, _ = QuadGK.quadgk!(
            integrand!, zeros(T, count), limits...;
            order = quadrature_order, rtol = sqrt(eps(T))
        )
        # Positive K moments can span very different scales.
        scales = map(c -> iszero(c) ? one(T) : abs(c), approximate)
        norm = v -> maximum(i -> abs(v[i]) / scales[i], eachindex(v))
    end
    coefficients, error = QuadGK.quadgk!(
        integrand!, zeros(T, count), limits...;
        order = quadrature_order, rtol = 8eps(T), norm
    )
    error <= 64eps(T) * norm(coefficients) ||
        throw(TaylorError("Scaled Bessel derivative quadrature did not converge"))
    coefficients ./= first(coefficients)
    return coefficients
end

# DLMF 10.32.3: differentiating e^-x I_n(x) inserts (cos(theta)-1)^k.
# u=theta*sqrt(x) resolves the peak, including at huge centers. A breakpoint
# beyond the highest moment's peak keeps a long, negligible tail from hiding it.
# Normalizing the moments cancels the common Jacobian factor 1/sqrt(x).
function scaled_besseli_series(n, a, c0::T) where {T <: AbstractFloat}
    x = coefficient_convert(T, constant_term(a))
    root = sqrt(abs(x))
    count = working_order(a.algebra) + 1
    function integrand!(values, u)
        theta = u / root
        decay = 2(root * sin(theta / 2))^2
        if !isfinite(decay)
            fill!(values, zero(T))
            return nothing
        end
        values[1] = exp(-decay) * cos(T(n) * theta)
        for k in 1:(count - 1)
            values[k + 1] = -decay * values[k] / k
        end
        return nothing
    end
    upper = T(π) * root
    boundary = 2sqrt(T(count + precision(T)))
    limits = boundary < upper ? (zero(T), boundary, upper) : (zero(T), upper)
    # I moments can cancel through zero. An absolute norm avoids forcing
    # relative accuracy on those components and accumulating roundoff by refinement.
    coefficients = bessel_moments(integrand!, T, count, limits; norm = v -> maximum(abs, v))
    result = series(normalized_argument(a, x), coefficients; filter = false)
    return scaled(result, c0)
end

# Bessel's second-order differential equation supplies the singular Y/K solutions.
# With y=exp(s*x)*z, scaled I/K remain scaled throughout (no overflow-prone unscaling).
function bessel_series(fn::F, n::Integer, a::TaylorPolynomial, q::Int, s) where {F}
    typemin(Int32) < n <= typemax(Int32) || throw(ArgumentError("Bessel order out of range"))
    x = constant_term(a)
    !invertible_constant(x) &&
        fn in (
        SpecialFunctions.bessely, SpecialFunctions.besselk, SpecialFunctions.besselkx,
    ) &&
        throw(TaylorError("Bessel expansion is singular at zero"))
    if fn === SpecialFunctions.besselix && !invertible_constant(x) &&
            (a.len > 1 || !coefficient_isequal(x, coefficient_primal(x)))
        throw(TaylorError("Scaled Bessel I is not analytic at zero"))
    end
    c0 = domain_call(x -> scalar_bessel(fn, n, x), x)
    a.len == 1 && return TaylorPolynomial{typeof(float(c0))}(valid(a), c0)
    fn === SpecialFunctions.besselix && x isa AbstractFloat && isfinite(x) &&
        abs(x) >= 2 && abs(n) <= max(2, sqrt(abs(x))) &&
        return scaled_besseli_series(n, a, c0)
    if fn in (SpecialFunctions.besselj, SpecialFunctions.besseli, SpecialFunctions.besselix)
        return regular_bessel_series(fn, n, a, c0, s)
    end
    fn === SpecialFunctions.besselkx && x isa AbstractFloat && isfinite(x) && x >= 2 &&
        return large_besselk_series(n, a, c0)
    cutoff = working_order(a.algebra)
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
