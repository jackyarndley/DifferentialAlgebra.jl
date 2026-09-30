# Scalar seeds missing from SpecialFunctions' BigFloat API. These private
# methods preserve the caller's precision and never change global MPFR settings.
scalar_psi(n::Integer, x::Real) = n == 0 ? SpecialFunctions.digamma(x) : SpecialFunctions.polygamma(n, x)

function scalar_psi(n::Integer, x::BigFloat)
    n == 0 && return SpecialFunctions.digamma(x)
    n > 0 || throw(ArgumentError("Polygamma order must be nonnegative"))
    isfinite(x) || return x == Inf ? zero(x) : oftype(x, NaN)
    x <= 0 && isinteger(x) && throw(DomainError(x, "Polygamma has a pole at nonpositive integers"))
    if x < 0
        # Reflection, DLMF 5.15.6. c[k+1] is the Taylor coefficient of
        # cot(pi*(x+t)); c' = -pi*(1+c^2), avoiding numerical differentiation.
        piT = BigFloat(π)
        c = zeros(BigFloat, n + 1)
        c[1] = cospi(x) / sinpi(x)
        for k in 0:(n - 1)
            value = k == 0 ? one(x) : zero(x)
            for j in 0:k
                value += c[j + 1] * c[k - j + 1]
            end
            c[k + 2] = -piT * value / (k + 1)
        end
        return (-1)^n * scalar_psi(n, 1 - x) - piT * c[end] * factorial(big(n))
    end
    # psi^(n)(x) = (-1)^(n+1) n! zeta(n+1,x), DLMF 25.11.12.
    # Euler-Maclaurin after shifting to a well-conditioned positive argument.
    s = n + 1
    z = x
    value = zero(x)
    target = precision(BigFloat) / 2 + s
    while z < target
        value += z^(-s)
        z += 1
    end
    value += z^(1 - s) / (s - 1) + z^(-s) / 2
    pi2 = (2BigFloat(π))^2
    factor = BigFloat(s) * z^(-s - 1) / pi2
    previous = BigFloat(Inf)
    for k in 1:precision(BigFloat)
        term = 2SpecialFunctions.zeta(BigFloat(2k)) * factor
        abs(term) < previous || throw(TaylorError("Polygamma asymptotic sum did not converge"))
        value += isodd(k) ? term : -term
        abs(term) <= eps(BigFloat) * abs(value) && return (isodd(n) ? value : -value) * factorial(big(n))
        previous = abs(term)
        factor *= BigFloat(s + 2k - 1) * (s + 2k) / (z * z * pi2)
    end
    throw(TaylorError("Polygamma sum did not converge"))
end

scalar_bessel(fn, n::Integer, x::Real) = fn(n, x)
# AMOS promotes Float32 arguments. Round scalar seeds back to their requested
# coefficient precision rather than promoting the entire polynomial to Float64.
scalar_bessel(fn, n::Integer, x::Float32) = Float32(fn(n, x))
function scalar_bessel(fn, n::Integer, x::BigFloat)
    fn in (SpecialFunctions.besselj, SpecialFunctions.bessely) && return fn(n, x)
    n = abs(n)
    scaled = fn in (SpecialFunctions.besselix, SpecialFunctions.besselkx)
    if fn in (SpecialFunctions.besseli, SpecialFunctions.besselix)
        return big_besseli(n, x, scaled)
    end
    return big_besselk(n, x, scaled)
end

function big_besseli(n::Integer, x::BigFloat, scaled::Bool)
    ax = abs(x)
    iszero(ax) && return n == 0 ? one(x) : zero(x)
    isnan(ax) && return ax
    isinf(ax) && return scaled ? zero(x) : (isodd(n) && x < 0 ? -ax : ax)
    sign = x < 0 && isodd(n) ? -1 : 1
    if ax > max(precision(BigFloat), 2big(n)^2)
        # Scaled large-argument expansion, DLMF 10.40.1 and 10.17.1.
        value = term = one(x)
        for k in 1:precision(BigFloat)
            term *= (big(2k - 1)^2 - 4big(n)^2) / (8ax * k)
            value += term
            abs(term) <= eps(BigFloat) * abs(value) && break
            k == precision(BigFloat) && throw(TaylorError("Bessel I expansion did not converge"))
        end
        value *= sign / sqrt(2BigFloat(π) * ax)
        return scaled ? value : exp(ax) * value
    end
    # Positive-term convergent series, DLMF 10.25.2 (no cancellation).
    term = (ax / 2)^n / factorial(big(n))
    value = term
    k = 0
    while true
        k += 1
        term *= (ax / 2)^2 / (BigFloat(k) * (n + k))
        value += term
        abs(term) <= eps(BigFloat) * abs(value) && break
    end
    return sign * (scaled ? exp(-ax) * value : value)
end

function big_besselk(n::Integer, x::BigFloat, scaled::Bool)
    x < 0 && throw(DomainError(x, "Real Bessel K requires a positive argument"))
    iszero(x) && return BigFloat(Inf)
    !isfinite(x) && return isnan(x) ? x : zero(x)
    # Positive integral, DLMF 10.32.9. Center at the peak and rescale its width
    # so quadrature also resolves large x/order. This avoids the cancellation
    # of the integer-order power series and the exp(x)*K(x) overflow problem.
    peak = asinh(BigFloat(n) / x)
    width = sqrt(max(hypot(x, BigFloat(n)), one(x)))
    logcosh(z) = abs(z) + log1p(exp(-2abs(z))) - log(BigFloat(2))
    function exponent(t)
        decay = 2x * sinh(t / 2)^2
        return isinf(decay) ? -decay : -decay + logcosh(n * t)
    end
    shift = exponent(peak)
    integrand(u) = exp(exponent(peak + u / width) - shift) / width
    # Increase the quadrature order with precision: order seven needs excessive
    # subdivision when the requested error is around 1e-75.
    value, error = QuadGK.quadgk(
        integrand, -peak * width, zero(x), BigFloat(Inf);
        rtol = 64eps(BigFloat), order = max(15, cld(precision(BigFloat), 8))
    )
    error <= 128eps(BigFloat) * abs(value) || throw(TaylorError("Bessel K quadrature did not converge"))
    return exp(scaled ? shift : shift - x) * value
end
