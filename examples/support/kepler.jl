# Scalar elliptic Kepler solve shared by the implicit-map, continuation and ADS
# tutorials. Bracketing keeps Newton steps on the unique real elliptic branch.
function kepler_root(mean, eccentricity)
    m, e = promote(float(mean), float(eccentricity))
    isfinite(m) && isfinite(e) && 0 <= e < 1 ||
        throw(DomainError((mean, eccentricity), "Use a finite mean anomaly and 0 ≤ e < 1"))
    lower, upper = m - e, m + e
    anomaly = m
    for _ in 1:80
        sine, cosine = sincos(anomaly)
        residual = anomaly - e * sine - m
        iszero(residual) && return anomaly
        if residual < 0
            lower = anomaly
        else
            upper = anomaly
        end
        correction = residual / (1 - e * cosine)
        candidate = anomaly - correction
        if abs(correction) <= 4eps(one(m)) * max(one(m), abs(anomaly))
            return candidate
        end
        anomaly = lower < candidate < upper ? candidate : lower + (upper - lower) / 2
    end
    error("Kepler root did not converge")
end

# Central gravity in two or three dimensions; state ordering is [position; velocity].
function two_body!(du, u, μ, t)
    dimensions = length(u) ÷ 2
    radius2 = sum(x -> x^2, @view(u[1:dimensions]))
    factor = -μ / (radius2 * sqrt(radius2))
    for i in 1:dimensions
        du[i] = u[i + dimensions]
        du[i + dimensions] = factor * u[i]
    end
    return nothing
end

function two_body(u, μ, t)
    dimensions = length(u) ÷ 2
    position = @view(u[1:dimensions])
    radius2 = sum(x -> x^2, position)
    factor = -μ / (radius2 * sqrt(radius2))
    return vcat(u[(dimensions + 1):end], factor .* position)
end

# Cartesian state from an elliptic orbit's eccentric anomaly. The rotation is
# Rz(Ω) Rx(i) Rz(ω); all arguments may carry Taylor or interval uncertainty.
function element_state(a, e, i, Ω, ω, E, μ)
    sine, cosine = sincos(E)
    beta = sqrt(one(e) - e^2)
    speed = sqrt(μ / a) / (one(e) - e * cosine)
    sw, cw = sincos(ω)
    si, ci = sincos(i)
    sn, cn = sincos(Ω)
    function rotate(x, y)
        along = cw * x - sw * y
        across = sw * x + cw * y
        return [cn * along - sn * ci * across, sn * along + cn * ci * across, si * across]
    end
    return vcat(
        rotate(a * (cosine - e), a * beta * sine),
        rotate(-speed * sine, speed * beta * cosine)
    )
end
