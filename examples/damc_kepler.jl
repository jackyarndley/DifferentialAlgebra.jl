# # Differential algebra Monte Carlo for Kepler motion
#
# Propagate the same uncertain orbit as DACE.jl's DAMC example: 10,000 samples,
# 30 nominal revolutions, and polynomial orders 2, 4 and 8. The method follows
# Valli et al. (2013), DOI: 10.2514/1.58068. The initial covariance is diagonal.
#
# A Taylor map is propagated once at each order, then evaluated for every sample.
# Ordinary vectors hold the Cartesian states.

using DifferentialAlgebra
using LinearAlgebra, Random, Statistics
using CairoMakie

# Solve the elliptic or hyperbolic Kepler equation for the anomaly increment.
# First converge the scalar center, then lift that root in the polynomial algebra.
# Each Newton step doubles the number of correct Taylor orders.
function kepler_residual(anomaly, mean, sigma, eccentric_term, hyperbolic)
    if hyperbolic
        s, c = sinh(anomaly), cosh(anomaly)
        return -mean - anomaly + sigma * (c - 1) + eccentric_term * s,
            -1 + sigma * s + eccentric_term * c
    end
    s, c = sincos(anomaly)
    return -mean + anomaly + sigma * (1 - c) - eccentric_term * s,
        1 + sigma * s - eccentric_term * c
end

function kepler_increment(mean, sigma, eccentric_term, hyperbolic)
    m, s, e = constant_term.((mean, sigma, eccentric_term))
    anomaly = hyperbolic ? zero(m) : m
    converged = false
    for _ in 1:40
        residual, derivative = kepler_residual(anomaly, m, s, e, hyperbolic)
        correction = residual / derivative
        anomaly -= correction
        if abs(correction) <= 8eps(typeof(m)) * max(one(m), abs(anomaly))
            converged = true
            break
        end
    end
    converged || error("Kepler iteration did not converge")
    mean isa TaylorPolynomial || return anomaly
    expansion = zero(mean) + anomaly
    for _ in 1:ceil(Int, log2(truncation_order() + 1))
        residual, derivative = kepler_residual(expansion, mean, sigma, eccentric_term, hyperbolic)
        expansion -= residual / derivative
    end
    return expansion
end

# Lagrange f and g coefficients propagate Cartesian position and velocity.
# The formulas cover both elliptic and hyperbolic orbits, away from the
# parabolic limit. Time and gravitational parameter use consistent units.
function lagrange_propagator(state, Δt, μ)
    r0 = state[1:3]
    v0 = state[4:6]
    radius0 = norm(r0)
    a = μ / (2μ / radius0 - dot(v0, v0))
    sigma0 = dot(r0, v0) / sqrt(μ)
    hyperbolic = constant_term(a) < 0
    scale = hyperbolic ? sqrt(-a) : sqrt(a)
    mean = Δt * sqrt(μ) / (scale^3)
    anomaly = kepler_increment(mean, sigma0 / scale, 1 - radius0 / a, hyperbolic)
    s, c = hyperbolic ? (sinh(anomaly), cosh(anomaly)) : sincos(anomaly)
    f = 1 - a / radius0 * (1 - c)
    g = a * sigma0 / sqrt(μ) * (1 - c) + radius0 * scale / sqrt(μ) * s
    position = f * r0 + g * v0
    radius = norm(position)
    fdot = -sqrt(μ) * scale / (radius * radius0) * s
    gdot = 1 - a / radius * (1 - c)
    velocity = fdot * r0 + gdot * v0
    return vcat(position, velocity)
end

# Initial state and covariance from the original example.
x0 = [-0.68787, -0.39713, 0.28448, -0.51331, 0.98266, 0.37611]
variance = [1.0e-7, 1.0e-7, 1.0e-7, 1.0e-9, 1.0e-9, 1.0e-9]
μ, Δt = 1.0, 30 * 2π
scale, nsamples = 3.0, 10_000
rng = Xoshiro(2026)
samples = randn(rng, 6, nsamples)
initial_samples = x0 .+ sqrt.(variance) .* samples
nominal = lagrange_propagator(x0, Δt, μ)
energy(state) = dot(state[4:6], state[4:6]) / 2 - μ / norm(state[1:3])
@assert abs(energy(nominal) - energy(x0)) < 1.0e-12
monte_carlo = stack(lagrange_propagator(x, Δt, μ) for x in eachcol(initial_samples));

# The independent coordinates are scaled by three standard deviations.
perturbations = variables(6; order = 8)
initial_map = x0 + scale * sqrt.(variance) .* perturbations
approximations = Dict{Int, Matrix{Float64}}()
errors = Dict{Int, Float64}()
for order in (2, 4, 8)
    propagated = with_order(order) do
        lagrange_propagator(initial_map, Δt, μ)
    end
    @assert maximum(abs, constant_term(propagated) - nominal) < 1.0e-10
    map = CompiledMap(propagated)
    values = Matrix{Float64}(undef, 6, nsamples)
    work = Vector{Float64}(undef, degree(map) + 1)
    point = Vector{Float64}(undef, 6)
    for (j, sample) in enumerate(eachcol(samples))
        point .= sample ./ scale
        evaluate!(@view(values[:, j]), map, point, work)
    end
    approximations[order] = values
    errors[order] = sqrt(mean(abs2, values - monte_carlo))
    println("Order ", order, " RMS error: ", errors[order])
end
@assert errors[8] < errors[4] < errors[2]
@assert truncation_order() == 8

# Check the hyperbolic branch by reversing a short propagation.
escape_state = [1.0, 0.0, 0.0, 0.0, 2.0, 0.0]
escaped = lagrange_propagator(escape_state, 0.25, μ)
@assert norm(lagrange_propagator(escaped, -0.25, μ) - escape_state) < 1.0e-12

# Compare the final position distributions and save the same PDF as DACE.jl.
figure = Figure()
axis = Axis(figure[1, 1]; xlabel = "x [-]", ylabel = "y [-]")
scatter!(axis, monte_carlo[1, :], monte_carlo[2, :]; marker = :diamond, color = :black, label = "MC")
for (order, marker, color) in ((2, :dtriangle, "#666666"), (4, :circle, "#999999"), (8, :utriangle, "#CCCCCC"))
    values = approximations[order]
    scatter!(axis, values[1, :], values[2, :]; marker, color, label = "DAMC-$order")
end
axislegend(axis; position = :rt)
save("damc_kepler.pdf", figure)
save("damc_kepler.png", figure)

# ![Monte Carlo and polynomial propagation of the final position distribution.](damc_kepler.png)
#
# [Download the figure as a PDF.](damc_kepler.pdf)
