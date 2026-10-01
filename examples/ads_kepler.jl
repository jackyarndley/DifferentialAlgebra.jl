# # Automatic domain splitting for a Kepler map
#
# Propagate a box of elliptic orbits from periapsis for three nominal periods.
# The uncertain parameters are semimajor axis `a` and eccentricity `e`, with
# gravitational parameter μ = 1. A single Taylor expansion becomes inaccurate
# across this box because mean motion depends on `a`. Automatic domain splitting
# (ADS) replaces it with smaller, locally expanded maps.
#
# ADS originates in the uncertainty-propagation work of [Wittig2015](@citet),
# DOI: 10.1007/s10569-015-9618-3. Here the complete Kepler map is recomputed on
# each child domain using the analytic flow.

using DifferentialAlgebra
using Random

# Solve E - e sin(E) = M. Newton iteration first finds the scalar root, then
# lifts it to the active polynomial order. The same function handles numeric
# validation points and Taylor inputs.
function eccentric_anomaly(mean, eccentricity)
    m, e = constant_term(mean), constant_term(eccentricity)
    anomaly = m
    converged = false
    for _ in 1:40
        s, c = sincos(anomaly)
        step = (anomaly - e * s - m) / (1 - e * c)
        anomaly -= step
        if abs(step) <= 8eps(typeof(m)) * max(one(m), abs(anomaly))
            converged = true
            break
        end
    end
    converged || error("Kepler iteration did not converge")
    mean isa TaylorPolynomial || return anomaly
    expansion = zero(mean) + anomaly
    for _ in 1:ceil(Int, log2(truncation_order() + 1))
        s, c = sincos(expansion)
        expansion -= (expansion - eccentricity * s - mean) / (1 - eccentricity * c)
    end
    return expansion
end

# The planar state is [x, y, vx, vy]. Length and time use consistent units.
function kepler_map(elements; time = 6π)
    a, e = elements
    mean = time / (a * sqrt(a))
    anomaly = eccentric_anomaly(mean, e)
    s, c = sincos(anomaly)
    beta = sqrt(1 - e^2)
    radius = a * (1 - e * c)
    return [a * (c - e), a * beta * s, -sqrt(a) * s / radius, sqrt(a) * beta * c / radius]
end

# `adaptive_map` takes physical bounds. Every patch has its own center and
# scale; calling the result selects the appropriate patch automatically.
lower, upper = [0.95, 0.2], [1.05, 0.5]
order, tolerance = 5, 1.0e-6
single = adaptive_map(kepler_map, lower, upper; order, atol = tolerance, max_depth = 0, strict = false)
split = adaptive_map(kepler_map, lower, upper; order, atol = tolerance, names = (:a, :e))
@assert split.converged && length(split.patches) > 1

# Check an independent seeded sample and a regular grid, including the box
# boundary. ADS estimates are heuristic, so validation is separate from its
# acceptance rule. The energy check also verifies the Kepler state formulas.
rng = Xoshiro(2026)
samples = [lower + (upper - lower) .* rand(rng, 2) for _ in 1:1000]
append!(samples, [[a, e] for a in range(lower[1], upper[1]; length = 31), e in range(lower[2], upper[2]; length = 31)])
reference = kepler_map.(samples)
single_error = maximum(maximum(abs.(single(x) - y)) for (x, y) in zip(samples, reference))
split_error = maximum(maximum(abs.(split(x) - y)) for (x, y) in zip(samples, reference))
@assert split_error < tolerance
@assert split_error < single_error / 100
for (elements, state) in zip(samples, reference)
    radius = hypot(state[1], state[2])
    energy = (state[3]^2 + state[4]^2) / 2 - 1 / radius
    @assert abs(energy + 1 / (2elements[1])) < 1.0e-13
end
println((patches = length(split.patches), single_error = single_error, split_error = split_error))

# ## Inspect the adaptive map
#
# Plot the partition in parameter space and its image in the orbital plane.
# The colored curves are images of patch boundaries, not trajectory segments.
# The lower panels compare a single expansion with the ADS map on the same grid
# and color scale. Errors include all four state components, not just position.
using CairoMakie

semimajor_axes = range(lower[1], upper[1]; length = 81)
eccentricities = range(lower[2], upper[2]; length = 61)
map_error(map, elements) = maximum(abs, map(elements) - kepler_map(elements))
single_errors = [map_error(single, [a, e]) for a in semimajor_axes, e in eccentricities]
split_errors = [map_error(split, [a, e]) for a in semimajor_axes, e in eccentricities]
@assert maximum(split_errors) < tolerance

fig = Figure(size = (1150, 850), fontsize = 15)
domain_axis = Axis(fig[1, 1]; xlabel = "Semimajor axis a", ylabel = "Eccentricity e", title = "ADS: $(length(split.patches)) patches")
image_axis = Axis(fig[1, 2]; xlabel = "x", ylabel = "y", title = "Propagated uncertainty", aspect = DataAspect())
colors = Makie.wong_colors()
for (i, patch) in enumerate(split.patches)
    a0, e0 = patch.lower
    a1, e1 = patch.upper
    color = colors[mod1(i, length(colors))]
    lines!(domain_axis, [a0, a1, a1, a0, a0], [e0, e0, e1, e1, e0]; color, linewidth = 1)
    boundary = vcat(
        [[a, e0] for a in range(a0, a1; length = 20)],
        [[a1, e] for e in range(e0, e1; length = 20)],
        [[a, e1] for a in range(a1, a0; length = 20)],
        [[a0, e] for e in range(e1, e0; length = 20)],
    )
    states = split.(boundary)
    lines!(image_axis, first.(states), getindex.(states, 2); color, linewidth = 0.8)
end
single_axis = Axis(fig[2, 1]; xlabel = "Semimajor axis a", ylabel = "Eccentricity e", title = "Single map: max error $(round(maximum(single_errors); sigdigits = 2))")
split_axis = Axis(fig[2, 2]; xlabel = "Semimajor axis a", ylabel = "Eccentricity e", title = "ADS map: max error $(round(maximum(split_errors); sigdigits = 2))")
color_limits = (-12.0, max(-6.0, ceil(log10(maximum(single_errors)))))
heatmap!(single_axis, semimajor_axes, eccentricities, log10.(max.(single_errors, 1.0e-12)); colormap = :viridis, colorrange = color_limits)
heat = heatmap!(split_axis, semimajor_axes, eccentricities, log10.(max.(split_errors, 1.0e-12)); colormap = :viridis, colorrange = color_limits)
Colorbar(fig[2, 3], heat; label = "log₁₀ maximum component error")

fig
#
# Errors below 10⁻¹² share the color scale's lower limit. The estimates that drive
# splitting are heuristic; the plotted errors come from independent pointwise
# evaluations. Tightening `atol` generally produces more patches.
#
