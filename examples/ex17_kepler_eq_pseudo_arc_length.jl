# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex17_kepler_eq_pseudo_arc_length.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # Pseudo-arclength continuation
#
# Parameterize the Kepler solution curve with a local coordinate s.
# At the nominal root, its tangent satisfies fE tE + fe te = 0.
# Normalize that tangent, then impose tE δE + te δe = s.
# This is pseudo-arclength: projection on the fixed nominal unit tangent,
# not exact arc length along the curved solution branch.
using DifferentialAlgebra
using LinearAlgebra
using CairoMakie

function nominal_root(mean, e)
    E = mean
    for _ in 1:30
        step = (E - e * sin(E) - mean) / (1 - e * cos(E))
        E -= step
        abs(step) < 2.0e-15 && return E
    end
    error("Nominal root did not converge")
end

δE, δe, s = variables((:δE, :δe, :s); order = 10)
mean, e0 = π / 2, 0.3
E0 = nominal_root(mean, e0)
slope = sin(E0) / (1 - e0 * cos(E0))
tangent = [slope, 1.0] / hypot(slope, 1.0)
E, e = E0 + δE, e0 + δe
constraints = [E - e * sin(E) - mean, tangent[1] * δE + tangent[2] * δe - s, s]
inverse = invert(constraints)
increments = evaluate(inverse, [zero(s), zero(s), s])
Ecurve, ecurve = E0 + increments[1], e0 + increments[2]
residual = coefficient_norm(Ecurve - ecurve * sin(Ecurve) - mean)
@assert residual < 1.0e-12
@assert norm(tangent) ≈ 1
@assert [coefficient(Ecurve, [0, 0, 1]), coefficient(ecurve, [0, 0, 1])] ≈ tangent
println("E(s) = ", Ecurve)
println("e(s) = ", ecurve)
println((unit_tangent = tangent, coefficient_residual = residual))

# Coefficient-based radii are estimates. Back off until sampled residuals
# meet the requested tolerance, including both ends of the continuation interval.
function continuation_interval(Ecurve, ecurve, mean; tolerance = 1.0e-10)
    radius = min(
        0.3, DifferentialAlgebra.convergence_radius(Ecurve, tolerance),
        DifferentialAlgebra.convergence_radius(ecurve, tolerance)
    )
    for _ in 1:30
        points = range(-radius, radius; length = 101)
        anomalies = [Ecurve([0.0, 0.0, t]) for t in points]
        eccentricities = [ecurve([0.0, 0.0, t]) for t in points]
        errors = abs.(anomalies - eccentricities .* sin.(anomalies) .- mean)
        maximum(errors) <= tolerance && return (; points, anomalies, eccentricities, errors, radius)
        radius /= 2
    end
    error("Continuation interval did not meet the residual tolerance")
end
curve = continuation_interval(Ecurve, ecurve, mean)
println((half_width = curve.radius, maximum_residual = maximum(curve.errors)))

#-
fig = Figure(size = (1000, 390), fontsize = 15)
ax = Axis(fig[1, 1]; xlabel = "Eccentricity", ylabel = "Eccentric anomaly (rad)", title = "Local continuation curve")
lines!(ax, curve.eccentricities, curve.anomalies; linewidth = 3)
scatter!(ax, [e0], [E0]; color = :black, label = "Expansion center")
axislegend(ax; position = :lt)
err = Axis(fig[1, 2]; xlabel = "Pseudo-arclength s", ylabel = "Kepler residual", yscale = log10)
lines!(err, curve.points, max.(curve.errors, eps(Float64)); linewidth = 2)
fig
