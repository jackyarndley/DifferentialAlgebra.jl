# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # Newton lifting and fixed-point iteration
#
# Both algorithms solve Kepler's implicit equation in a truncated algebra.
# Newton iteration doubles the number of correct orders near a scalar root;
# fixed-point iteration E ← M + e sin(E) converges linearly here.
# Monitor every polynomial coefficient, not only the nominal residual.
using DifferentialAlgebra
using CairoMakie

set_theme!(palette = (color = Makie.to_colormap(:tab10),))

function kepler_iterations(mean, e; newton, tolerance = 1.0e-12)
    E = copy(mean)
    errors = Float64[]
    for _ in 1:100
        residual = E - e * sin(E) - mean
        push!(errors, coefficient_norm(residual))
        errors[end] <= tolerance && return E, errors
        E = newton ? E - residual / (1 - e * cos(E)) : mean + e * sin(E)
    end
    error("Polynomial iteration did not converge")
end

δa, δe = variables((:δa, :δe); order = 5)
a, e = 1 + δa, 0.3 + δe
period = 2π * a * sqrt(a)
mean_motion = 1 / (a * sqrt(a))

# A quarter of each orbit's own period gives M=π/2 for every a.
# At a fixed physical time t=π/2, M instead depends on a.
fractional_mean = mean_motion * (period / 4)
@assert coefficient_norm(fractional_mean - π / 2) < 1.0e-13
fractional, _ = kepler_iterations(fractional_mean, e; newton = true)
@assert coefficient_norm(differentiate(fractional, 1)) < 1.0e-12
mean = mean_motion * (π / 2)
newton, newton_errors = kepler_iterations(mean, e; newton = true)
fixed, fixed_errors = kepler_iterations(mean, e; newton = false)
@assert coefficient_norm(newton - fixed) < 2.0e-12
println("Quarter-period solution: ", fractional)
println("Fixed-time solution: ", newton)
println(
    (
        newton_iterations = length(newton_errors) - 1, fixed_point_iterations = length(fixed_errors) - 1,
        coefficient_difference = coefficient_norm(newton - fixed),
    )
)

#-
fig = Figure(size = (760, 420), fontsize = 15)
ax = Axis(fig[1, 1]; xlabel = "Iteration", ylabel = "Maximum coefficient residual", yscale = log10, title = "Convergence in the polynomial algebra")
scatterlines!(ax, 0:(length(newton_errors) - 1), max.(newton_errors, eps(Float64)); linewidth = 2, label = "Newton")
scatterlines!(ax, 0:(length(fixed_errors) - 1), max.(fixed_errors, eps(Float64)); linewidth = 2, label = "Fixed point")
axislegend(ax; position = :rt)
fig

# ## Fixed time and a fixed fraction of the period are different maps
# Evaluate both solutions in physical (a,e) coordinates. The residual panel
# substitutes the stored order-five map into the original scalar equation;
# its samples are a numerical check rather than a uniform error certificate.
axes_a, axes_e = range(0.9, 1.1; length = 65), range(0.1, 0.5; length = 65)
fixed_time = [newton([a - 1, e - 0.3]) for a in axes_a, e in axes_e]
period_fraction = [fractional([a - 1, e - 0.3]) for a in axes_a, e in axes_e]
residual = [abs(fixed_time[i, j] - e * sin(fixed_time[i, j]) - (π / 2) / (a * sqrt(a))) for (i, a) in enumerate(axes_a), (j, e) in enumerate(axes_e)]
parameter_fig = Figure(size = (1210, 420), fontsize = 14)
for (column, values, title, palette) in ((1, fixed_time, "Anomaly at fixed physical time", :viridis), (3, fixed_time - period_fraction, "Fixed time minus quarter-period", :balance), (5, log10.(max.(residual, eps(Float64))), "log₁₀ implicit-equation residual", :magma))
    ax = Axis(parameter_fig[1, column]; xlabel = "semimajor axis a", ylabel = "eccentricity e", title)
    heat = heatmap!(ax, axes_a, axes_e, values; colormap = palette)
    Colorbar(parameter_fig[1, column + 1], heat)
end
parameter_fig
