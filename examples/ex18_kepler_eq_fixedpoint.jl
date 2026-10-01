# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex18_kepler_eq_fixedpoint.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # Newton lifting and fixed-point iteration
#
# Both algorithms solve Kepler's implicit equation in a truncated algebra.
# Newton iteration doubles the number of correct orders near a scalar root;
# fixed-point iteration E ← M + e sin(E) converges linearly here.
# Monitor every polynomial coefficient, not only the nominal residual.
using DifferentialAlgebra
using CairoMakie

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
