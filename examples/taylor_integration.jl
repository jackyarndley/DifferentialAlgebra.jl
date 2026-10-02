# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # Taylor integration in time
#
# Solve an ODE by expanding its solution about the current epoch, evaluating
# the expansion at a small time offset, and repeating at the new epoch.
# Unlike a Runge–Kutta step, each step constructs a local polynomial in time.
# Its coefficients can themselves be polynomials in uncertain initial data.
using DifferentialAlgebra
using OrdinaryDiffEqVerner
using SciMLBase: successful_retcode
using CairoMakie

set_theme!(palette = (color = Makie.to_colormap(:tab10),))

# ## Coefficient recurrence
#
# Substitute ``u(t₀+τ)=\sum_{k=0}^p u_k τ^k`` into ``u'=f(u,t)``.
# Matching coefficients gives ``u_{k+1}=[τ^k]f(u,t₀+τ)/(k+1)``.
# This is coefficient-by-coefficient Picard iteration: degree k of the right
# hand side only depends on the state coefficients through degree k.
# [`taylor_expand`](@ref) implements the recurrence, without finite differences
# or symbolic differentiation of the right-hand side.
exponential = only(taylor_expand((u, p, t) -> u, [1.0], 0.0; order = 12))
println("First five coefficients of exp(τ): ", exponential.coefficients[1:5])
println("Taylor step at τ = 0.2: ", exponential(0.2))
@assert isapprox(exponential(0.2), exp(0.2); rtol = 1.0e-14)

# ## Adaptive steps through SciML
#
# Loading OrdinaryDiffEq enables `TaylorMethod`, so time expansions can use
# SciML's `ODEProblem`, `solve`, interpolation and callback interfaces.
# The last two retained time coefficients supply a local error estimate.
# Polynomial coefficients use their sum of absolute values, which bounds
# their magnitude on the normalized uncertainty box. The SciML controller
# scales this estimate by `abstol + reltol*size(u)`, accepts or rejects the
# step, and proposes the next step size. This is a truncation estimate,
# not a rigorous remainder bound; independent trajectories below check it.
#
# ## An eccentric Kepler orbit
#
# Use normalized units with μ=1 and nominal semimajor axis a=1. The orbital
# period is 2π. The out-of-place right-hand side is shared by Taylor integration
# and an independent Vern9 integration; it accepts numbers or time series.
function kepler(u, μ, t)
    factor = -μ / (u[1]^2 + u[2]^2)^(3 // 2)
    return [u[3], u[4], factor * u[1], factor * u[2]]
end
periapsis(a, e) = [a * (1 - e), zero(a), zero(a), sqrt((1 + e) / (a * (1 - e)))]
initial = periapsis(1.0, 0.4)
period = 2π
problem = ODEProblem(kepler, initial, (0.0, period), 1.0)
trajectory = solve(problem, TaylorMethod(18); abstol = 1.0e-13, reltol = 1.0e-13)
@assert successful_retcode(trajectory)
reference = solve(problem, Vern9(); abstol = 1.0e-13, reltol = 1.0e-13)
@assert successful_retcode(reference)
errors = [maximum(abs, u - reference(t)) for (t, u) in zip(trajectory.t, trajectory.u)]
closure_error = maximum(abs, last(trajectory.u) - initial)
@assert maximum(errors) < 2.0e-10
@assert closure_error < 2.0e-10
reverse_trajectory = solve(
    remake(problem; u0 = last(trajectory.u), tspan = (period, 0.0)),
    TaylorMethod(18); abstol = 1.0e-13, reltol = 1.0e-13
)
@assert successful_retcode(reverse_trajectory)
@assert maximum(abs, last(reverse_trajectory.u) - initial) < 4.0e-10
println("Numeric Taylor steps per orbit: ", length(trajectory.t) - 1)
println("Maximum difference from Vern9: ", maximum(errors))
println("One-period closure error: ", closure_error)

# ## Independent time and uncertainty orders
#
# Keep degree two in (δa, δe), but degree eighteen in time. A TimeSeries stores
# nineteen coefficients, each with at most six uncertainty monomials. Terms
# such as ``τ^{18}δa^2`` survive: time degree is not part of the uncertainty
# degree. No additional variables or high-dimensional monomial tables are
# needed. The multivariate algebra itself still uses a total-degree limit.
δa, δe = variables((:δa, :δe); order = 2)
map_initial = periapsis(1 + 2.0e-4δa, 0.4 + 1.0e-4δe)
flow = solve(remake(problem; u0 = map_initial), TaylorMethod(18); abstol = 1.0e-13, reltol = 1.0e-13)
@assert successful_retcode(flow)
final_map = last(flow.u)
@assert max_order() == 2
println("Uncertainty order: ", max_order())
println("Time order: ", 18)
println("Polynomial Taylor steps per orbit: ", length(flow.t) - 1)
println("Final x map: ", final_map[1])
println("Final y map: ", final_map[2])

# Validate across the uncertainty box with separate numeric integrations.
# These errors include both the time integration and the degree-two map error.
grid = range(-1, 1; length = 9)
map_errors = [
    let
        u0 = periapsis(1 + 2.0e-4a, 0.4 + 1.0e-4e)
        numeric = solve(
            ODEProblem(kepler, u0, (0.0, period), 1.0), Vern9();
            abstol = 1.0e-13, reltol = 1.0e-13, save_everystep = false
        )
        @assert successful_retcode(numeric)
        maximum(abs, evaluate(final_map, [a, e]) - last(numeric.u))
    end for a in grid, e in grid
]
@assert maximum(map_errors) < 2.0e-7
println("Maximum sampled uncertainty-map error: ", maximum(map_errors))

#-
fig = Figure(size = (960, 740), fontsize = 15)
orbit_axis = Axis(fig[1, 1]; xlabel = "x", ylabel = "y", aspect = DataAspect(), title = "Taylor integration: eccentric Kepler orbit")
orbit = trajectory(range(0, period; length = 501))
lines!(orbit_axis, orbit[1, :], orbit[2, :]; linewidth = 2)
scatter!(orbit_axis, [0.0], [0.0]; color = :black, markersize = 10)
step_axis = Axis(fig[1, 2]; xlabel = "Time / period", ylabel = "Step size", title = "Adaptive time steps")
stairs!(step_axis, trajectory.t[1:(end - 1)] ./ period, diff(trajectory.t); linewidth = 3, label = "Numeric state")
stairs!(step_axis, flow.t[1:(end - 1)] ./ period, diff(flow.t); linewidth = 2, linestyle = :dash, label = "Degree-two map")
axislegend(step_axis; position = :ct)
error_axis = Axis(fig[2, 1]; xlabel = "Time / period", ylabel = "Maximum state difference", yscale = log10, title = "Independent Vern9 comparison")
lines!(error_axis, trajectory.t ./ period, max.(errors, eps(Float64)); linewidth = 2)
map_axis = Axis(fig[2, 2]; xlabel = "Normalized δa", ylabel = "Normalized δe", title = "Uncertainty-map error")
heat = heatmap!(map_axis, grid, grid, log10.(max.(map_errors, 1.0e-15)); colormap = :viridis)
Colorbar(fig[2, 3], heat; label = "log₁₀ maximum state error")
fig
