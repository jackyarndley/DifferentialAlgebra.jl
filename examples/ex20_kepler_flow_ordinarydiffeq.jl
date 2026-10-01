# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex20_kepler_flow_ordinarydiffeq.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # Kepler flow with OrdinaryDiffEq
#
# Propagate the same eccentric orbit and uncertainty box as the adaptive
# Runge–Kutta tutorial, using the maintained Vern9 and Vern7 solvers.
# Compare every retained coefficient and test finite perturbations against
# independent Float64 integrations.
using DifferentialAlgebra
using OrdinaryDiffEqVerner
using SciMLBase: successful_retcode
using CairoMakie

function kepler_ode!(du, u, μ, time)
    factor = -μ * (u[1]^2 + u[2]^2)^(-3 / 2)
    du[1], du[2] = u[3], u[4]
    du[3], du[4] = factor * u[1], factor * u[2]
    return nothing
end

δx, δy = variables((:δx, :δy); order = 3)
initial = [1 + 0.01δx, 0.01δy, zero(δx), 1.1 + zero(δx)]
period = 2π / (2 - 1.1^2)^(3 / 2)
times = range(0, period; length = 9)
problem = ODEProblem(kepler_ode!, initial, (0.0, period), 1.0)
solution = solve(problem, Vern9(); abstol = 1.0e-12, reltol = 1.0e-12, saveat = times)
comparison = solve(problem, Vern7(); abstol = 1.0e-12, reltol = 1.0e-12, save_everystep = false)
@assert successful_retcode(solution) && successful_retcode(comparison)
coefficient_errors = coefficient_norm.(solution.u[end] - comparison.u[end])
@assert maximum(coefficient_errors) < 1.0e-8
println("Final Vern9 Taylor map: ", solution.u[end])
println("Maximum coefficient difference per component (Vern9 versus Vern7): ", coefficient_errors)

# The solver's default norm controls constant parts. The independent solver
# and pointwise checks above and below also inspect the uncertainty terms.
nominal_problem = remake(problem; u0 = [1.0, 0.0, 0.0, 1.1])
nominal = solve(
    nominal_problem, Vern9(); abstol = 1.0e-12, reltol = 1.0e-12,
    saveat = range(0, period; length = 241)
)
@assert successful_retcode(nominal)
angles = range(0, 2π; length = 65)
boundary = [[cos(θ), sin(θ)] for θ in angles]
errors, half_width_errors = Float64[], Float64[]
for u in boundary
    reference = solve(
        remake(nominal_problem; u0 = [1 + 0.01u[1], 0.01u[2], 0.0, 1.1]),
        Vern9(); abstol = 1.0e-12, reltol = 1.0e-12, save_everystep = false
    )
    @assert successful_retcode(reference)
    push!(errors, maximum(abs, evaluate(solution.u[end], u) - reference.u[end]))
    narrower = solve(
        remake(nominal_problem; u0 = [1 + 0.005u[1], 0.005u[2], 0.0, 1.1]),
        Vern9(); abstol = 1.0e-12, reltol = 1.0e-12, save_everystep = false
    )
    @assert successful_retcode(narrower)
    push!(half_width_errors, maximum(abs, evaluate(solution.u[end], u / 2) - narrower.u[end]))
end
@assert all(isfinite, errors)
@assert maximum(half_width_errors) < maximum(errors) / 10
@assert maximum(abs, nominal.u[end] - nominal.u[1]) < 1.0e-8
println("Maximum finite-perturbation error: ", maximum(errors))
println("Maximum error after halving the uncertainty radius: ", maximum(half_width_errors))
println("Nominal orbit closure error: ", maximum(abs, nominal.u[end] - nominal.u[1]))

# A cubic map has leading fourth-order truncation error. Halving the initial
# uncertainty should reduce that error by about a factor of sixteen.
#-
fig = Figure(size = (1100, 470), fontsize = 15)
ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "y", title = "Flow of an uncertainty circle", aspect = DataAspect())
lines!(ax, first.(nominal.u), getindex.(nominal.u, 2); color = :black, linewidth = 2)
for (i, map) in enumerate(solution.u)
    states = evaluate.(Ref(compile(map)), boundary)
    lines!(ax, first.(states), getindex.(states, 2); color = (i - 1) / (length(times) - 1), colorrange = (0, 1), colormap = :viridis, linewidth = 2)
end
scatter!(ax, [0.0], [0.0]; color = :orange, markersize = 12)
Colorbar(fig[2, 1]; limits = (0, 1), colormap = :viridis, vertical = false, label = "Time / nominal period")
err = Axis(fig[1, 2]; xlabel = "Initial uncertainty angle (rad)", ylabel = "Maximum state error", yscale = log10, title = "Taylor truncation at finite perturbations")
lines!(err, angles, max.(errors, eps(Float64)); linewidth = 2, label = "Radius 0.01")
lines!(err, angles, max.(half_width_errors, eps(Float64)); linewidth = 2, label = "Radius 0.005")
axislegend(err; position = :lb)
fig
