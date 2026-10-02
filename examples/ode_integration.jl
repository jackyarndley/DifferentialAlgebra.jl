# # Orbit integration and state transition matrix
#
# Propagate a circular orbit and its second-order Taylor expansion through one
# revolution. All six initial position and velocity components are independent
# variables. Install the example environment as described in examples/README.md.

using DifferentialAlgebra
using OrdinaryDiffEqVerner
using SciMLBase: successful_retcode
using CairoMakie

set_theme!(palette = (color = Makie.to_colormap(:tab10),))

# Normalized Kepler equations: acceleration = -μ r / |r|³.
function kepler_ode!(du, u, μ, _)
    radius2 = u[1]^2 + u[2]^2 + u[3]^2
    factor = -μ / (radius2 * sqrt(radius2))
    for i in 1:3
        du[i] = u[i + 3]
        du[i + 3] = factor * u[i]
    end
    return nothing
end

μ = 1.0
initial = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
timespan = (0.0, 2π)
problem = ODEProblem(kepler_ode!, initial, timespan, μ)
nominal = solve(
    problem, Vern9(); abstol = 1.0e-12, reltol = 1.0e-12,
    saveat = range(timespan...; length = 161), save_everystep = false
)

perturbed = initial .+ variables(6; order = 2)

# Adaptive error control uses constant parts; it is not an error bound for every
# Taylor coefficient; use convergence checks when controlling higher-order terms.
solution = solve(
    remake(problem; u0 = perturbed), Vern9();
    abstol = 1.0e-12, reltol = 1.0e-12, save_everystep = false
)
final = solution.u[end]
comparison = solve(
    remake(problem; u0 = perturbed), Vern7();
    abstol = 1.0e-12, reltol = 1.0e-12, save_everystep = false
)
@assert successful_retcode(nominal) && successful_retcode(solution) && successful_retcode(comparison)
coefficient_error = maximum(coefficient_norm, final - comparison.u[end])
@assert coefficient_error < 1.0e-6
println("Maximum coefficient difference, Vern9 versus Vern7: ", coefficient_error)
constants = constant_term.(final)
@assert maximum(abs.(constants - nominal.u[end])) < 1.0e-9
@assert maximum(abs.(constants - initial)) < 1.0e-9

# The Jacobian at zero perturbation is the state transition matrix.
stm = constant_term.(jacobian(final))
@assert size(stm) == (6, 6)
@assert all(isfinite, stm)
println("Maximum nominal orbit error: ", maximum(abs.(constants - initial)))
println("State transition matrix:")
show(stdout, MIME"text/plain"(), stm); println()

# ## Linear and nonlinear sensitivity
#
# A first-order prediction uses the STM: Φ δx₀. Evaluating the second-order
# flow map also includes quadratic terms. Compare both with independent numeric
# integrations after small changes to the initial radial position.

displacements = 10.0 .^ range(-4, -2; length = 17)
linear_errors, quadratic_errors = Float64[], Float64[]
for displacement in displacements
    delta = [displacement, 0.0, 0.0, 0.0, 0.0, 0.0]
    reference_solution = solve(
        remake(problem; u0 = initial + delta), Vern9();
        abstol = 1.0e-12, reltol = 1.0e-12, save_everystep = false
    )
    @assert successful_retcode(reference_solution)
    reference = reference_solution.u[end]
    push!(linear_errors, maximum(abs, constants + stm * delta - reference))
    push!(quadratic_errors, maximum(abs, evaluate(final, delta) - reference))
end
@assert all(quadratic_errors .< linear_errors)
@assert quadratic_errors[end] < linear_errors[end] / 10

# The STM heatmap includes position and velocity components in the normalized
# units of this problem. The error plot shows the benefit of quadratic terms
# without needing to reintegrate the Taylor map for each initial condition.

fig = Figure(size = (1350, 430), fontsize = 15)
orbit_axis = Axis(fig[1, 1]; xlabel = "x", ylabel = "y", title = "One Kepler revolution", aspect = DataAspect())
lines!(orbit_axis, getindex.(nominal.u, 1), getindex.(nominal.u, 2); color = Makie.to_colormap(:tab10)[1], linewidth = 3)
scatter!(orbit_axis, [0.0], [0.0]; color = :black, markersize = 14, label = "Central body")
scatter!(orbit_axis, [initial[1]], [initial[2]]; color = Makie.to_colormap(:tab10)[2], markersize = 12, label = "Initial / final state")
axislegend(orbit_axis; position = :lb, labelsize = 11)
labels = ["x", "y", "z", "vx", "vy", "vz"]
stm_axis = Axis(fig[1, 2]; xlabel = "Initial component", ylabel = "Final component", title = "State transition matrix", xticks = (1:6, labels), yticks = (1:6, labels), yreversed = true, aspect = DataAspect())
limit = maximum(abs, stm)
heat = heatmap!(stm_axis, 1:6, 1:6, transpose(stm); colormap = :balance, colorrange = (-limit, limit))
Colorbar(fig[1, 3], heat)
error_axis = Axis(fig[1, 4]; xlabel = "Initial radial displacement", ylabel = "Maximum final-state error", title = "Nonlinear terms improve the map", xscale = log10, yscale = log10, xticks = ([1.0e-4, 1.0e-3, 1.0e-2], ["10⁻⁴", "10⁻³", "10⁻²"]), yticks = ([1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2], ["10⁻⁸", "10⁻⁶", "10⁻⁴", "10⁻²"]))
scatterlines!(error_axis, displacements, linear_errors; color = Makie.to_colormap(:tab10)[1], linewidth = 2, label = "STM (first order)")
scatterlines!(error_axis, displacements, quadratic_errors; color = Makie.to_colormap(:tab10)[2], linewidth = 2, label = "Second-order map")
axislegend(error_axis; position = :lt, labelsize = 11)

fig
#
