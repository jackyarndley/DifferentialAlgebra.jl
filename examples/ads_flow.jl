# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # Domain splitting during orbit propagation
#
# Compare a final-time map with checkpointed and accepted-step flow splitting.
# Both return the same PiecewiseTaylorMap interface in initial-position
# coordinates. Vern9 chooses its internal ODE steps adaptively; ADS checks
# the initial state and the supplied checkpoint times.
# IntervalBound, polygon geometry and C0/C1/C2 blends are covered in the static
# ADS tutorials. These propagated maps use heuristic spatial indicators and
# an ordinary ODE solver; they do not certify time integration error.
using DifferentialAlgebra
using OrdinaryDiffEqVerner
using SciMLBase: successful_retcode, DiscreteCallback, terminate!
using CairoMakie

set_theme!(palette = (color = Makie.to_colormap(:tab20),))

function kepler_ode!(du, u, μ, time)
    factor = -μ * (u[1]^2 + u[2]^2)^(-3 / 2)
    du[1], du[2] = u[3], u[4]
    du[3], du[4] = factor * u[1], factor * u[2]
    return nothing
end

initial(x) = [x[1], x[2], zero(x[1]), 1.1 + zero(x[1])]
function advance(state, span)
    problem = ODEProblem(kepler_ode!, state, span, 1.0)
    solution = solve(problem, Vern9(); abstol = 1.0e-12, reltol = 1.0e-12, save_everystep = false)
    successful_retcode(solution) || error("Orbit integration failed")
    return solution.u[end]
end

lower, upper = [0.98, -0.02], [1.02, 0.02]
period = 2π / (2 - 1.1^2)^(3 / 2)
times = range(0, period; length = 9)
tolerance = 1.0e-5
terminal(x) = advance(initial(x), (0.0, period))
endpoint_map = adaptive_map(terminal, lower, upper; order = 4, atol = tolerance)
flow_map = adaptive_flow(advance, initial, lower, upper, times; order = 4, atol = tolerance)

# For online checks, the three-argument propagator calls the monitor after
# each accepted step, stops when requested and returns its state and time.
# A tuple selects this interface; a vector of times selects checkpoints.
function advance(state, span, monitor)
    condition(u, t, integrator) = monitor(u, t)
    callback = DiscreteCallback(condition, terminate!; save_positions = (false, true))
    problem = ODEProblem(kepler_ode!, state, span, 1.0)
    solution = solve(problem, Vern9(); callback, abstol = 1.0e-12, reltol = 1.0e-12, save_everystep = false)
    successful_retcode(solution) || error("Orbit integration failed")
    return (; state = solution.u[end], time = solution.t[end])
end
online_map = adaptive_flow(
    advance, initial, lower, upper, (0.0, period);
    order = 4, atol = tolerance, check_points = false
)
@assert endpoint_map.converged && flow_map.converged && online_map.converged

# An unsplit patch reuses its propagated state at the next interval, including
# guard coefficients. If a checkpoint fails, both children restart from their
# own initial conditions. This recovers information lost by the parent map.
# The checkpoint calculation includes independent numeric trajectories.
# Here online monitoring uses guard coefficients alone to avoid scalar solves
# at every accepted step. Both need independent validation of spatial accuracy.
points = [
    [x, y] for x in range(lower[1], upper[1]; length = 13),
        y in range(lower[2], upper[2]; length = 13)
]
reference = terminal.(points)
comparisons = (("Final-time", endpoint_map), ("Checkpoints", flow_map), ("Accepted steps", online_map))
validation_errors = [[maximum(abs, map(x) - truth) for (x, truth) in zip(points, reference)] for (_, map) in comparisons]
for ((name, map), errors) in zip(comparisons, validation_errors)
    @assert maximum(errors) < tolerance
    @assert sum(prod(p.upper - p.lower) for p in map.patches) ≈ prod(upper - lower)
    println((method = name, patches = length(map.patches), maximum_error = maximum(errors)))
end
println("Nominal orbit closure error: ", maximum(abs, flow_map([1.0, 0.0]) - initial([1.0, 0.0])))

# Each patch is indexed by its initial-position box, even though its map
# evaluates the final orbital state. The boundary images show nonlinear shear.
fig = Figure(size = (1250, 800), fontsize = 14)
for (column, (name, map)) in enumerate(comparisons)
    domain_axis = Axis(fig[1, column]; xlabel = "initial x", ylabel = "initial y", title = "$name: $(length(map.patches)) leaves", aspect = DataAspect())
    image_axis = Axis(fig[2, column]; xlabel = "final x", ylabel = "final y", title = "Mapped patch boundaries\nAxes scaled independently", aspect = 1)
    for (i, patch) in enumerate(map.patches)
        x0, y0 = patch.lower
        x1, y1 = patch.upper
        color = Makie.to_colormap(:tab20)[mod1(i, 20)]
        lines!(domain_axis, [x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0]; color, linewidth = 1)
        edge = range(0, 1; length = 21)
        boundary = vcat(
            [[x0 + t * (x1 - x0), y0] for t in edge], [[x1, y0 + t * (y1 - y0)] for t in edge],
            [[x1 - t * (x1 - x0), y1] for t in edge], [[x0, y1 - t * (y1 - y0)] for t in edge]
        )
        states = map.(boundary)
        lines!(image_axis, first.(states), getindex.(states, 2); color, linewidth = 1)
    end
end
fig

# Compare every method on the same independently propagated grid. Pointwise
# discrepancy includes numerical solver effects; it is not a uniform bound.
# The empirical CDF retains all grid errors instead of displaying only a maximum.
# Errors below machine epsilon sit at the log-scale display floor.
accuracy_fig = Figure(size = (1250, 410), fontsize = 14)
ax = Axis(accuracy_fig[1, 1]; xlabel = "maximum state-component error", ylabel = "fraction of validation points", xscale = log10, title = "Independent numeric trajectories")
for ((name, _), errors) in zip(comparisons, validation_errors)
    values = sort(vec(max.(errors, eps(Float64))))
    lines!(ax, values, collect(eachindex(values)) ./ length(values); label = name)
end
vlines!(ax, [tolerance]; color = :black, linestyle = :dash, label = "Spatial atol (heuristic)")
axislegend(ax; position = :rb, labelsize = 10)
ax = Axis(accuracy_fig[1, 2]; xticks = (1:3, collect(first.(comparisons))), ylabel = "leaf count", title = "Partition size, not integration cost")
barplot!(ax, 1:3, [length(m.patches) for (_, m) in comparisons]; color = Makie.to_colormap(:tab10)[1:3])
ax = Axis(accuracy_fig[1, 3]; xlabel = "initial x (initial y=0)", ylabel = "sampled maximum component error", yscale = log10, title = "Shared physical-coordinate section")
for ((name, _), errors) in zip(comparisons, validation_errors)
    section = reshape(errors, 13, 13)[:, 7]
    lines!(ax, range(lower[1], upper[1]; length = 13), max.(section, eps(Float64)); label = name)
end
axislegend(ax; position = :rt, labelsize = 10)
accuracy_fig
