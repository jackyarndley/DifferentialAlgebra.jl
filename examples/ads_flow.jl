# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # Domain splitting during orbit propagation
#
# Compare a final-time map with checkpointed and accepted-step flow splitting.
# Both return the same PiecewiseTaylorMap interface in initial-position
# coordinates. Vern9 chooses its internal ODE steps adaptively; ADS checks
# the initial state and the supplied checkpoint times.
using DifferentialAlgebra
using OrdinaryDiffEqVerner
using SciMLBase: successful_retcode, DiscreteCallback, terminate!
using CairoMakie

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
for (name, map) in (("Final-time splitting", endpoint_map), ("Checkpoint splitting", flow_map), ("Accepted-step splitting", online_map))
    errors = [maximum(abs, map(x) - truth) for (x, truth) in zip(points, reference)]
    @assert maximum(errors) < tolerance
    @assert sum(prod(p.upper - p.lower) for p in map.patches) ≈ prod(upper - lower)
    println((method = name, patches = length(map.patches), maximum_error = maximum(errors)))
end
println("Nominal orbit closure error: ", maximum(abs, flow_map([1.0, 0.0]) - initial([1.0, 0.0])))

# Each patch is indexed by its initial-position box, even though its map
# evaluates the final orbital state. The boundary images show nonlinear shear.
fig = Figure(size = (1080, 480), fontsize = 15)
domain = Axis(fig[1, 1]; xlabel = "Initial x", ylabel = "Initial y", title = "Accepted-step ADS partition", aspect = DataAspect())
image = Axis(fig[1, 2]; xlabel = "Final x", ylabel = "Final y", title = "Uncertainty after one orbit", aspect = 1)
for (i, patch) in enumerate(online_map.patches)
    x0, y0 = patch.lower
    x1, y1 = patch.upper
    color = Makie.wong_colors()[mod1(i, 7)]
    lines!(domain, [x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0]; color, linewidth = 1)
    edge = range(0, 1; length = 21)
    boundary = vcat(
        [[x0 + t * (x1 - x0), y0] for t in edge], [[x1, y0 + t * (y1 - y0)] for t in edge],
        [[x1 - t * (x1 - x0), y1] for t in edge], [[x0, y1 - t * (y1 - y0)] for t in edge]
    )
    states = online_map.(boundary)
    lines!(image, first.(states), getindex.(states, 2); color, linewidth = 1)
end
fig
