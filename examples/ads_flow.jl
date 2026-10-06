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
    domain_axis = Axis(fig[1, column]; xlabel = "initial x", ylabel = "initial y", title = "$name: $(length(map.patches)) patches", aspect = DataAspect())
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
ax = Axis(accuracy_fig[1, 2]; xticks = (1:3, collect(first.(comparisons))), ylabel = "patch count", title = "Partition size, not integration cost")
barplot!(ax, 1:3, [length(m.patches) for (_, m) in comparisons]; color = Makie.to_colormap(:tab10)[1:3])
ax = Axis(accuracy_fig[1, 3]; xlabel = "initial x (initial y=0)", ylabel = "sampled maximum component error", yscale = log10, title = "Shared physical-coordinate section")
for ((name, _), errors) in zip(comparisons, validation_errors)
    section = reshape(errors, 13, 13)[:, 7]
    lines!(ax, range(lower[1], upper[1]; length = 13), max.(section, eps(Float64)); label = name)
end
axislegend(ax; position = :rt, labelsize = 10)
accuracy_fig

# ## Interval ADS on initial-epoch orbit diagnostics
# An ordinary ODE callback cannot acquire a validated time-error bound by
# changing its ADS estimator. Instead compare all four estimators on related
# explicit two-body diagnostics of the same initial-position box. These are
# initial energy, periapsis, apoapsis and period/(2π), with μ=1. The 1.1 literal
# denotes its stored binary value throughout, including validated calculations.
using IntervalArithmetic
function orbit_diagnostics(v)
    x, y = v
    speed = zero(x) + 1.1
    radius = sqrt(x^2 + y^2)
    energy = speed^2 / 2 - 1 / radius
    h = x * speed
    a = -1 / (2energy)
    e = sqrt(1 + 2energy * h^2)
    return [energy, a * (1 - e), a * (1 + e), a * sqrt(a)]
end
diagnostic_methods = (("GuardedTail", GuardedTail()), ("ExtrapolatedTail", ExtrapolatedTail()), ("LastTerms", LastTerms()), ("IntervalBound", IntervalBound()))
diagnostic_maps = [adaptive_map(orbit_diagnostics, lower, upper; order = 3, atol = 1 // 100000, estimator) for (_, estimator) in diagnostic_methods]
@assert all(m -> m.converged, diagnostic_maps)
@assert maximum(p -> maximum(sup.(abs.(p.error_bounds))), last(diagnostic_maps).patches) <= 1 / 100000
println("Static orbit diagnostic patch counts: ", [(label, length(m.patches)) for ((label, _), m) in zip(diagnostic_methods, diagnostic_maps)])

# Green intervals bound the original apoapsis expression on complete x cells
# at y=0. Every patch verifies the square-root and reciprocal assumptions on
# its whole cover. Numerical curves are illustrations; these certificates
# cover the static diagnostics and do not bound any propagated state above.
diagnostics_fig = Figure(size = (1320, 700), fontsize = 13)
diagnostic_edges = collect(range(lower[1], upper[1]; length = 65))
diagnostic_x = [x for k in 1:(length(diagnostic_edges) - 1) for x in (diagnostic_edges[k], diagnostic_edges[k + 1])]
for (column, ((label, estimator), m)) in enumerate(zip(diagnostic_methods, diagnostic_maps))
    certified = estimator isa IntervalBound
    color = certified ? :darkgreen : :dodgerblue
    ax = Axis(diagnostics_fig[1, column]; xlabel = "initial x", ylabel = "initial y", title = "$label: $(length(m.patches)) patches")
    for p in m.patches
        lo, hi = certified ? (inf.(domain(p)), sup.(domain(p))) : (p.lower, p.upper)
        poly!(ax, Rect2f(lo[1], lo[2], hi[1] - lo[1], hi[2] - lo[2]); color = (color, 0.1), strokecolor = color, strokewidth = 0.8)
    end
    ax = Axis(diagnostics_fig[2, column]; xlabel = "initial x (y=0)", ylabel = "initial osculating apoapsis", title = certified ? "Original-function cell enclosure" : "Ordinary polynomial fit")
    if certified
        cells = [enclose(m, [interval(diagnostic_edges[k], diagnostic_edges[k + 1]), interval(0)])[3] for k in 1:(length(diagnostic_edges) - 1)]
        @assert all(isguaranteed, cells)
        band!(ax, diagnostic_x, [inf(v) for v in cells for _ in 1:2], [sup(v) for v in cells for _ in 1:2]; color = (:green, 0.35))
    else
        lines!(ax, diagnostic_edges, [m([x, 0])[3] for x in diagnostic_edges]; color, linewidth = 3)
    end
    lines!(ax, diagnostic_edges, [orbit_diagnostics([x, 0])[3] for x in diagnostic_edges]; color = :black, linestyle = :dash)
end
diagnostics_fig
