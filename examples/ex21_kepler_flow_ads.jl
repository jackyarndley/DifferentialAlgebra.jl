# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex21_kepler_flow_ads.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # Automatic domain splitting of an integrated Kepler flow
#
# Build a piecewise order-four map from uncertain initial position to orbital
# state. Every child callback integrates from its own initial box, so splitting
# recomputes missing information rather than rescaling an already truncated map.
# The ODE tolerances are tighter than the ADS target.
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

function flow(position, time)
    x, y = position
    initial = [x, y, zero(x), 1.1 + zero(x)]
    iszero(time) && return initial
    problem = ODEProblem(kepler_ode!, initial, (0.0, time), 1.0)
    solution = solve(problem, Vern9(); abstol = 1.0e-12, reltol = 1.0e-12, save_everystep = false)
    successful_retcode(solution) || error("Orbit integration failed")
    return solution.u[end]
end

lower, upper = [0.98, -0.02], [1.02, 0.02]
period = 2π / (2 - 1.1^2)^(3 / 2)
times = range(0, period; length = 5)
tolerance = 1.0e-5
maps = [adaptive_map(x -> flow(x, t), lower, upper; order = 4, atol = tolerance) for t in times]
@assert all(map -> map.converged, maps)
@assert length(maps[end].patches) > 1

# Check the initial coordinate transformation, domain edges, center and an
# independent grid. These also cover the upstream ADS validation scripts.
points = [
    [x, y] for x in range(lower[1], upper[1]; length = 11),
        y in range(lower[2], upper[2]; length = 11)
]
initial_error = maximum(maximum(abs, maps[1](x) - flow(x, 0.0)) for x in points)
errors = [maximum(abs, maps[end](x) - flow(x, period)) for x in points]
@assert initial_error < 1.0e-14
@assert maximum(errors) < tolerance
println("Patch counts at t/T = 0, ¼, ½, ¾, 1: ", length.(getproperty.(maps, :patches)))
println((initial_coordinate_error = initial_error, maximum_flow_error = maximum(errors)))
println("Nominal closure error: ", maximum(abs, maps[end]([1.0, 0.0]) - [1, 0, 0, 1.1]))

# Patch coordinates are physical. Compiled patch maps accept normalized local
# coordinates (position - center) ./ radius; the piecewise map does this for us.
patch = first(maps[end].patches)
@assert sum(prod(p.upper - p.lower) for p in maps[end].patches) ≈ prod(upper - lower)
@assert maps[end](patch.center) ≈ patch.map(zeros(2))
println(
    (
        lower = patch.lower, upper = patch.upper, center = patch.center,
        radius = patch.radius, outputs = noutputs(patch.map), status = patch.status,
    )
)
println("State at the first patch's center: ", maps[end](patch.center))
println("State at the physical lower corner: ", maps[end](lower))

# Visualize physical patch coordinates and the images of their boundaries.
# The image uses independent axis scales to make the narrow mapped set visible.
fig = Figure(size = (1080, 480), fontsize = 15)
domain = Axis(fig[1, 1]; xlabel = "Initial x", ylabel = "Initial y", title = "Initial-domain partition", aspect = DataAspect())
image = Axis(fig[1, 2]; xlabel = "Final x", ylabel = "Final y", title = "ADS uncertainty after one orbit", aspect = 1)
for (i, patch) in enumerate(maps[end].patches)
    x0, y0 = patch.lower
    x1, y1 = patch.upper
    color = Makie.wong_colors()[mod1(i, 7)]
    lines!(domain, [x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0]; color, linewidth = 1)
    edge = range(0, 1; length = 21)
    boundary = vcat(
        [[x0 + t * (x1 - x0), y0] for t in edge], [[x1, y0 + t * (y1 - y0)] for t in edge],
        [[x1 - t * (x1 - x0), y1] for t in edge], [[x0, y1 - t * (y1 - y0)] for t in edge]
    )
    states = maps[end].(boundary)
    lines!(image, first.(states), getindex.(states, 2); color, linewidth = 1)
end
fig
