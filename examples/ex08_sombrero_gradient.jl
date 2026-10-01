# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex08_sombrero_gradient.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # Gradient of the sombrero function
#
# Tutorial 1, exercise 8 from DACE.jl: expand sin(r)/r about (2, 3) and
# extract its gradient using first-order differential algebra.

using DifferentialAlgebra
using LinearAlgebra
using CairoMakie

function sombrero(x)
    r = sqrt(x[1]^2 + x[2]^2)
    return iszero(r) ? one(r) : sin(r) / r
end

dx, dy = variables((:δx, :δy); order = 1)
x = [2.0 + dx, 3.0 + dy]
z = sombrero(x)
grad_z = gradient(z)

println("Coordinates about (2, 3): ", x)
println("Sombrero function: ", z)
println("Gradient: ", grad_z)

# Check against the analytic gradient at the expansion center.
r = sqrt(13.0)
expected = (r * cos(r) - sin(r)) / r^3 .* [2.0, 3.0]
@assert constant_term(grad_z) ≈ expected

# ## The gradient defines a tangent approximation
#
# Follow the gradient direction through the expansion center. The first-order
# Taylor map gives the tangent line to this section of the surface. It matches
# both the function value and directional derivative at the center.

center = [2.0, 3.0]
direction = expected / norm(expected)
distance = range(-1.5, 1.5; length = 201)
section = [sombrero(center + t * direction) for t in distance]
tangent = [z(t * direction) for t in distance]
@assert z([0.0, 0.0]) ≈ sombrero(center)
@assert dot(constant_term(grad_z), direction) ≈ norm(expected)

grid = range(-5, 5; length = 151)
height = [sombrero([a, b]) for a in grid, b in grid]
fig = Figure(size = (1100, 440), fontsize = 15)
surface_axis = Axis(fig[1, 1]; xlabel = "x", ylabel = "y", title = "Sombrero function", aspect = DataAspect())
contours = contourf!(surface_axis, grid, grid, height; levels = 24, colormap = :viridis)
Colorbar(fig[1, 2], contours; label = "sin(r) / r")
path = [center + t * direction for t in distance]
lines!(surface_axis, first.(path), last.(path); color = :white, linewidth = 3, label = "Gradient direction")
scatter!(surface_axis, [center[1]], [center[2]]; color = :white, strokecolor = :black, strokewidth = 1, markersize = 12, label = "Expansion center")
axislegend(surface_axis; position = :lb, labelsize = 11, backgroundcolor = (:black, 0.65), labelcolor = :white)
section_axis = Axis(fig[1, 3]; xlabel = "Distance along gradient", ylabel = "Function value", title = "First-order tangent map")
lines!(section_axis, distance, section; color = :black, linewidth = 3, label = "Exact section")
lines!(section_axis, distance, tangent; color = Makie.wong_colors()[1], linewidth = 2, linestyle = :dash, label = "Taylor map")
scatter!(section_axis, [0.0], [sombrero(center)]; color = :black, markersize = 10)
axislegend(section_axis; position = :lt)

fig
#
