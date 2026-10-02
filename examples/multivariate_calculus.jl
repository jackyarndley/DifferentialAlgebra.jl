# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # Multivariate derivatives and removable singularities
#
# sin(r)/r is smooth at r=0 although evaluating its quotient there fails.
# With q=x²+y², use the entire series Σ (-q)ᵏ/(2k+1)!.
# Each q raises the degree by two; retain every term through order ten.
using DifferentialAlgebra
using CairoMakie

set_theme!(palette = (color = Makie.to_colormap(:tab10),))

x, y = variables((:x, :y); order = 10)
q = x^2 + y^2
p = one(q)
term = one(q)
for k in 1:5
    global term *= -q / ((2k) * (2k + 1))
    global p += term
end
@assert coefficient(p, [10, 0]) ≈ -1 / factorial(11)
@assert p([0.0, 0.0]) == 1
println("Sombrero expansion at the origin: ", p)

# A radial slice checks the removable singularity and truncation error.
radius = range(-2, 2; length = 201)
reference = sinc.(radius ./ π)
approximation = [p([r, 0.0]) for r in radius]
error = maximum(abs, approximation - reference)
@assert error < 7.0e-7
println("Maximum radial error on [-2, 2]: ", error)

#-
fig = Figure(size = (950, 390), fontsize = 15)
ax = Axis(fig[1, 1]; xlabel = "Radius", ylabel = "sin(r)/r", title = "Smooth through the origin")
lines!(ax, radius, reference; color = :black, linewidth = 3, label = "Exact")
lines!(ax, radius, approximation; linestyle = :dash, linewidth = 2, label = "Order 10")
axislegend(ax; position = :lb)
err = Axis(fig[1, 2]; xlabel = "Radius", ylabel = "Absolute error", yscale = log10)
lines!(err, radius, max.(abs.(approximation - reference), eps(Float64)); linewidth = 2)
fig

# ## Gradient of the sombrero function
#
# Expand sin(r)/r about (2, 3) and
# extract its gradient using first-order differential algebra.

using LinearAlgebra

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
lines!(section_axis, distance, tangent; color = Makie.to_colormap(:tab10)[1], linewidth = 2, linestyle = :dash, label = "Taylor map")
scatter!(section_axis, [0.0], [sombrero(center)]; color = :black, markersize = 10)
axislegend(section_axis; position = :lt)

fig
#
