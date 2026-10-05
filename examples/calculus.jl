# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # Differentiation and integration
#
# Differentiation loses one known order. Integration chooses a zero constant
# and discards terms beyond the active order. Compare like orders and constants.
using DifferentialAlgebra
using CairoMakie

x, = variables((:x,); order = 20)
p = sin(x)
derivative = differentiate(p, 1)
primitive = integrate(p, 1)
expected_derivative = with_order(() -> cos(x), 19)
@assert coefficient_norm(derivative - expected_derivative) < 1.0e-15
@assert coefficient_norm(primitive - (1 - cos(x))) < 1.0e-15
println("d sin(x) / dx = ", derivative)
println("∫ sin(x) dx, with value zero at x=0: ", primitive)
println("Derivative residual through order 19: ", coefficient_norm(derivative - expected_derivative))
println("Primitive residual: ", coefficient_norm(primitive - (1 - cos(x))))

# Differentiation/integration act on the stored finite polynomial. Their
# agreement with the original function away from zero is a separate question.
points = range(-3, 3; length = 301)
calculus_fig = Figure(size = (1050, 650), fontsize = 14)
for (column, computed, exact, label) in ((1, derivative, cos, "Derivative"), (2, primitive, t -> 1 - cos(t), "Primitive, zero constant"))
    ax = Axis(calculus_fig[1, column]; xlabel = "x", ylabel = "value", title = label)
    lines!(ax, points, exact.(points); color = :black, linewidth = 3, label = "Original-function formula")
    lines!(ax, points, computed.(points); linestyle = :dash, label = "Polynomial calculus")
    axislegend(ax; position = :lb, labelsize = 11)
    ax = Axis(calculus_fig[2, column]; xlabel = "x", ylabel = "sampled absolute error", yscale = log10)
    lines!(ax, points, max.(abs.(computed.(points) - exact.(points)), eps(Float64)))
end
calculus_fig

# ## Gaussian probability by polynomial integration
#
# Expand the standard normal density at zero, integrate its polynomial,
# and subtract the primitive at the endpoints. Compare with the independent
# error-function formula for the probability in [-1, 1].
using SpecialFunctions: erf
using CairoMakie

set_theme!(palette = (color = Makie.to_colormap(:tab10),))

x, = variables((:x,); order = 24)
density = exp(-x^2 / 2) / sqrt(2π)
primitive = integrate(density, 1)
probability = primitive(1.0) - primitive(-1.0)
reference = erf(1 / sqrt(2))
@assert abs(probability - reference) < 1.0e-13
println("Density polynomial: ", density)
println((probability = probability, reference = reference, error = abs(probability - reference)))

# The shaded area is the integrated probability.
points = range(-3, 3; length = 301)
inside = range(-1, 1; length = 101)
fig = Figure(size = (1050, 400), fontsize = 14)
ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "Probability density", title = "Standard normal probability in [-1, 1]")
lines!(ax, points, exp.(-points .^ 2 ./ 2) ./ sqrt(2π); color = :black, linewidth = 2)
band!(ax, inside, zeros(length(inside)), density.(inside); color = (Makie.to_colormap(:tab10)[1], 0.4))
ax = Axis(fig[1, 2]; xlabel = "symmetric limit a", ylabel = "P(-a ≤ X ≤ a)", title = "Integrated polynomial probability")
limits = range(0, 2.5; length = 151)
lines!(ax, limits, erf.(limits ./ sqrt(2)); color = :black, linewidth = 3, label = "Error-function formula")
lines!(ax, limits, [primitive(a) - primitive(-a) for a in limits]; linestyle = :dash, label = "Integrated density polynomial")
scatter!(ax, [1], [probability]; color = :dodgerblue, label = "Tested interval [-1,1]")
axislegend(ax; position = :rb, labelsize = 11)
fig
