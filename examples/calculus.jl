# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # Differentiation and integration
#
# Differentiation loses one known order. Integration chooses a zero constant
# and discards terms beyond the active order. Compare like orders and constants.
using DifferentialAlgebra

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
fig = Figure(size = (760, 400), fontsize = 15)
ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "Probability density", title = "Standard normal probability in [-1, 1]")
lines!(ax, points, exp.(-points .^ 2 ./ 2) ./ sqrt(2π); color = :black, linewidth = 2)
band!(ax, inside, zeros(length(inside)), density.(inside); color = (Makie.to_colormap(:tab10)[1], 0.4))
fig
