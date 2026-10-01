# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex06_gaussian_integral.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # Gaussian probability by polynomial integration
#
# Expand the standard normal density at zero, integrate its polynomial,
# and subtract the primitive at the endpoints. Compare with the independent
# error-function formula for the probability in [-1, 1].
using DifferentialAlgebra
using SpecialFunctions: erf
using CairoMakie

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
band!(ax, inside, zeros(length(inside)), density.(inside); color = (Makie.wong_colors()[1], 0.4))
fig
