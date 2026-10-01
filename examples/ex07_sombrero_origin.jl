# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex07_sombrero_origin.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # The sombrero function at its removable singularity
#
# sin(r)/r is smooth at r=0 although evaluating its quotient there fails.
# With q=x²+y², use the entire series Σ (-q)ᵏ/(2k+1)!.
# Each q raises the degree by two; retain every term through order ten.
using DifferentialAlgebra
using CairoMakie

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
