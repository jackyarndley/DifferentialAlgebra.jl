# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex09_sinx_inversion.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # Polynomial inversion
#
# Invert a polynomial map and check that its composition is the identity.
# Install the example environment as described in examples/README.md.

using DifferentialAlgebra
using CairoMakie

x, = variables((:x,); order = 10)
map = [sin(x)]
inverse = invert(map)

# The inverse of sin near zero is asin, through the initialized order.
@assert DifferentialAlgebra.coefficient_norm(inverse[1] - asin(x)) < 1.0e-14
composed = evaluate(map[1], inverse)
@assert DifferentialAlgebra.coefficient_norm(composed - x) < 1.0e-14
println("Inverse Taylor polynomial of sin(x):")
println(inverse[1])
println("asin(x):")
println(asin(x))

# ## A local inverse
#
# Polynomial composition is exact up to the truncation order, but evaluating
# the inverse away from its expansion center still incurs truncation error.
# The derivatives of asin diverge at ±1, so a fixed-order polynomial loses
# accuracy as its argument approaches those endpoints.

points = range(-0.95, 0.95; length = 401)
exact = asin.(points)
orders = (1, 3, 5, 9)
approximations = [with_order(() -> invert([sin(x)])[1], order).(points) for order in orders]
errors = [abs.(values - exact) for values in approximations]
@assert maximum(errors[end]) < maximum(errors[1])
@assert isapprox(inverse[1](0.1), asin(0.1); atol = 3.0e-13)

fig = Figure(size = (1050, 420), fontsize = 15)
value_axis = Axis(fig[1, 1]; xlabel = "y", ylabel = "Inverse value", title = "Inverting sin(x) near zero")
error_axis = Axis(fig[1, 2]; xlabel = "y", ylabel = "Absolute error", yscale = log10, title = "The inverse is local")
lines!(value_axis, points, exact; color = :black, linewidth = 3, label = "asin(y)")
colors = Makie.wong_colors()
for (i, order) in enumerate(orders)
    color = colors[i]
    lines!(value_axis, points, approximations[i]; color, linewidth = 2, label = "Order $order")
    lines!(error_axis, points, max.(errors[i], eps(Float64)); color, linewidth = 2)
end
axislegend(value_axis; position = :lt, labelsize = 12)

fig
#
