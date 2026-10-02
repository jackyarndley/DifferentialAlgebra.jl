# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # Taylor expansions and truncation
#
# Construct a Taylor expansion, inspect its coefficients and evaluate it.
# Install the example environment as described in examples/README.md.

using DifferentialAlgebra
using CairoMakie

set_theme!(palette = (color = Makie.to_colormap(:tab10),))

# Initialize a 20th-order algebra with one variable.
x, = variables((:x,); order = 20)
p = sin(x)

# Coefficients multiply ordinary monomials, so the cubic coefficient is -1/6.
@assert coefficient(p, [3]) ≈ -1 / 6
value = p(1.0)
@assert isapprox(value, sin(1.0); atol = 1.0e-14)
println("Taylor approximation of sin(1): ", value)
println("x = ", x)
println("sin(x) = ", p)
println("Base.sin(1) = ", sin(1.0))

# ## Arithmetic in a truncated algebra
#
# Trigonometric identities hold coefficient by coefficient up to rounding.
# A polynomial whose minimum degree exceeds the working order is zero.
identity_residual = coefficient_norm(sin(x)^2 + cos(x)^2 - 1)
nilpotent = (cos(x) - 1)^11
@assert identity_residual < 1.0e-14
@assert iszero(nilpotent)
println("Identity coefficient residual: ", identity_residual)
println("(cos(x) - 1)¹¹ at order 20: ", nilpotent)

# ## Approximation across an interval
#
# The same expansion can be recomputed at a lower working order without changing
# the algebra. Increasing the order extends the useful interval about zero.
# Sine has only odd powers, so orders 19 and 20 give the same polynomial.

points = range(-π, π; length = 401)
exact = sin.(points)
orders = (3, 7, 11, 19)
approximations = [with_order(() -> sin(x), order).(points) for order in orders]
errors = [abs.(values - exact) for values in approximations]
@assert maximum(errors[end]) < 6.0e-10
@assert maximum(errors[end]) < maximum(errors[1])

# The right panel uses a logarithmic scale. Errors below machine epsilon are
# placed at the plotting floor; this is a display convention, not an error bound.
fig = Figure(size = (1050, 420), fontsize = 15)
value_axis = Axis(fig[1, 1]; xlabel = "x", ylabel = "sin(x)", title = "Taylor expansions about zero")
error_axis = Axis(fig[1, 2]; xlabel = "x", ylabel = "Absolute error", yscale = log10, title = "Accuracy improves with order")
lines!(value_axis, points, exact; color = :black, linewidth = 3, label = "sin(x)")
colors = Makie.to_colormap(:tab10)
for (i, order) in enumerate(orders)
    color = colors[i]
    lines!(value_axis, points, approximations[i]; color, linewidth = 2, label = "Order $order")
    lines!(error_axis, points, max.(errors[i], eps(Float64)); color, linewidth = 2)
end
axislegend(value_axis; position = :lt, labelsize = 12)

fig
#
