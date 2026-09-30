# # Gradient of a radial function
#
# Polynomial differentiation gives the gradient of sin(r)/r at (2, 3).
# Install the example environment as described in examples/README.md.

using DifferentialAlgebra

dx, dy = variables(2; order=1)
x, y = 2.0 + dx, 3.0 + dy
r = sqrt(x*x + y*y)
p = sin(r)/r

g = constant_term.(gradient(p))
r0 = sqrt(13.0)
expected = (r0*cos(r0) - sin(r0))/r0^3 .* [2.0, 3.0]
@assert isapprox(g, expected; atol=1e-14)
println("Gradient at (2, 3): ", g)
