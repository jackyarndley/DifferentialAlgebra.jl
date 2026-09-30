# # Polynomial inversion
#
# Invert a polynomial map and check that its composition is the identity.
# Install the example environment as described in examples/README.md.

using DifferentialAlgebra

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
