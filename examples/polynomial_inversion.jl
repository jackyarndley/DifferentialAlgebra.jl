# # Polynomial inversion
#
# Invert a polynomial map and check that its composition is the identity.
# Install the example environment as described in examples/README.md.

using DifferentialAlgebra

DifferentialAlgebra.init(10, 1)
x = DifferentialAlgebra.variable(1)
map = AlgebraicVector([sin(x)])
inverse = DifferentialAlgebra.invert(map)

# The inverse of sin near zero is asin, through the initialized order.
@assert DifferentialAlgebra.norm(inverse[1] - asin(x)) < 1e-14
composed = DifferentialAlgebra.evaluate(map[1], inverse)
@assert DifferentialAlgebra.norm(composed - x) < 1e-14
println("Inverse Taylor polynomial of sin(x):")
println(inverse[1])
