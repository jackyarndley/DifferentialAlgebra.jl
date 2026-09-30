# # Sine function
#
# Construct a Taylor expansion, inspect its coefficients and evaluate it.
# Install the example environment as described in examples/README.md.

using DifferentialAlgebra

# Initialize a 20th-order algebra with one variable.
DifferentialAlgebra.init(20, 1)
x = DifferentialAlgebra.variable(1)
p = sin(x)

# Coefficients multiply ordinary monomials, so the cubic coefficient is -1/6.
@assert DifferentialAlgebra.getCoefficient(p, [3]) ≈ -1/6
value = DifferentialAlgebra.evaluate(p, [1.0])
@assert isapprox(value, sin(1.0); atol=1e-14)
println("Taylor approximation of sin(1): ", value)
