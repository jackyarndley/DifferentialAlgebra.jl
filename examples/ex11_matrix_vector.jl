# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex11_matrix_vector.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # Polynomial matrices and vectors
#
# Ordinary Julia arrays support polynomial entries. No special vector or
# matrix wrapper is needed. Matrix products use truncated scalar arithmetic.
using DifferentialAlgebra

x, y = variables((:x, :y); order = 3)
A = [1 + 0.2x y; -y 1 - 0.1x]
B = [x y; x * y y^2 + 0.5]
v = [x + y, x - y]
C = A * B
u = C * v
println("A B = ", C)
println("(A B) v = ", u)
println("Evaluation at (0.1, -0.2): ", evaluate(u, [0.1, -0.2]))

# Associativity holds in the truncated algebra. Evaluating a product at a
# finite point need not equal multiplying numeric evaluations: terms above
# degree three have already been discarded.
@assert maximum(coefficient_norm, u - A * (B * v)) < 1.0e-14
println("Associativity coefficient residual: ", maximum(coefficient_norm, u - A * (B * v)))
