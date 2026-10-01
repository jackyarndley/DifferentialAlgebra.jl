# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex12_linearalgebra.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # Linear algebra with polynomial entries
#
# Factorization, solving, determinants and traces use Julia's LinearAlgebra.
# The constant matrix must be nonsingular for a local polynomial inverse.
using DifferentialAlgebra
using LinearAlgebra

x, y = variables((:x, :y); order = 3)
A = [1 + 0.2x x - y; 0.5y 2 - 0.1x + y]
B = [x one(x); y x + y]
b = [1 + x, -0.5 + y]
F = lu(A)
z = F \ b
Ainverse = F \ Matrix{typeof(x)}(I, 2, 2)
solve_error = maximum(coefficient_norm, A * z - b)
inverse_error = maximum(coefficient_norm, A * Ainverse - Matrix{typeof(x)}(I, 2, 2))
@assert max(solve_error, inverse_error) < 1.0e-13
println("Solution of A z = b: ", z)
println("trace(A B) = ", tr(A * B))
println("det(A) = ", det(A))
println("Frobenius norm of A B: ", sqrt(sum(abs2, A * B)))
println((solve_error = solve_error, inverse_error = inverse_error))

# Compose each matrix entry with a shifted coordinate system.
shifted = [evaluate(p, [x + 0.1, y - 0.2]) for p in A]
println("A(x+0.1, y-0.2) = ", shifted)
