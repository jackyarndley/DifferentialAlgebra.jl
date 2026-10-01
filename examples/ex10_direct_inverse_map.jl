# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex10_direct_inverse_map.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # Direct and inverse maps
#
# A linear map is its own simplest local inversion problem.
# Composition of f(x)=3x with its inverse must give the identity polynomial.
using DifferentialAlgebra

x, = variables((:x,); order = 2)
forward = [3x]
backward = invert(forward)
composition = evaluate(forward, backward)
@assert coefficient_norm(only(composition) - x) < 1.0e-15
println("Forward map: ", forward)
println("Inverse map: ", backward)
println("Composition: ", composition)

# Nonlinear maps work the same way. A nonsingular linear part is required.
x, y = variables((:x, :y); order = 5)
map = [2x + y + x * y, x + 3y + x^2]
inverse = invert(map)
error = maximum(coefficient_norm, evaluate(map, inverse) - [x, y])
@assert error < 1.0e-13
println("Nonlinear inverse: ", inverse)
println("Maximum composition coefficient residual: ", error)
