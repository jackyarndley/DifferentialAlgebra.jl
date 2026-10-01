# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex05_diff_integral.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # Differentiation and integration
#
# Differentiation loses one known order. Integration chooses a zero constant
# and discards terms beyond the active order. Compare like orders and constants.
using DifferentialAlgebra

x, = variables((:x,); order = 20)
p = sin(x)
derivative = differentiate(p, 1)
primitive = integrate(p, 1)
expected_derivative = with_order(() -> cos(x), 19)
@assert coefficient_norm(derivative - expected_derivative) < 1.0e-15
@assert coefficient_norm(primitive - (1 - cos(x))) < 1.0e-15
println("d sin(x) / dx = ", derivative)
println("∫ sin(x) dx, with value zero at x=0: ", primitive)
println("Derivative residual through order 19: ", coefficient_norm(derivative - expected_derivative))
println("Primitive residual: ", coefficient_norm(primitive - (1 - cos(x))))
