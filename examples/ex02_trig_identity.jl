# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex02_trig_identity.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # Trigonometric identity
#
# Polynomial arithmetic takes place modulo terms above the working order.
# Thus sin²(x) + cos²(x) = 1 coefficient by coefficient, up to rounding.
using DifferentialAlgebra

x, = variables((:x,); order = 20)
identity = sin(x)^2 + cos(x)^2
residual = coefficient_norm(identity - 1)
@assert residual < 1.0e-14
println("sin²(x) + cos²(x) = ", identity)
println("Maximum coefficient residual: ", residual)
