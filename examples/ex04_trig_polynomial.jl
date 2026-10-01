# Adapted for DifferentialAlgebra.jl from DACEjl; algorithms, API, plots and checks modified. #src
# Source: https://github.com/arma1978/DACEjl/blob/c5d062d277a28b02c89e587e4eed098fd3331036/examples/ex04_trig_polynomial.jl #src
# SPDX-License-Identifier: Apache-2.0; see LICENSE-DACEjl and NOTICE-DACEjl. #src
# # Truncation and nilpotent polynomials
#
# cos(x) - 1 starts at degree two. Its eleventh power starts at degree 22,
# so it vanishes in an order-20 algebra. This is an algebraic truncation,
# not a claim that the underlying function is identically zero.
using DifferentialAlgebra

x, = variables((:x,); order = 20)
p = cos(x) - 1
q = p^11
@assert iszero(q)
println("cos(x) - 1 = ", p)
println("(cos(x) - 1)¹¹, truncated at order 20 = ", q)
