# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # A rational function near a nonzero point
#
# The variable δx represents a displacement, not the physical coordinate.
# For f(x) = x/(x²+1), f(3) = 0.3 and f′(3) = -0.08.
using DifferentialAlgebra

δx, = variables((:δx,); order = 1)
x = 3 + δx
f = 1 / (x + 1 / x)
@assert constant_term(f) ≈ 0.3
@assert coefficient(f, [1]) ≈ -0.08
println("x = ", x)
println("f(3 + δx) = ", f)
println("Linear prediction at x = 3.01: ", f(0.01))
println("Direct value: ", 3.01 / (3.01^2 + 1))
