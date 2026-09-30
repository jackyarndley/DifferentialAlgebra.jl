# # Gradient of the sombrero function
#
# Tutorial 1, exercise 8 from DACE.jl: expand sin(r)/r about (2, 3) and
# extract its gradient using first-order differential algebra.

using DifferentialAlgebra

function sombrero(x)
    r = sqrt(x[1]^2 + x[2]^2)
    return sin(r) / r
end

dx, dy = variables((:δx, :δy); order = 1)
x = [2.0 + dx, 3.0 + dy]
z = sombrero(x)
grad_z = gradient(z)

println("Coordinates about (2, 3): ", x)
println("Sombrero function: ", z)
println("Gradient: ", grad_z)

# Check against the analytic gradient at the expansion center.
r = sqrt(13.0)
expected = (r * cos(r) - sin(r)) / r^3 .* [2.0, 3.0]
@assert constant_term(grad_z) ≈ expected
