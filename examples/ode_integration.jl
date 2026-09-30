# # Orbit integration and state transition matrix
#
# Propagate a circular orbit and its second-order Taylor expansion through one
# revolution. All six initial position and velocity components are independent
# variables. Install the example environment as described in examples/README.md.

using DifferentialAlgebra
using OrdinaryDiffEqVerner

# Normalized Kepler equations: acceleration = -μ r / |r|³.
function kepler_ode!(du, u, μ, _)
    radius2 = u[1]^2 + u[2]^2 + u[3]^2
    factor = -μ / (radius2 * sqrt(radius2))
    for i in 1:3
        du[i] = u[i + 3]
        du[i + 3] = factor * u[i]
    end
    return nothing
end

μ = 1.0
initial = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
timespan = (0.0, 2π)
problem = ODEProblem(kepler_ode!, initial, timespan, μ)
nominal = solve(
    problem, Vern9(); abstol = 1.0e-12, reltol = 1.0e-12,
    save_everystep = false
)

perturbed = initial .+ variables(6; order = 2)

# Adaptive error control uses constant parts; it is not an error bound for every
# Taylor coefficient; use convergence checks when controlling higher-order terms.
solution = solve(
    remake(problem; u0 = perturbed), Vern9();
    abstol = 1.0e-12, reltol = 1.0e-12, save_everystep = false
)
final = solution.u[end]
constants = constant_term.(final)
@assert maximum(abs.(constants - nominal.u[end])) < 1.0e-9
@assert maximum(abs.(constants - initial)) < 1.0e-9

# The Jacobian at zero perturbation is the state transition matrix.
stm = constant_term.(jacobian(final))
@assert size(stm) == (6, 6)
@assert all(isfinite, stm)
println("Maximum nominal orbit error: ", maximum(abs.(constants - initial)))
println("State transition matrix:")
display(stm)
