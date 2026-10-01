# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # A Kepler flow map with an explicit time variable
#
# Propagate six initial-state perturbations through an eccentric Earth orbit.
# A seventh variable τ then represents time measured from the final epoch.
# Picard iteration Uₖ₊₁(τ)=Ufinal+∫₀ᵗᵃᵘ f(Uₖ(s)) ds determines one
# further time order per iteration, including mixed state/time coefficients.
using DifferentialAlgebra
using OrdinaryDiffEqVerner
using SciMLBase: successful_retcode
using CairoMakie

function kepler_rhs(u, μ)
    factor = -μ * (u[1]^2 + u[2]^2 + u[3]^2)^(-3 / 2)
    return [u[4], u[5], u[6], factor * u[1], factor * u[2], factor * u[3]]
end
function kepler_ode!(du, u, μ, time)
    du .= kepler_rhs(u, μ)
    return nothing
end

order = 4
δ = variables((:δx, :δy, :δz, :δvx, :δvy, :δvz, :τ); order)
μ, r0, eccentricity = 398600.0, 6678.0, 0.5
nominal = [r0, 0.0, 0.0, 0.0, sqrt(μ / r0) * sqrt(1 + eccentricity), 0.0]
period = 2π * sqrt((r0 / (1 - eccentricity))^3 / μ)
problem = ODEProblem(kepler_ode!, nominal + δ[1:6], (0.0, period), μ)
solution = solve(problem, Vern9(); abstol = 1.0e-11, reltol = 1.0e-12, save_everystep = false)
@assert successful_retcode(solution)
final = solution.u[end]
time_map = copy(final)
for _ in 1:order
    global time_map = final + integrate(kepler_rhs(time_map, μ), 7)
end

# Differentiate and compare through order three; the fourth time derivative
# would require a fifth-order time expansion.
residual = differentiate(time_map, 7) - kepler_rhs(time_map, μ)
residual = DifferentialAlgebra.trim(residual, 0, order - 1)
@assert maximum(coefficient_norm, residual) < 1.0e-9
@assert maximum(coefficient_norm, DifferentialAlgebra.substitute(time_map, 7, 0) - final) < 1.0e-12
println("Orbital period (s): ", period)
println("Maximum Picard coefficient residual through order three: ", maximum(coefficient_norm, residual))
for (i, component) in enumerate(time_map)
    println("Time-expanded state component ", i, ": ", component)
end

# Test a nonzero state perturbation together with positive and negative time
# offsets against fresh integrations, without using the time-expanded map.
perturbation = [0.01, -0.02, 0.01, 1.0e-5, -2.0e-5, 1.0e-5]
offsets = range(-10, 10; length = 41)
errors = Float64[]
for τ in offsets
    reference = solve(
        remake(problem; u0 = nominal + perturbation, tspan = (0.0, period + τ)),
        Vern9(); abstol = 1.0e-11, reltol = 1.0e-12, save_everystep = false
    )
    @assert successful_retcode(reference)
    prediction = evaluate(time_map, vcat(perturbation, τ))
    push!(errors, maximum(abs, prediction[1:3] - reference.u[end][1:3]))
end
@assert maximum(errors) < 1.0e-4
println("Maximum position error (km): ", maximum(errors))

#-
fig = Figure(size = (780, 410), fontsize = 15)
ax = Axis(fig[1, 1]; xlabel = "Time offset τ (s)", ylabel = "Maximum position error (km)", yscale = log10, title = "Fourth-order state/time Kepler map")
lines!(ax, offsets, max.(errors, eps(Float64)); linewidth = 2)
fig
