# # Long-term CR3BP propagation
#
# Propagate a planar libration near the Earth–Moon triangular point L₄ for
# **200 revolutions of the primaries**. Compare two Taylor time orders with
# Vern9 using the same requested tolerances. Examine both conservation of the
# Jacobi constant and trajectory error: a small invariant drift alone does not
# establish an accurate trajectory.
using DifferentialAlgebra
using OrdinaryDiffEqVerner
using SciMLBase: successful_retcode
using LinearAlgebra
using CairoMakie

set_theme!(palette = (color = Makie.to_colormap(:tab10),))

# ## Rotating-frame equations and normalization
#
# Distance between the primaries, their total mass, and angular speed are one.
# Thus one primary revolution takes 2π time units. The Earth is at (-μ,0),
# the Moon at (1-μ,0), and ``u=(x,y,v_x,v_y)`` contains rotating-frame velocities.
# The normalization and mass ratio follow [JPL's CR3BP convention](https://ssd-api.jpl.nasa.gov/doc/periodic_orbits.html).
# This is the ideal circular model, rather than an Earth–Moon ephemeris.
#
# With ``r_1^2=(x+\mu)^2+y^2`` and ``r_2^2=(x-1+\mu)^2+y^2``,
# ``\ddot{x}=2\dot{y}+x-(1-\mu)(x+\mu)/r_1^3-\mu(x-1+\mu)/r_2^3``
# and ``\ddot{y}=-2\dot{x}+y-(1-\mu)y/r_1^3-\mu y/r_2^3``.
function cr3bp(u, μ, t)
    x, y, vx, vy = u
    dx1, dx2 = x + μ, x - 1 + μ
    r13 = (dx1^2 + y^2)^(3 // 2)
    r23 = (dx2^2 + y^2)^(3 // 2)
    return [
        vx, vy,
        2vy + x - (1 - μ) * dx1 / r13 - μ * dx2 / r23,
        -2vx + y - (1 - μ) * y / r13 - μ * y / r23,
    ]
end
μ = 0.01215058560962404
l4 = [0.5 - μ, sqrt(3) / 2, 0.0, 0.0]
initial = l4 + [0.01, 0, 0, -0.01]
@assert maximum(abs, cr3bp(l4, μ, 0.0)) < 1.0e-14
revolutions = 200
final_time = revolutions * 2π
times = range(0, final_time; length = 8001)
problem = ODEProblem(cr3bp, initial, (0.0, final_time), μ)

# ## Which energy is conserved?
#
# The [Jacobi constant](https://farside.ph.utexas.edu/teaching/celestial/Celestial/node83.html)
# is ``C=x^2+y^2+2(1-\mu)/r_1+2\mu/r_2-v_x^2-v_y^2``.
# The rotating-frame energy is ``H=-C/2``; inertial Kepler energy is not
# conserved. We report ΔC, so the corresponding energy drift is -ΔC/2.
function jacobi(u, μ)
    x, y, vx, vy = u
    r1 = sqrt((x + μ)^2 + y^2)
    r2 = sqrt((x - 1 + μ)^2 + y^2)
    return x^2 + y^2 + 2 * ((1 - μ) / r1 + μ / r2) - vx^2 - vy^2
end

# All integrations below use Float64. Evaluating only the diagnostic with
# 128-bit arithmetic avoids hiding small drifts in subtraction of two
# nearly equal, rounded Float64 Jacobi constants. It does not improve the
# trajectories themselves.
function jacobi_drift(states, initial, μ)
    return setprecision(128) do
        mass = BigFloat(μ)
        initial_constant = jacobi(BigFloat.(initial), mass)
        [Float64(jacobi(BigFloat.(u), mass) - initial_constant) for u in states]
    end
end

# ## Independent trajectory checks
#
# Use a tighter Vern9 solution as the reference, and cross-check it with a
# tighter, order-28 Taylor solution. Neither solution is assumed exact.
# `saveat` samples each solver's dense interpolant on the same time grid;
# the adaptive step sequence remains free to differ between algorithms.
reference = solve(problem, Vern9(); abstol = 2.0e-14, reltol = 2.0e-14, saveat = times, dense = false)
crosscheck = solve(problem, TaylorMethod(28); abstol = 2.0e-14, reltol = 2.0e-14, saveat = times, dense = false)
@assert successful_retcode(reference) && successful_retcode(crosscheck)
reference_difference = maximum(maximum(abs, u - v) for (u, v) in zip(reference.u, crosscheck.u))
@assert reference_difference < 1.0e-9
println("Primary revolutions: ", revolutions)
println("Reference cross-check, maximum state difference: ", reference_difference)
println("Initial Jacobi constant: ", jacobi(initial, μ))

algorithms = ["Taylor order 12" => TaylorMethod(12), "Taylor order 20" => TaylorMethod(20), "Vern9" => Vern9()]
solutions = [solve(problem, alg; abstol = 1.0e-10, reltol = 1.0e-10, saveat = times, dense = false) for (_, alg) in algorithms]
drifts = [jacobi_drift(sol.u, initial, μ) for sol in solutions]
state_errors = [[maximum(abs, u - ref) for (u, ref) in zip(sol.u, reference.u)] for sol in solutions]
for (i, (name, _)) in enumerate(algorithms)
    sol = solutions[i]
    @assert successful_retcode(sol)
    @assert last(sol.t) == final_time
    @assert maximum(norm(u[1:2] - l4[1:2]) for u in sol.u) < 0.12
    @assert maximum(abs, drifts[i]) < 1.0e-8
    @assert maximum(state_errors[i]) < 2.0e-7
    println(
        name, ": accepted steps = ", sol.stats.naccept,
        ", maximum |ΔC| = ", maximum(abs, drifts[i]),
        ", final signed ΔC = ", last(drifts[i]),
        ", maximum state difference = ", maximum(state_errors[i])
    )
end

# TaylorMethod is not symplectic and does not enforce C. These plots describe
# this smooth, bounded trajectory and these tolerances, not a general ordering
# of the methods. Accepted step counts also do not measure wall-clock cost.
# The companion [CR3BP uncertainty tutorial](taylor_cr3bp_maps.md) propagates
# a polynomial flow map over the same 200-revolution horizon.

#-
fig = Figure(size = (1120, 780), fontsize = 15)
colors = Makie.to_colormap(:tab10)
orbit_axis = Axis(fig[1, 1]; xlabel = "x - x(L₄)", ylabel = "y - y(L₄)", aspect = DataAspect(), title = "200 primary revolutions near L₄")
lines!(orbit_axis, getindex.(reference.u, 1) .- l4[1], getindex.(reference.u, 2) .- l4[2]; color = (colors[1], 0.45), linewidth = 0.7)
scatter!(orbit_axis, [0.0], [0.0]; color = :black, marker = :star5, markersize = 15, label = "L₄")
scatter!(orbit_axis, [initial[1] - l4[1]], [initial[2] - l4[2]]; color = colors[2], markersize = 10, label = "Initial state")
axislegend(orbit_axis; position = :lb, labelsize = 11)
drift_axis = Axis(fig[1, 2]; xlabel = "Primary revolutions", ylabel = "Absolute Jacobi drift |ΔC|", yscale = log10, title = "Conservation at sampled times")
error_axis = Axis(fig[2, 1]; xlabel = "Primary revolutions", ylabel = "Maximum state difference", yscale = log10, title = "Comparison with tighter Vern9")
for (i, (name, _)) in enumerate(algorithms)
    lines!(drift_axis, times ./ (2π), max.(abs.(drifts[i]), 1.0e-17); color = colors[i], linewidth = 1.5, label = name)
    lines!(error_axis, times ./ (2π), max.(state_errors[i], 1.0e-16); color = colors[i], linewidth = 1.5, label = name)
end
axislegend(drift_axis; position = :lt, labelsize = 11)
axislegend(error_axis; position = :lt, labelsize = 11)
step_axis = Axis(fig[2, 2]; ylabel = "Accepted steps", xticks = (1:3, first.(algorithms)), title = "Adaptive step counts")
barplot!(step_axis, 1:3, [sol.stats.naccept for sol in solutions]; color = colors[1:3])
fig
