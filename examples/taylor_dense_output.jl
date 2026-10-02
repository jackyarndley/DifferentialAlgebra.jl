# # Taylor dense output and events
#
# A driven oscillator has an explicit time dependence and an exact solution.
# Use it to check a Taylor solver's values between steps, derivatives of its
# dense output, and the event times found by SciML's continuous callbacks.
using DifferentialAlgebra
using OrdinaryDiffEqVerner
using SciMLBase: successful_retcode, ContinuousCallback, u_modified!
using CairoMakie

set_theme!(palette = (color = Makie.to_colormap(:tab10),))

# ## A nonautonomous equation
#
# For ``q''+\omega^2q=A\cos(\Omega t)`` with zero initial position and
# velocity, ``q=A(\cos(\Omega t)-\cos(\omega t))/(\omega^2-\Omega^2)``.
# The forcing is nonresonant: ``\omega\ne\Omega``. This solution supplies
# an independent reference without integrating another ODE.
function driven_oscillator(u, p, t)
    ω, Ω, amplitude = p
    return [u[2], -ω^2 * u[1] + amplitude * cos(Ω * t)]
end
parameters = (1.0, 0.7, 0.2)
function exact_state(t, p)
    ω, Ω, amplitude = p
    scale = amplitude / (ω^2 - Ω^2)
    return scale * [cos(Ω * t) - cos(ω * t), -Ω * sin(Ω * t) + ω * sin(ω * t)]
end
problem = ODEProblem(driven_oscillator, [0.0, 0.0], (0.0, 12π), parameters)

# ## Locate zero crossings with dense output
#
# A callback finds roots of q(t) inside accepted steps. Both crossing directions
# are recorded, and the state is unchanged. ``\cos a-\cos b`` factors into two
# sines, so the positive crossing times are the union of
# ``2\pi k/(\omega+\Omega)`` and ``2\pi k/(\omega-\Omega)``.
# No two positive roots coincide within this integration interval.
crossings = Float64[]
function record_crossing!(integrator)
    integrator.t > 1.0e-8 && push!(crossings, integrator.t)
    u_modified!(integrator, false)
    return nothing
end
callback = ContinuousCallback(
    (u, t, integrator) -> u[1], record_crossing!, record_crossing!;
    save_positions = (false, false), abstol = 1.0e-13, reltol = 0
)
solution = solve(problem, TaylorMethod(18); abstol = 1.0e-12, reltol = 1.0e-12, callback)
@assert successful_retcode(solution)
ω, Ω, amplitude = parameters
expected_crossings = sort(
    vcat(
        [
            [2π * k / frequency for k in 1:floor(Int, last(problem.tspan) * frequency / (2π))]
                for frequency in (ω + Ω, ω - Ω)
        ]...
    )
)
@assert length(crossings) == length(expected_crossings)
event_errors = abs.(crossings - expected_crossings)
@assert maximum(event_errors) < 1.0e-9
println("Accepted Taylor steps: ", solution.stats.naccept)
println("Detected zero crossings: ", length(crossings))
println("Maximum event-time error: ", maximum(event_errors))

# ## Values and derivatives between steps
#
# `solution(t, Val{1})` differentiates the stored time polynomial. It does not
# estimate derivatives with finite differences. Check both state components
# at 801 uniformly spaced times, including times between accepted steps.
times = range(first(problem.tspan), last(problem.tspan); length = 801)
reference = stack(exact_state(t, parameters) for t in times)
values = Array(solution(times))
state_errors = [maximum(abs, solution(t) - exact_state(t, parameters)) for t in times]
derivative_errors = [
    maximum(abs, solution(t, Val{1}) - driven_oscillator(exact_state(t, parameters), parameters, t))
        for t in times
]
@assert maximum(state_errors) < 1.0e-10
@assert maximum(derivative_errors) < 1.0e-9
println("Maximum dense state error: ", maximum(state_errors))
println("Maximum dense derivative error: ", maximum(derivative_errors))

# ## Time-order convergence at a fixed step size
#
# Fixing the step size separates the effect of polynomial order from adaptive
# step selection. At high orders, Float64 roundoff eventually limits accuracy.
orders = [6, 10, 14, 18]
order_errors = map(orders) do order
    fixed = solve(problem, TaylorMethod(order); adaptive = false, dt = 0.75)
    @assert successful_retcode(fixed)
    maximum(maximum(abs, fixed(t) - exact_state(t, parameters)) for t in times)
end
@assert last(order_errors) < 1.0e-10
@assert first(order_errors) > 100last(order_errors)
println("Time orders: ", orders)
println("Maximum fixed-step errors: ", order_errors)

#-
fig = Figure(size = (1080, 760), fontsize = 15)
state_axis = Axis(fig[1, 1]; xlabel = "Time", ylabel = "Position q", title = "A driven oscillator")
lines!(state_axis, times, values[1, :]; linewidth = 2, label = "Taylor solution")
lines!(state_axis, times, reference[1, :]; color = :black, linestyle = :dash, label = "Exact solution")
scatter!(state_axis, crossings, zeros(length(crossings)); color = Makie.to_colormap(:tab10)[2], markersize = 8, label = "Located crossing")
axislegend(state_axis; position = :lb, labelsize = 11)
error_axis = Axis(fig[1, 2]; xlabel = "Time", ylabel = "Maximum component error", yscale = log10, title = "Dense values and derivatives")
error_axis.yticks = (10.0 .^ (-16:-14), ["10⁻¹⁶", "10⁻¹⁵", "10⁻¹⁴"])
ylims!(error_axis, 1.0e-16, 1.0e-14)
lines!(error_axis, times, max.(state_errors, eps(Float64)); label = "State", linewidth = 2)
lines!(error_axis, times, max.(derivative_errors, eps(Float64)); label = "Time derivative", linewidth = 2)
axislegend(error_axis; position = :lb)
event_axis = Axis(fig[2, 1]; xlabel = "Crossing time", ylabel = "Absolute time error", yscale = log10, title = "Continuous callback accuracy")
event_axis.yticks = (10.0 .^ (-16:-13), ["10⁻¹⁶", "10⁻¹⁵", "10⁻¹⁴", "10⁻¹³"])
ylims!(event_axis, 1.0e-16, 1.0e-13)
scatter!(event_axis, expected_crossings, max.(event_errors, eps(Float64)); markersize = 10)
order_axis = Axis(fig[2, 2]; xlabel = "Time polynomial order", ylabel = "Maximum state error", yscale = log10, xticks = orders, title = "Fixed steps: Δt = 0.75")
scatterlines!(order_axis, orders, order_errors; linewidth = 2, markersize = 10)
fig
