# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # Adaptive Runge–Kutta propagation of a Taylor map
#
# An eccentric planar orbit starts at [1, 0, 0, 1.1], with uncertainty
# ±0.01 in its initial position. Integrate its order-three Taylor map.
# This self-contained Dormand–Prince 5(4) implementation exposes the stages,
# embedded error estimate, rejection rule, and step-size controller.
using DifferentialAlgebra
using CairoMakie

function kepler_rhs(u)
    factor = -(u[1]^2 + u[2]^2)^(-3 / 2)
    return [u[3], u[4], factor * u[1], factor * u[2]]
end

magnitude(x::Real) = abs(x)
magnitude(x::TaylorPolynomial) = coefficient_norm(x)

function dopri5(rhs, initial, times; tolerance = 1.0e-11)
    ## Each row contains the explicit stage weights. The final row gives
    ## the fifth-order solution; the embedded weights give order four.
    rows = (
        (), (1 / 5,), (3 / 40, 9 / 40), (44 / 45, -56 / 15, 32 / 9),
        (19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729),
        (9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656),
        (35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84),
    )
    low = (5179 / 57600, 0.0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40)
    state = copy(initial)
    snapshots = [copy(state)]
    time, step = first(times), 0.01
    accepted = rejected = 0
    for target in times[2:end]
        attempts = 0
        while time < target
            attempts += 1
            attempts <= 100000 || error("Too many integration steps")
            h = min(step, target - time)
            time + h > time || error("Step-size underflow")
            stages = [rhs(state)]
            trial = state
            for i in 2:7
                trial = copy(state)
                for j in 1:(i - 1)
                    iszero(rows[i][j]) || (trial += (h * rows[i][j]) * stages[j])
                end
                push!(stages, rhs(trial))
            end
            lower = copy(state)
            for j in 1:7
                lower += (h * low[j]) * stages[j]
            end
            ## Control all retained coefficients, including uncertainty terms.
            ratio = maximum(magnitude, trial - lower) /
                (tolerance * (1 + maximum(magnitude, trial)))
            isfinite(ratio) || error("Nonfinite embedded error")
            if ratio <= 1
                state = trial
                time = h == target - time ? target : time + h
                accepted += 1
            else
                rejected += 1
            end
            step = h * (iszero(ratio) ? 5.0 : clamp(0.9ratio^(-1 / 5), 0.2, 5.0))
        end
        push!(snapshots, copy(state))
    end
    return (; snapshots, accepted, rejected)
end

δx, δy = variables((:δx, :δy); order = 3)
initial = [1 + 0.01δx, 0.01δy, zero(δx), 1.1 + zero(δx)]
period = 2π / (2 - 1.1^2)^(3 / 2)
times = collect(range(0, period; length = 9))
result = dopri5(kepler_rhs, initial, times)
nominal = dopri5(kepler_rhs, [1.0, 0.0, 0.0, 1.1], collect(range(0, period; length = 241)))
closure = maximum(abs, constant_term.(result.snapshots[end]) - [1, 0, 0, 1.1])
@assert closure < 1.0e-8
tighter = dopri5(kepler_rhs, initial, times; tolerance = 1.0e-12)
coefficient_change = maximum(coefficient_norm, result.snapshots[end] - tighter.snapshots[end])
@assert coefficient_change < 1.0e-8
println(
    (
        period = period, accepted = result.accepted, rejected = result.rejected,
        closure_error = closure, coefficient_change = coefficient_change,
    )
)
println("Final flow map: ", result.snapshots[end])

# Map the boundary of the normalized uncertainty square at each snapshot.
edge = range(-1, 1; length = 21)
boundary = vcat(
    [[x, -1.0] for x in edge], [[1.0, y] for y in edge],
    [[x, 1.0] for x in reverse(edge)], [[-1.0, y] for y in reverse(edge)]
)
fig = Figure(size = (780, 620), fontsize = 15)
ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "y", title = "Adaptive Taylor flow: one eccentric orbit", aspect = DataAspect())
lines!(ax, first.(nominal.snapshots), getindex.(nominal.snapshots, 2); color = :black, linewidth = 2, label = "Nominal orbit")
for (i, map) in enumerate(result.snapshots)
    states = evaluate.(Ref(compile(map)), boundary)
    lines!(ax, first.(states), getindex.(states, 2); color = (i - 1) / (length(times) - 1), colorrange = (0, 1), colormap = :viridis, linewidth = 2)
end
scatter!(ax, [0.0], [0.0]; color = :orange, markersize = 14)
Colorbar(fig[1, 2]; limits = (0, 1), colormap = :viridis, label = "Time / nominal period")
fig
