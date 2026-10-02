# # CR3BP flow maps and uncertainty propagation
#
# Propagate uncertain initial position and velocity through **200 primary
# revolutions** near Earth–Moon L₄. Build a fifth-order flow map with adaptive
# Taylor integration, evaluate thousands of samples, and compute the mean and
# covariance of the polynomial map. Check the predictions against fresh
# numerical trajectories, including the corners of the uncertainty box.
using DifferentialAlgebra
using OrdinaryDiffEqVerner
using SciMLBase: successful_retcode
using LinearAlgebra, Random
using CairoMakie

set_theme!(palette = (color = Makie.to_colormap(:tab10),))

# ## A two-dimensional uncertainty box
#
# The model and nominal state match [Long-term CR3BP propagation](taylor_cr3bp.md).
# The equations are repeated here so this script runs independently. Positions
# and rotating-frame velocities are nondimensional, following [JPL's convention](https://ssd-api.jpl.nasa.gov/doc/periodic_orbits.html).
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
function jacobi(u, μ)
    x, y, vx, vy = u
    return x^2 + y^2 + 2 * ((1 - μ) / sqrt((x + μ)^2 + y^2) + μ / sqrt((x - 1 + μ)^2 + y^2)) - vx^2 - vy^2
end
μ = 0.01215058560962404
nominal = [0.5 - μ + 0.01, sqrt(3) / 2, 0.0, -0.01]
width = 5.0e-4
initial_state(δ) = nominal + [width * δ[1], 0, 0, width * δ[2]]
final_time = 200 * 2π
checkpoints = range(0, final_time; length = 2001)
δx, δvy = variables((:δx, :δvy); order = 5)
initial_map = initial_state([δx, δvy])
problem = ODEProblem(cr3bp, initial_map, (0.0, final_time), μ)

# δx and δvy are independent uniforms on [-1,1]; y and vx are fixed initially.
# Their initial standard deviations are width/√3, not width. Time order 20
# is independent of uncertainty order 5. Saving only checkpoints avoids
# retaining every adaptive step's dense time expansion.
flow = solve(problem, TaylorMethod(20); abstol = 1.0e-13, reltol = 1.0e-13, saveat = checkpoints, dense = false)
@assert successful_retcode(flow)
@assert last(flow.t) == final_time
final_map = last(flow.u)
orders = [1, 3, 5]
maps = [CompiledMap(DifferentialAlgebra.trim.(final_map, 0, order)) for order in orders]
println("Time order: ", 20, "; uncertainty order: ", max_order())
println("Accepted polynomial steps: ", flow.stats.naccept)
println("Final x map: ", final_map[1])
println("Final y map: ", final_map[2])
# The linear coefficients differentiate with respect to normalized δ. Divide
# by width to obtain physical sensitivities to initial x and vy. This is a
# 4×2 submatrix of the state transition matrix, not the full 4×4 STM.
println("Sensitivity to initial (x, vy): ", linear_part(final_map) / width)

# ## Independent sample validation
#
# Every reference below is a new Vern9 solve with a numeric initial state.
# Check all nine points of the {-1,0,1}² grid and 64 seeded random points.
# The errors contain both time-integration error and uncertainty truncation.
rng = Xoshiro(2026)
corners = stack([x, y] for x in (-1.0, 0.0, 1.0) for y in (-1.0, 0.0, 1.0))
validation_points = hcat(corners, 2rand(rng, 2, 64) .- 1)
function numeric_endpoint(δ)
    sol = solve(
        remake(problem; u0 = initial_state(δ)), Vern9();
        abstol = 2.0e-14, reltol = 2.0e-14, save_everystep = false, dense = false
    )
    @assert successful_retcode(sol) && last(sol.t) == final_time
    return last(sol.u)
end
reference = stack(numeric_endpoint(point) for point in eachcol(validation_points))
predictions = [stack(map(point) for point in eachcol(validation_points)) for map in maps]
errors = [vec(maximum(abs.(values - reference); dims = 1)) for values in predictions]
maximum_errors = maximum.(errors)
@assert maximum_errors[3] < maximum_errors[2] < maximum_errors[1]
@assert last(maximum_errors) < 1.0e-7
for (order, error) in zip(orders, maximum_errors)
    println("Map order ", order, ": maximum validation error = ", error)
end

# Conservation can also be checked coefficient by coefficient through degree
# five. This residual is an invariant of the truncated algebra, not an error
# bound for finite perturbations. Evaluating C on numeric map predictions
# additionally exposes terms discarded by polynomial truncation.
invariant_residual = maximum(coefficient_norm(jacobi(map, μ) - jacobi(initial_map, μ), 1) for map in flow.u)
sampled_drift = maximum(abs(jacobi(predictions[3][:, j], μ) - jacobi(initial_state(point), μ)) for (j, point) in enumerate(eachcol(validation_points)))
@assert invariant_residual < 1.0e-10
@assert sampled_drift < 1.0e-8
println("Maximum polynomial Jacobi residual, coefficient 1-norm: ", invariant_residual)
println("Maximum sampled Jacobi drift of the fifth-order map: ", sampled_drift)

# ## Moments without Monte Carlo sampling error
#
# For an independent uniform coordinate, E[δᵏ] is zero for odd k and 1/(k+1)
# for even k. Multiply these moments across coordinates and sum monomials.
# For covariance, form products of **coefficient lists**: multiplying two
# degree-five polynomials in the current algebra would discard degrees 6–10
# and produce an incorrect second moment.
uniform_moment(powers) = any(isodd, powers) ? 0.0 : prod(k -> 1 / (k + 1), powers)
function uniform_statistics(map)
    mean = [sum(m.coefficient * uniform_moment(m.exponents) for m in monomials(p); init = 0.0) for p in map]
    centered = monomials.(map - mean)
    covariance = [
        sum(
            a.coefficient * b.coefficient * uniform_moment(a.exponents + b.exponents)
                for a in centered[i], b in centered[j]; init = 0.0
        ) for i in eachindex(map), j in eachindex(map)
    ]
    return mean, covariance
end
position_moments = [uniform_statistics(map[1:2]) for map in flow.u]
means = stack(first.(position_moments))
deviations = stack(sqrt.(diag(covariance)) for (_, covariance) in position_moments)
initial_mean, initial_covariance = uniform_statistics(first(flow.u))
final_mean, final_covariance = uniform_statistics(final_map)
@assert isapprox(initial_mean, nominal; atol = 1.0e-14)
@assert isapprox(diag(initial_covariance), [width^2 / 3, 0, 0, width^2 / 3]; atol = 1.0e-18)
@assert minimum(eigvals(Symmetric(final_covariance))) >= -1.0e-16
println("Final polynomial mean: ", final_mean)
println("Final polynomial covariance:")
show(stdout, MIME"text/plain"(), final_covariance)
println()

# Validate the moments with a 6×6 Gauss–Legendre rule applied to independently
# integrated trajectories. Six points integrate polynomials through degree 11
# in each coordinate, so the rule exactly integrates the degree-five map's
# mean and covariance in exact arithmetic. For the true nonlinear flow it is
# a quadrature approximation, not a proof of a distribution-wide error bound.
quadrature_order = 6
quadrature = eigen(SymTridiagonal(zeros(quadrature_order), [k / sqrt(4k^2 - 1) for k in 1:(quadrature_order - 1)]))
nodes, weights = quadrature.values, quadrature.vectors[1, :] .^ 2
quadrature_points = [[x, y] for x in nodes for y in nodes]
joint_weights = [wx * wy for wx in weights for wy in weights]
# First check the moment implementation against quadrature of the polynomial.
map_values = stack(maps[3](point) for point in quadrature_points)
map_mean = map_values * joint_weights
map_centered = map_values .- map_mean
map_covariance = (map_centered .* transpose(joint_weights)) * transpose(map_centered)
@assert maximum(abs, map_mean - final_mean) < 1.0e-13
@assert norm(map_covariance - final_covariance) / norm(final_covariance) < 1.0e-11
numeric_values = stack(numeric_endpoint(point) for point in quadrature_points)
quadrature_mean = numeric_values * joint_weights
centered_values = numeric_values .- quadrature_mean
quadrature_covariance = (centered_values .* transpose(joint_weights)) * transpose(centered_values)
mean_error = maximum(abs, final_mean - quadrature_mean)
covariance_error = norm(final_covariance - quadrature_covariance) / norm(quadrature_covariance)
@assert mean_error < 1.0e-8
@assert covariance_error < 1.0e-4
println("Mean difference from trajectory quadrature: ", mean_error)
println("Relative covariance difference from trajectory quadrature: ", covariance_error)

# ## A sampled cloud and its evolution
#
# Evaluating the compiled map requires no additional ODE solves. All 5,000
# samples stay inside the stated uniform box. The mean and standard-deviation
# histories below come from analytic polynomial moments, not from this cloud.
# They are sampled ten times per primary revolution. Bands show each
# revolution's sampled minimum and maximum; they are not confidence intervals
# or bounds between sampling times.
# A single local map eventually loses accuracy as the domain or horizon grows;
# increase its order, shrink its domain, or use automatic domain splitting,
# and repeat independent validation. The sampled checks do not bound the box.
samples = 2rand(rng, 2, 5000) .- 1
linear_cloud = stack(maps[1](point) for point in eachcol(samples))
nonlinear_cloud = stack(maps[3](point) for point in eachcol(samples))
nominal_history = stack(constant_term.(map[1:2]) for map in flow.u)
mean_shift = means - nominal_history
revolution_times = (1:200) .- 0.5
intervals = [((10k - 9):(10k + 1)) for k in 1:200]

#-
fig = Figure(size = (1140, 800), fontsize = 15)
colors = Makie.to_colormap(:tab10)
cloud_axis = Axis(fig[1, 1]; xlabel = "Final x - nominal x", ylabel = "Final y - nominal y", title = "Propagated uniform uncertainty", aspect = DataAspect())
scatter!(cloud_axis, linear_cloud[1, :] .- nominal_history[1, end], linear_cloud[2, :] .- nominal_history[2, end]; color = (colors[2], 0.25), markersize = 3, label = "Linear map")
scatter!(cloud_axis, nonlinear_cloud[1, :] .- nominal_history[1, end], nonlinear_cloud[2, :] .- nominal_history[2, end]; color = (colors[1], 0.25), markersize = 3, label = "Fifth-order map")
scatter!(cloud_axis, reference[1, :] .- nominal_history[1, end], reference[2, :] .- nominal_history[2, end]; color = :black, marker = :cross, markersize = 7, label = "Independent trajectories")
axislegend(
    cloud_axis,
    [
        MarkerElement(color = colors[2], marker = :circle, markersize = 9),
        MarkerElement(color = colors[1], marker = :circle, markersize = 9),
        MarkerElement(color = :black, marker = :cross, markersize = 9),
    ],
    ["Linear map", "Fifth-order map", "Independent trajectories"];
    position = :lb, labelsize = 11
)
spread_axis = Axis(fig[1, 2]; xlabel = "Primary revolutions", ylabel = "Position standard deviation", title = "Spread: per-revolution range")
bias_axis = Axis(fig[2, 1]; xlabel = "Primary revolutions", ylabel = "Mean minus nominal position", title = "Mean shift: per-revolution range")
for (i, name) in enumerate(("x", "y"))
    for (axis, history) in ((spread_axis, deviations), (bias_axis, mean_shift))
        lower = [minimum(history[i, indices]) for indices in intervals]
        upper = [maximum(history[i, indices]) for indices in intervals]
        band!(axis, revolution_times, lower, upper; color = (colors[i], 0.15))
        lines!(axis, revolution_times, lower; color = colors[i], linewidth = 1)
        lines!(axis, revolution_times, upper; color = colors[i], label = name, linewidth = 1)
    end
end
axislegend(spread_axis; position = :lt)
axislegend(bias_axis; position = :lt)
error_axis = Axis(fig[2, 2]; xlabel = "Maximum state error per validation point", ylabel = "Fraction of validation points", xscale = log10, title = "Independent map validation")
for (i, order) in enumerate(orders)
    sorted = sort(errors[i])
    stairs!(error_axis, max.(sorted, eps(Float64)), (1:length(sorted)) ./ length(sorted); color = colors[i], label = "Degree $order", linewidth = 2)
end
axislegend(error_axis; position = :lt)
fig
