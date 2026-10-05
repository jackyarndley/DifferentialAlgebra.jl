# # Optimization with a C2 oriented ADS surrogate
# A smooth objective can be useful to gradient and Newton methods. Here a
# nonlinear convex function has diagonal sensitivity, so we compare box ADS
# with an oriented polygon partition, then optimize a C2 overlapping blend.
# IntervalBound certifies function errors; this is not verified optimization.

using DifferentialAlgebra, IntervalArithmetic, ForwardDiff, LinearAlgebra, CairoMakie

function objective(v)
    u = v[1] + v[2] - 1 // 10
    w = v[1] - v[2] - 3 // 10
    return u^2 + 4w^2 + (exp(u) + exp(-u) - 2) / 10
end
# u=w=0 gives the exact unique minimizer (1/5,-1/10), with value zero.
# The Hessian is positive definite because the quadratic part is, and the
# exponential term is convex. This independent analytical solution is used
# for comparison, not as a certificate for the surrogate's derivatives.
box = fill(interval(-1, 1), 2)
tolerance = 1 // 10000
boxes = adaptive_map(objective, box; estimator = IntervalBound(), order = 3, atol = tolerance)
polygons = adaptive_map(objective, box; estimator = IntervalBound(), splitter = :oriented, directions = [1 1; -1 1], order = 3, atol = tolerance)
smooth = continuous_map(objective, polygons; continuity = :c2, order = 4, atol = 1 // 1000)
println("Box / oriented leaves: ", length.((boxes.patches, polygons.patches)))
println("C2 surrogate uniform function-error interval: ", only(smooth.error_bounds))

# Fresh overlap fits are computed at order four. This separate retained order
# can improve derivative accuracy; selecting C2 alone does not do so. The
# wrapper checks its new bound against atol and never silently inherits the
# tighter tolerance of the order-three source.
function damped_newton(f, initial; max_steps = 30, tolerance = 1.0e-8)
    x = copy(initial)
    path = [copy(x)]
    gradients = Float64[]
    for _ in 1:max_steps
        g = ForwardDiff.gradient(f, x)
        push!(gradients, norm(g))
        norm(g) <= tolerance && break
        H = ForwardDiff.hessian(f, x)
        factor = cholesky(Symmetric(H); check = false)
        direction = issuccess(factor) ? -(factor \ g) : -g
        step = 1.0
        accepted = false
        for _ in 1:40
            candidate = x + step * direction
            ## Keep every line-search query inside the declared validity box.
            if all(-1 .<= candidate .<= 1) && f(candidate) <= f(x) + 1.0e-4 * step * dot(g, direction)
                x = candidate
                push!(path, copy(x))
                accepted = true
                break
            end
            step /= 2
        end
        accepted || error("Line search did not find a decreasing in-domain step")
    end
    return (; point = x, path, gradients)
end
solution = damped_newton(smooth, [0.8, 0.5])
println("Numeric surrogate minimizer: ", solution.point)
println("Distance to analytical minimizer: ", norm(solution.point - [1 / 5, -1 / 10]))
println("Original objective / enclosure there: ", (objective(solution.point), enclose(smooth, solution.point)))
@assert norm(ForwardDiff.gradient(smooth, solution.point)) < 1.0e-8
@assert norm(solution.point - [1 / 5, -1 / 10]) < 1.0e-3

# # Physical partitions and optimization trajectory
grid = collect(range(-1, 1; length = 81))
levels = [objective([x, y]) for x in grid, y in grid]
fig = Figure(size = (1100, 490), fontsize = 14)
ax = Axis(fig[1, 1]; xlabel = "physical x", ylabel = "physical y", title = "Box ADS: $(length(boxes.patches)) leaves", aspect = DataAspect())
for p in boxes.patches
    x, y = domain(p)
    poly!(ax, Rect2f(inf(x), inf(y), diam(x), diam(y)); color = (:lightblue, 0.6), strokecolor = :steelblue, strokewidth = 1)
end
contour!(ax, grid, grid, levels; levels = [0.1, 1, 4, 10, 20], color = :gray, linewidth = 1)
ax = Axis(fig[1, 2]; xlabel = "physical x", ylabel = "physical y", title = "C2 oriented fit: $(length(polygons.patches)) leaves", aspect = DataAspect())
for p in polygons.patches
    poly!(ax, [Point2f(Float64.(v)) for v in polygon_vertices(domain(p))]; color = (:palegreen, 0.6), strokecolor = :darkgreen, strokewidth = 1)
end
contour!(ax, grid, grid, levels; levels = [0.1, 1, 4, 10, 20], color = :gray, linewidth = 1)
lines!(ax, first.(solution.path), last.(solution.path); color = :orange, linewidth = 3, label = "Damped Newton")
scatter!(ax, first.(solution.path), last.(solution.path); color = :orange, markersize = 9)
scatter!(ax, [1 / 5], [-1 / 10]; color = :black, marker = :star5, markersize = 16, label = "Analytical minimizer")
axislegend(ax; position = :lt, labelsize = 11)
fig

# Original-function intervals at each iterate retain the full overlap remainder.
# Their widths are reported separately from the numerical surrogate gradient.
# Function-value intervals do not prove that the iterates minimize the function.
certificates = [enclose(smooth, p) for p in solution.path]
values_fig = Figure(size = (1360, 390), fontsize = 14)
ax = Axis(values_fig[1, 1]; xlabel = "Newton iterate", ylabel = "objective value", title = "Original function with retained interval bounds")
iterations = 0:(length(solution.path) - 1)
rangebars!(ax, iterations, inf.(certificates), sup.(certificates); color = :dodgerblue, linewidth = 3, whiskerwidth = 10, label = "Original-function enclosure")
scatter!(ax, iterations, objective.(solution.path); color = :black, label = "Original value (rounded)")
scatter!(ax, iterations, smooth.(solution.path); color = :orange, marker = :cross, label = "C2 surrogate")
axislegend(ax; position = :rt, labelsize = 11)
ax = Axis(values_fig[1, 2]; xlabel = "Newton iterate", ylabel = "original function enclosure", title = "Final iterates: interval detail")
final_iterations = max(0, length(solution.path) - 2):(length(solution.path) - 1)
final_indices = collect(final_iterations) .+ 1
rangebars!(ax, final_iterations, inf.(certificates[final_indices]), sup.(certificates[final_indices]); color = :dodgerblue, linewidth = 3, whiskerwidth = 10)
scatter!(ax, final_iterations, objective.(solution.path[final_indices]); color = :black)
hlines!(ax, [0]; color = (:gray, 0.5), linestyle = :dash)
ax = Axis(values_fig[1, 3]; xlabel = "Newton iterate", ylabel = "surrogate gradient norm", yscale = log10, title = "Numerical convergence of the smooth objective")
scatterlines!(ax, 0:(length(solution.gradients) - 1), max.(solution.gradients, eps(Float64)); color = :orange)
values_fig

# Exact clipping avoids holes in the original partition; the compact weight
# proof gives C2 at its faces and junctions. There is no claim that oriented ADS
# is always faster. See benchmark/continuous_ads.jl for time, allocations,
# and enclosure widths reported separately, including a six-variable box case.
if abspath(PROGRAM_FILE) == @__FILE__
    directory = joinpath(@__DIR__, "..", "results")
    mkpath(directory)
    for (name, figure) in (("ads_optimization.png", fig), ("ads_optimization_intervals.png", values_fig))
        path = joinpath(directory, name)
        save(path, figure; px_per_unit = 2)
        println("Saved optimization plot: ", abspath(path))
    end
end
