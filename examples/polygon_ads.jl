# # Oriented polygon ADS and interval fitting bands
# Error method and domain geometry are independent options. IntervalBound
# validates the uniform fitting error; the other methods estimate it from
# coefficients and optional point checks. Automatic directions are heuristic.

using DifferentialAlgebra, IntervalArithmetic, CairoMakie

f(v) = (v[1] + v[2])^2
box = fill(interval(Float64, -1, 1), 2)
tolerance = 1 // 16
axis = adaptive_map(f, box; estimator = IntervalBound(), order = 1, atol = tolerance)
oriented = adaptive_map(f, box; estimator = IntervalBound(), splitter = :oriented, order = 1, atol = tolerance)
manual = adaptive_map(f, box; estimator = IntervalBound(), splitter = :oriented, directions = [1 1; -1 1], order = 1, atol = tolerance)
@assert axis.converged && oriented.converged && manual.converged
@assert sum(p -> domain_area(domain(p)), oriented.patches) == 4
println("Automatically chosen projection rows: ", split_directions(oriented))
println("Axis / automatic / manual leaf counts: ", length.((axis.patches, oriented.patches, manual.patches)))

# The exact analytical identity is f=z1² for z1=x+y. In a leaf with
# z1=c1+r1*xi1, the omitted term is r1²*xi1², enclosed by [0,r1²].
# This proves a uniform error bound. Sample curves below illustrate the fit;
# they are not used to decide inclusion or acceptance.
for (label, map) in (("Axis", axis), ("Oriented", oriented))
    println((label, largest_uniform_error = maximum(p -> sup(abs(only(p.error_bounds))), map.patches), enclosure_width = diam(enclose(map))))
end

# # Different polygons and ordinary estimators
# The same geometry also works with ordinary ADS. Ordinary queries return
# approximate numbers; interval queries retain coefficients, domain and remainder.
for estimator in (GuardedTail(), ExtrapolatedTail(), LastTerms())
    fit = adaptive_map(f, [-1, -1], [1, 1]; estimator, splitter = :oriented, order = 1, atol = tolerance)
    println((estimator, polygons = length(fit.patches), value = fit([1 // 2, 1 // 4])))
end
triangle = ConvexPolygon([(0, 0), (1, 0), (0, 1)])
triangular = adaptive_map(v -> exp(v[1] + v[2]), triangle; estimator = IntervalBound(), order = 3, atol = 1 // 1000, directions = [1 1; -1 1])
println("Convex triangle leaves: ", length(triangular.patches))

# # Plot physical domains and whole-cell interval bands
# Both partitions cover exactly the same square. Polygon clipping and inverse
# coordinates use exact rationals. Interval arithmetic uses outward rounding
# on each polygon's bounding parallelogram. Covers may extend outside a polygon.
fig = Figure(size = (1150, 800), fontsize = 14)
ax = Axis(fig[1, 1]; xlabel = "physical x", ylabel = "physical y", title = "Axis-aligned: $(length(axis.patches)) boxes", aspect = DataAspect())
for patch in axis.patches
    x, y = domain(patch)
    poly!(ax, Rect2f(inf(x), inf(y), diam(x), diam(y)); color = :lightblue, strokecolor = :steelblue, strokewidth = 1)
end
ax = Axis(fig[1, 2]; xlabel = "physical x", ylabel = "physical y", title = "Automatic diagonal: $(length(oriented.patches)) polygons", aspect = DataAspect())
for patch in oriented.patches
    points = [Point2f(Float64.(v)) for v in polygon_vertices(domain(patch))]
    poly!(ax, points; color = :palegreen, strokecolor = :darkgreen, strokewidth = 1)
end

edges = collect(range(-1, 1; length = 129))
cells = [[interval(edges[i], edges[i + 1]), interval(0)] for i in 1:(length(edges) - 1)]
axis_values = [enclose(axis, b) for b in cells]
oriented_values = [enclose(oriented, b) for b in cells]
plot_x = [x for i in 1:(length(edges) - 1) for x in (edges[i], edges[i + 1])]
lower(values) = [inf(v) for v in values for _ in 1:2]
upper(values) = [sup(v) for v in values for _ in 1:2]
ax = Axis(fig[2, 1]; xlabel = "physical x (y=0)", ylabel = "function value", title = "Enclosures on complete plotted cells")
band!(ax, plot_x, lower(axis_values), upper(axis_values); color = (:dodgerblue, 0.45), label = "Axis-aligned interval ADS")
band!(ax, plot_x, lower(oriented_values), upper(oriented_values); color = (:green, 0.35), label = "Oriented interval ADS")
lines!(ax, edges, edges .^ 2; color = :black, label = "Original function (samples)")
axislegend(ax; position = :ct, labelsize = 11)
ax = Axis(fig[2, 2]; xlabel = "leaf index (different local fits)", ylabel = "uniform fit-error interval", title = "Certified splitting criterion")
for (map, color, offset, label) in ((axis, :dodgerblue, -0.15, "Axis"), (oriented, :green, 0.15, "Oriented"))
    errors = [only(p.error_bounds) for p in map.patches]
    rangebars!(ax, collect(eachindex(errors)) .+ offset, inf.(errors), sup.(errors); color, label, linewidth = 2, whiskerwidth = 5)
end
hlines!(ax, [-Float64(tolerance), Float64(tolerance)]; color = :black, linestyle = :dash, label = "Tolerance")
axislegend(ax; position = :rb, labelsize = 11)
fig

# General convex input domains are supported too. The closed polygon domain
# is retained even though each local expansion uses a parallelogram cover.
polygon_fig = Figure(size = (1200, 450), fontsize = 14)
ax = Axis(polygon_fig[1, 1]; xlabel = "physical x", ylabel = "physical y", title = "exp(x+y) on a convex triangle", aspect = DataAspect())
for patch in triangular.patches
    poly!(ax, [Point2f(Float64.(v)) for v in polygon_vertices(domain(patch))]; color = :palegreen, strokecolor = :darkgreen, strokewidth = 1)
end
grid = range(0, 1; length = 65)
point_widths = [x + y <= 1 ? diam(triangular([x, y])) : NaN for x in grid, y in grid]
ax = Axis(polygon_fig[1, 2]; xlabel = "physical x", ylabel = "physical y", title = "Point enclosure width inside triangle", aspect = DataAspect())
heat = heatmap!(ax, grid, grid, point_widths; colormap = :viridis)
Colorbar(polygon_fig[1, 3], heat; label = "Original-function interval width")
ax = Axis(polygon_fig[1, 4]; xlabel = "polygon leaf", ylabel = "uniform fit-error interval", title = "Fresh leaf certificates")
triangle_errors = [only(p.error_bounds) for p in triangular.patches]
rangebars!(ax, eachindex(triangle_errors), inf.(triangle_errors), sup.(triangle_errors); color = :darkgreen, linewidth = 3, whiskerwidth = 6)
hlines!(ax, [-1 / 1000, 1 / 1000]; color = :black, linestyle = :dash)
polygon_fig

# ## Choosing directions: axis, oblique and diagonal frames
# Projection rows need not be perpendicular. The exact inverse defines each
# cover; it is not replaced by a transpose. These three frames solve the same
# quadratic problem with the same certified tolerance. Good alignment can
# reduce leaf counts and dependence overestimation, but does not guarantee speed.
frame_names = ("Axis projection", "Oblique projection", "Diagonal projection")
frame_fits = (
    adaptive_map(f, box; estimator = IntervalBound(), splitter = :oriented, directions = :axes, order = 1, atol = tolerance),
    adaptive_map(f, box; estimator = IntervalBound(), splitter = :oriented, directions = [1 1 // 2; 0 1], order = 1, atol = tolerance),
    manual,
)
@assert all(m -> sum(p -> domain_area(domain(p)), m.patches) == 4, frame_fits)
direction_fig = Figure(size = (1130, 720), fontsize = 14)
for (column, (name, fit)) in enumerate(zip(frame_names, frame_fits))
    ax = Axis(direction_fig[1, column]; xlabel = "physical x", ylabel = "physical y", title = "$name\n$(length(fit.patches)) polygons", aspect = DataAspect())
    for p in fit.patches
        poly!(ax, [Point2f(Float64.(v)) for v in polygon_vertices(domain(p))]; color = (:palegreen, 0.4), strokecolor = :darkgreen, strokewidth = 0.8)
    end
end
labels = ["Axis", "Oblique", "Diagonal"]
ax = Axis(direction_fig[2, 1]; xticks = (1:3, labels), ylabel = "leaf count", title = "Partition size")
barplot!(ax, 1:3, [length(m.patches) for m in frame_fits]; color = :darkgreen)
ax = Axis(direction_fig[2, 2]; xticks = (1:3, labels), ylabel = "full-domain enclosure width", title = "Dependence and cover overestimation")
barplot!(ax, 1:3, [diam(enclose(m)) for m in frame_fits]; color = :darkgreen)
hlines!(ax, [4]; color = :black, linestyle = :dash, label = "Analytical range width: [0,4]")
axislegend(ax; position = :rt, labelsize = 10)
ax = Axis(direction_fig[2, 3]; xticks = (1:3, labels), ylabel = "uniform error upper bound", title = "Each frame retains its certificate")
scatter!(ax, 1:3, [maximum(p -> sup(abs(only(p.error_bounds))), m.patches) for m in frame_fits]; color = :darkgreen, markersize = 14)
hlines!(ax, [Float64(tolerance)]; color = :black, linestyle = :dash, label = "Requested atol")
axislegend(ax; position = :lb, labelsize = 10)
direction_fig

# The standalone script also saves all three displayed figures as shareable PNGs.
if abspath(PROGRAM_FILE) == @__FILE__
    directory = joinpath(@__DIR__, "..", "results")
    mkpath(directory)
    for (name, figure) in (("polygon_ads.png", fig), ("polygon_triangle.png", polygon_fig), ("polygon_directions.png", direction_fig))
        path = joinpath(directory, name)
        save(path, figure; px_per_unit = 2)
        println("Saved interval/polygon plot: ", abspath(path))
    end
end

# A point inside the triangle's bounding square can still be outside the
# physical polygon. The query is rejected, even if its algebraic formula exists.
try
    triangular([3 // 4, 3 // 4])
catch err
    err isa DomainError || rethrow()
    println("Rejected outside-polygon query: ", err.msg)
end

# Having a valid center does not validate log on its full model domain.
try
    adaptive_map(v -> log(1 + v[1]), box; estimator = IntervalBound(), splitter = :oriented, directions = :axes)
catch err
    err isa DomainError || rethrow()
    println("Rejected invalid whole-domain logarithm: ", err.msg)
end

# Fewer leaves do not guarantee faster construction: exact clipping has a cost.
# Run benchmark/polygon_ads.jl for warmed timing, allocations and interval widths
# separately. Polygon geometry currently supports 2D static maps; six-variable
# uncertainty maps continue to use box ADS. Time integration is not validated.
