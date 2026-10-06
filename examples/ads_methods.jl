# SPDX-License-Identifier: Apache-2.0; modified examples, see ../THIRD_PARTY_NOTICES.md. #src
# # Choosing an automatic domain splitting method
#
# A Gaussian on an anisotropic box illustrates how the error estimator and
# split direction affect the partition. All methods use the same interface,
# polynomial order, tolerance and independent validation grid.
# A second comparison adds IntervalBound and oriented polygons on a diagonal
# nonlinear function. Estimator, geometry and optional continuity are distinct.
using DifferentialAlgebra
using CairoMakie
using IntervalArithmetic

set_theme!(palette = (color = Makie.to_colormap(:tab10),))

gaussian(x) = exp(-(4x[1]^2 + x[2]^2) / 2)
lower, upper = [-2.0, -2.0], [2.0, 2.0]
tolerance = 1.0e-4
methods = [
    ("Guard degrees", GuardedTail(), :tail),
    ("Coefficient decay", ExtrapolatedTail(), :tail),
    ("Last two degrees", LastTerms(), :tail),
    ("Longest relative side", GuardedTail(), :width),
    ("Interval bounds", IntervalBound(), :tail),
]
maps = [
    adaptive_map(gaussian, lower, upper; order = 5, atol = tolerance, estimator, splitter)
        for (_, estimator, splitter) in methods
]

# GuardedTail computes two extra degrees. ExtrapolatedTail fits an exponential
# decay to retained coefficient maxima, and LastTerms inspects the last two
# degrees without fitting. The :tail splitter favors directions that reduce
# large coefficients; :width provides a simple geometric alternative.
# The first four cases use heuristic indicators. IntervalBound instead bounds
# retained coefficient uncertainty and the absolute remainder of fresh child
# Taylor models. Only that fifth acceptance test is a uniform error certificate.
# Compare scalar polynomial values with interval midpoints on an independent
# grid; midpoint discrepancy is a numerical diagnostic, not the interval proof.
grid = range(-2, 2; length = 101)
reference = [gaussian([x, y]) for x in grid, y in grid]
grid_values = [[map([x, y]) for x in grid, y in grid] for map in maps]
@assert all(isguaranteed, grid_values[end])
plot_value(x) = x isa Interval ? mid(x) : x
errors = [abs.(plot_value.(values) - reference) for values in grid_values]
for ((name, _, _), map, error) in zip(methods, maps, errors)
    @assert map.converged
    @assert maximum(error) < tolerance
    println((method = name, patches = length(map.patches), maximum_error = maximum(error)))
end
@assert maximum(p -> sup(abs(only(p.error_bounds))), maps[end].patches) <= tolerance
println("Largest guaranteed point-enclosure width: ", maximum(diam, grid_values[end]))

# Existing partitions can be refined without merging their boundaries.
# Always provide the original function so new expansions recover missing terms.
refined = adaptive_map(gaussian, maps[1]; atol = tolerance / 10)
refined_error = maximum(abs(refined([x, y]) - gaussian([x, y])) for x in grid, y in grid)
@assert refined.converged && refined_error < tolerance / 10
println((refined_patches = length(refined.patches), maximum_error = refined_error))

# The upper row shows partitions in physical coordinates. The lower row uses
# the same color scale for numerical discrepancies in all five cases. The last
# column plots interval midpoint discrepancy; its certified error is reported
# separately above. Bounds and sampled discrepancies answer different questions.
fig = Figure(size = (1530, 710), fontsize = 12)
for (column, ((name, estimator, _), map, error)) in enumerate(zip(methods, maps, errors))
    partition_axis = Axis(
        fig[1, column]; xlabel = "x", ylabel = "y",
        title = "$name\n$(length(map.patches)) patches", aspect = DataAspect()
    )
    for patch in map.patches
        lo, hi = estimator isa IntervalBound ? (inf.(domain(patch)), sup.(domain(patch))) : (patch.lower, patch.upper)
        x0, y0 = lo
        x1, y1 = hi
        lines!(partition_axis, [x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0]; color = estimator isa IntervalBound ? :darkgreen : Makie.to_colormap(:tab10)[1], linewidth = 0.7)
    end
    accuracy = Axis(fig[2, column]; xlabel = "x", ylabel = "y", title = estimator isa IntervalBound ? "Grid midpoint discrepancy" : "Independent grid discrepancy", aspect = DataAspect())
    heatmap!(accuracy, grid, grid, log10.(max.(error, 1.0e-12)); colormap = :magma, colorrange = (-12, -4))
end
Colorbar(fig[3, 1:5]; colormap = :magma, limits = (-12, -4), vertical = false, label = "log₁₀ numerical absolute discrepancy")
fig

# ## The new interval method and non-axis-aligned geometry
# All four combinations below approximate the same exp(x+y) on the same box,
# at order three and atol=1/1000. GuardedTail remains heuristic. IntervalBound
# reevaluates the original expression with Taylor-model inputs on every child
# and bounds the retained-coefficient uncertainty plus absolute remainder.
# Choosing a direction affects efficiency, not that acceptance proof.
diagonal(v) = exp(v[1] + v[2])
settings = (
    ("GuardedTail / boxes", GuardedTail(), :tail, false),
    ("GuardedTail / polygons", GuardedTail(), :oriented, false),
    ("IntervalBound / boxes", IntervalBound(), :tail, true),
    ("IntervalBound / polygons", IntervalBound(), :oriented, true),
)
fits = [adaptive_map(diagonal, [-0.5, -0.5], [0.5, 0.5]; estimator, splitter, order = 3, atol = 1 // 1000) for (_, estimator, splitter, _) in settings]
@assert all(m -> m.converged, fits)
for ((label, _, _, certified), fit) in zip(settings, fits)
    indicator = certified ? maximum(p -> sup(abs(only(p.error_bounds))), fit.patches) : maximum(p -> maximum(p.error_estimate), fit.patches)
    println((label, patches = length(fit.patches), acceptance_indicator = indicator, certified))
end

# Upper panels show physical partitions; lower panels show the y=0 slice.
# Blue curves are ordinary polynomial approximations. Green bands enclose the
# original function throughout each complete plotted cell, including remainders.
# Sampling the black reference curve is for illustration only.
comparison_fig = Figure(size = (1280, 710), fontsize = 13)
edges = collect(range(-0.5, 0.5; length = 81))
plot_x = [x for i in 1:(length(edges) - 1) for x in (edges[i], edges[i + 1])]
for (column, ((label, _, splitter, certified), fit)) in enumerate(zip(settings, fits))
    color = certified ? :darkgreen : :dodgerblue
    ax = Axis(comparison_fig[1, column]; xlabel = "physical x", ylabel = "physical y", title = "$label\n$(length(fit.patches)) patches", aspect = DataAspect())
    for patch in fit.patches
        if splitter == :oriented
            poly!(ax, [Point2f(Float64.(v)) for v in polygon_vertices(domain(patch))]; color = (color, 0.12), strokecolor = color, strokewidth = 0.8)
        else
            lo, hi = certified ? (inf.(domain(patch)), sup.(domain(patch))) : (patch.lower, patch.upper)
            poly!(ax, Rect2f(lo[1], lo[2], hi[1] - lo[1], hi[2] - lo[2]); color = (color, 0.12), strokecolor = color, strokewidth = 0.8)
        end
    end
    ax = Axis(comparison_fig[2, column]; xlabel = "physical x (y=0)", ylabel = "function value", title = certified ? "Whole-cell original-function bounds" : "Ordinary polynomial fit")
    if certified
        cells = [enclose(fit, [interval(edges[i], edges[i + 1]), interval(0)]) for i in 1:(length(edges) - 1)]
        @assert all(isguaranteed, cells)
        band!(ax, plot_x, [inf(v) for v in cells for _ in 1:2], [sup(v) for v in cells for _ in 1:2]; color = (:green, 0.35), label = "Certified cell enclosure")
    else
        lines!(ax, edges, [fit([x, 0]) for x in edges]; color, linewidth = 3, label = "Polynomial approximation")
    end
    lines!(ax, edges, exp.(edges); color = :black, linestyle = :dash, label = "Original function (samples)")
    axislegend(ax; position = :lt, labelsize = 10)
end
comparison_fig

# ## What the partition count and acceptance indicator mean
# The shared tolerance has different evidence behind it: blue indicators are
# heuristic; green indicators are uniform error bounds. Patch counts measure
# partition size, not construction time or allocation cost. Separate benchmark
# scripts measure those costs. The C0/C1/C2 and optimization tutorials add smooth
# overlapping surrogates without treating the source tolerance as inherited.
# The logarithmic display places zero indicators at machine epsilon.
labels = ["Guard\nboxes", "Guard\npolygons", "Interval\nboxes", "Interval\npolygons"]
colors = [:dodgerblue, :dodgerblue, :darkgreen, :darkgreen]
summary_fig = Figure(size = (1050, 410), fontsize = 14)
ax = Axis(summary_fig[1, 1]; xticks = (1:4, labels), ylabel = "patch count", title = "Same function, order and tolerance")
barplot!(ax, 1:4, [length(m.patches) for m in fits]; color = colors)
ax = Axis(summary_fig[1, 2]; xticks = (1:4, labels), ylabel = "acceptance indicator", yscale = log10, title = "Bound only for interval methods (green)")
indicators = [certified ? maximum(p -> sup(abs(only(p.error_bounds))), m.patches) : maximum(p -> maximum(p.error_estimate), m.patches) for ((_, _, _, certified), m) in zip(settings, fits)]
scatter!(ax, 1:4, max.(indicators, eps(Float64)); color = colors, markersize = 14)
hlines!(ax, [1 / 1000]; color = :black, linestyle = :dash, label = "Requested atol")
axislegend(ax; position = :lb, labelsize = 11)
summary_fig
