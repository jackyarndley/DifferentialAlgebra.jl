# # Interval evaluation and Taylor-model enclosures
# The same nonlinear expression produces three different kinds of bounds.
# Install the optional IntervalArithmetic dependency in this example environment.

using DifferentialAlgebra, IntervalArithmetic
using CairoMakie

f(x, y) = exp(x + y) * cos(x * y) + log(one(x) + one(x) + x^2)
box = [interval(Float64, -1 // 4, 1 // 4), interval(Float64, -1 // 8, 1 // 8)]

# Ordinary DA computes a truncated polynomial. The interval bound below encloses
# exactly that stored polynomial, including its stored Float64 coefficients.
# It contains no proof of the missing terms of `f`.
x, y = variables((:x, :y); order = 3)
p = f(x, y)
stored_range = enclose(p, box)
println("Stored polynomial: ", stored_range)
ordinary_fit = compile(p)
edges = range(inf(box[1]), sup(box[1]); length = 161)
slice_boxes = [[interval(edges[i], edges[i + 1]), interval(0)] for i in 1:(length(edges) - 1)]
stored_slice = [enclose(p, b) for b in slice_boxes]

# Direct interval evaluation encloses the original expression, but repeated
# occurrences of a variable lose dependence information.
direct_range = f(box...)
println("Direct interval function evaluation: ", direct_range)

# Reevaluating the original expression with model inputs propagates a rigorously
# justified absolute remainder at every operation. The normalized polynomial
# and its remainder jointly enclose the original function throughout the box.
x, y = taylor_models(box; order = 3, names = (:x, :y))
m = f(x, y)
println("Taylor model: ", enclose(m))
println("Absolute remainder: ", remainder(m))
println("At the physical center: ", evaluate(m, [0, 0]))
println("On a physical subbox: ", enclose(m, [interval(Float64, 0, 1 // 4), interval(0)]))

# Widths describe different objects; a narrower stored-polynomial enclosure does
# not imply a more accurate enclosure of the original function.
println(
    "Widths (polynomial / direct / model): ",
    (diam(stored_range), diam(direct_range), diam(enclose(m)))
)

# # Certified fitting with automatic box splitting
# Each child reevaluates `f` with fresh model inputs. Acceptance bounds the error
# relative to the retained midpoint polynomial, including coefficient widths.
# Restricting `m` itself would keep its original remainder.
parent = compile(m)
tolerance = 1.0e-5
root = validated_adaptive_map(v -> f(v...), box; order = 3, atol = tolerance, max_depth = 0, strict = false)
ads = validated_adaptive_map(v -> f(v...), box; order = 3, atol = tolerance)
@assert ads.converged
println(
    (
        patches = length(ads.patches), tolerance,
        parent_error_bound = sup(abs(only(only(root.patches).error_bounds))),
        largest_child_error_bound = maximum(p -> sup(abs(only(p.error_bounds))), ads.patches),
    )
)

# These bands bound entire physical cells along y=0. The black curve samples
# the original expression for orientation only; it is not an inclusion oracle.
parent_slice = [enclose(parent, b) for b in slice_boxes]
ads_slice = [enclose(ads, b) for b in slice_boxes]
direct_slice = [f(b...) for b in slice_boxes]
plot_x = [x for i in 1:(length(edges) - 1) for x in (edges[i], edges[i + 1])]
lower_band(values) = [inf(v) for v in values for _ in 1:2]
upper_band(values) = [sup(v) for v in values for _ in 1:2]
centers = [(edges[i] + edges[i + 1]) / 2 for i in 1:(length(edges) - 1)]

fig = Figure(size = (1100, 790), fontsize = 14)
fit_axis = Axis(fig[1, 1]; xlabel = "physical x (y=0)", ylabel = "function value", title = "Uniform interval bands on each plotted cell")
band!(fit_axis, plot_x, lower_band(parent_slice), upper_band(parent_slice); color = (:orange, 0.45), label = "Parent Taylor model")
band!(fit_axis, plot_x, lower_band(ads_slice), upper_band(ads_slice); color = (:dodgerblue, 0.5), label = "Certified box ADS")
lines!(fit_axis, edges, f.(edges, 0); color = :black, label = "Original expression (samples)")
lines!(fit_axis, edges, [only(evaluate(ordinary_fit, [x, 0])) for x in edges]; color = :red, linestyle = :dash, label = "Ordinary DA polynomial")
axislegend(fit_axis; position = :lt, labelsize = 11)

width_axis = Axis(fig[1, 2]; xlabel = "physical x (y=0)", ylabel = "interval width", title = "Cell enclosure widths (different contracts)", yscale = log10)
for (label, values, color) in (
        ("Stored polynomial only", stored_slice, :red), ("Direct interval function", direct_slice, :gray),
        ("Parent Taylor model", parent_slice, :orange), ("Certified box ADS", ads_slice, :dodgerblue),
    )
    lines!(width_axis, centers, diam.(values); label, color)
end
axislegend(width_axis; position = :lt, labelsize = 11)

partition_axis = Axis(fig[2, 1]; xlabel = "physical x", ylabel = "physical y", title = "$(length(ads.patches)) certified child boxes", aspect = DataAspect())
for patch in ads.patches
    bx, by = domain(patch)
    poly!(partition_axis, Rect2f(inf(bx), inf(by), diam(bx), diam(by)); color = :lightblue, strokecolor = :steelblue, strokewidth = 1)
end

error_axis = Axis(fig[2, 2]; xlabel = "child patch", ylabel = "uniform absolute error bound", yscale = log10, title = "Split criterion includes coefficient uncertainty")
child_errors = [sup(abs(only(p.error_bounds))) for p in ads.patches]
scatter!(error_axis, eachindex(child_errors), child_errors; color = :dodgerblue, label = "Child error bounds")
hlines!(error_axis, [sup(abs(only(only(root.patches).error_bounds)))]; color = :orange, label = "Unsplit parent")
hlines!(error_axis, [tolerance]; color = :black, linestyle = :dash, label = "Requested tolerance")
axislegend(error_axis; position = :rt, labelsize = 11)
fig

# The smaller errors are easier to see on their own scale. These intervals
# bound f minus each patch's retained midpoint polynomial uniformly throughout
# its complete box; they are not errors inferred from sampled discrepancies.
error_fig = Figure(size = (1000, 380), fontsize = 14)
parent_error_axis = Axis(error_fig[1, 1]; xlabel = "physical x", ylabel = "uniform fit-error interval", title = "Unsplit parent (whole box)")
parent_error = only(only(root.patches).error_bounds)
band!(parent_error_axis, [inf(box[1]), sup(box[1])], fill(inf(parent_error), 2), fill(sup(parent_error), 2); color = (:orange, 0.4))
hlines!(parent_error_axis, [-tolerance, tolerance]; color = :black, linestyle = :dash)
child_error_axis = Axis(error_fig[1, 2]; xlabel = "child patch", ylabel = "uniform fit-error interval", title = "Fresh child expansions (different local fits)")
rangebars!(
    child_error_axis, 1:length(ads.patches), [inf(only(p.error_bounds)) for p in ads.patches],
    [sup(only(p.error_bounds)) for p in ads.patches]; color = :dodgerblue, linewidth = 3, whiskerwidth = 8
)
hlines!(child_error_axis, [-tolerance, tolerance]; color = :black, linestyle = :dash)
error_fig

# The standalone script also writes a shareable PNG. Literate embeds the figure.
if abspath(PROGRAM_FILE) == @__FILE__
    output = joinpath(@__DIR__, "..", "results", "interval_models.png")
    mkpath(dirname(output))
    save(output, fig; px_per_unit = 2)
    println("Saved interval plot: ", abspath(output))
    error_output = joinpath(dirname(output), "interval_fit_errors.png")
    save(error_output, error_fig; px_per_unit = 2)
    println("Saved fit-error plot: ", abspath(error_output))
end

# A positive constant term does not establish logarithm validity on a whole box.
# This fails because `1+x` includes negative values and zero over [-2,2].
wide, = taylor_models([interval(-2, 2)]; order = 3)
try
    log(1 + wide)
catch err
    err isa DomainError || rethrow()
    println("Rejected invalid logarithm domain: ", err.msg)
end

# Floating literals denote stored binary values. For a desired decimal real,
# use a rational or string interval such as `interval(1//10)` or `I"0.1"`.
# GuardedTail/ExtrapolatedTail/LastTerms and Taylor time integration remain heuristic.
# The certified box API above validates uncertainty-domain approximation only.
