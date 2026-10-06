# # C0, C1 and C2 continuity across ADS fits
# Independent ADS expansions generally have jumps at their shared faces.
# `continuous_map` is an optional partition-of-unity approximation: reevaluate
# the original function on overlapping covers, then blend the retained fits.
# The linear/cubic/quintic taper endpoints prove C0/C1/C2 regularity respectively;
# the plots illustrate their effects rather than proving continuity.

using DifferentialAlgebra, ForwardDiff, IntervalArithmetic, CairoMakie

f(v) = v[1]^3
raw = adaptive_map(f, [-1.0], [1.0]; order = 2, atol = 1 // 2)
@assert length(raw.patches) == 2
surrogates = [continuous_map(f, raw; continuity = k) for k in (:c0, :c1, :c2)]
println("Raw patches: ", length(raw.patches))
println("Recomputed overlap estimates: ", [s.error_estimate for s in surrogates])

# The two core boxes meet at x=0. An overlap of 1/4 extends their supports
# to x=1/4 and x=-1/4. Differentiability also matters at these support edges.
# Quadratic fits of a cubic deliberately expose the mismatch of local jets.
# Blending is not interpolation and can change derivative accuracy.
jet(g, x, d) = d == 0 ? g(x) : ForwardDiff.derivative(t -> jet(g, t, d - 1), x)
xs = sort!(unique(vcat(collect(range(-0.6, 0.6; length = 481)), [x for edge in (-0.25, 0.0, 0.25) for x in range(edge - 0.025, edge + 0.025; length = 201)])))
fig = Figure(size = (1240, 800), fontsize = 13)
for (column, label, fit) in zip(1:4, ("Independent fits", "C0 blend", "C1 blend", "C2 blend"), (raw, surrogates...))
    g(x) = fit([x])
    for derivative in 0:2
        ax = Axis(fig[derivative + 1, column]; xlabel = derivative == 2 ? "physical x" : "", ylabel = column == 1 ? ("Value", "First derivative", "Second derivative")[derivative + 1] : "", title = derivative == 0 ? label : "")
        lines!(ax, xs, [jet(x -> x^3, x, derivative) for x in xs]; color = :black, linestyle = :dash, label = "Original function")
        ## Draw separate segments at the interfaces so the line plot displays
        ## one-sided jumps rather than connecting them as if they were smooth.
        if derivative <= column - 2
            lines!(ax, xs, [jet(g, x, derivative) for x in xs]; color = :dodgerblue, label = "ADS surrogate")
        else
            for (j, (lo, hi)) in enumerate(((-0.6, -0.25), (-0.25, 0.0), (0.0, 0.25), (0.25, 0.6)))
                segment = filter(x -> lo < x < hi, xs)
                lines!(ax, segment, [jet(g, x, derivative) for x in segment]; color = :dodgerblue, label = j == 1 ? "ADS surrogate" : nothing)
            end
        end
        vlines!(ax, [-0.25, 0.0, 0.25]; color = (:gray, 0.5), linestyle = :dot)
        derivative == 0 && column == 1 && axislegend(ax; position = :lt, labelsize = 10)
    end
end
fig

# # Interval ADS and a continuous numeric surrogate
# IntervalBound remains an ADS error method. Its fresh overlap models retain
# their coefficients and absolute remainder. Midpoint polynomials define the
# explicitly requested numeric surrogate only; no model arithmetic uses them.
quadratic(v) = v[1]^2
certified = adaptive_map(quadratic, [interval(-1, 1)]; estimator = IntervalBound(), order = 1, atol = 1 // 64)
smooth = continuous_map(quadratic, certified; continuity = :c2, atol = 1 // 16)
println("Core error bound: ", maximum(p -> sup(abs(only(p.error_bounds))), certified.patches))
println("Blended surrogate uniform error interval: ", only(smooth.error_bounds))
@assert sup(abs(only(smooth.error_bounds))) <= 1 // 16

# Enlarged fits have a new accuracy contract: the original patch tolerance is
# not inherited. For x², Taylor's exact identity f-P_c=(x-c)² proves the error
# interval [0,r²] for each cover; nonnegative normalized blending preserves the
# hull of those intervals. The band below instead encloses the original
# function throughout each complete plotted cell, using `enclose(smooth, cell)`.
edges = collect(range(-1, 1; length = 129))
values = [enclose(smooth, [interval(edges[i], edges[i + 1])]) for i in 1:(length(edges) - 1)]
plot_x = [x for i in 1:(length(edges) - 1) for x in (edges[i], edges[i + 1])]
band_fig = Figure(size = (1040, 400), fontsize = 14)
ax = Axis(band_fig[1, 1]; xlabel = "physical x", ylabel = "function value", title = "Original function enclosed on complete cells")
band!(ax, plot_x, [inf(v) for v in values for _ in 1:2], [sup(v) for v in values for _ in 1:2]; color = (:dodgerblue, 0.35), label = "Original-function enclosure")
lines!(ax, edges, edges .^ 2; color = :black, label = "Original function (samples)")
lines!(ax, edges, [smooth([x]) for x in edges]; color = :orange, label = "C2 numeric surrogate")
axislegend(ax; position = :ct, labelsize = 11)
ax = Axis(band_fig[1, 2]; xlabel = "core patch / overlap fit", ylabel = "uniform error interval", title = "Overlap error is checked again")
errors = [only(p.error_bounds) for p in certified.patches]
rangebars!(ax, eachindex(errors), inf.(errors), sup.(errors); color = :dodgerblue, linewidth = 3, whiskerwidth = 8, label = "Original core fits")
hlines!(ax, [inf(only(smooth.error_bounds)), sup(only(smooth.error_bounds))]; color = :orange, label = "Blended uniform bound")
axislegend(ax; position = :ct, labelsize = 11)
band_fig

# Numeric ForwardDiff gradients/Hessians belong to the surrogate. Its function
# error bound does not certify them or a minimizer. Enclosure endpoints are
# interval bounds for the original function and need not themselves be smooth.
println("Surrogate value / gradient / Hessian at x=1/8: ", (smooth([1 / 8]), jet(x -> smooth([x]), 1 / 8, 1), jet(x -> smooth([x]), 1 / 8, 2)))
try
    smooth([2.0])
catch err
    err isa DomainError || rethrow()
    println("Rejected outside-domain objective query: ", err.msg)
end

# The standalone script saves the same figures embedded by Literate.
if abspath(PROGRAM_FILE) == @__FILE__
    directory = joinpath(@__DIR__, "..", "results")
    mkpath(directory)
    for (name, figure) in (("ads_continuity.png", fig), ("ads_continuity_intervals.png", band_fig))
        path = joinpath(directory, name)
        save(path, figure; px_per_unit = 2)
        println("Saved continuity plot: ", abspath(path))
    end
end
