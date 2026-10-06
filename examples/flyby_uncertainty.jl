# # Mars flyby uncertainty on the B-plane
# The B-plane is perpendicular to the incoming hyperbolic asymptote; see
# [JPL's definition](https://cneos.jpl.nasa.gov/glossary/b_plane.html).
# A point-mass flyby turns an incoming velocity into a curved outgoing cone.
# Compare all ADS estimators and bound periapsis altitude over complete leaves.
# These are explicit patched-conic formulas, not a verified encounter integrator.
using DifferentialAlgebra, IntervalArithmetic, CairoMakie, LinearAlgebra, Random

# Use illustrative Mars parameters μ=42,828 km³/s², v∞=6 km/s and R=3,390 km.
# The length unit μ/v∞² is exactly 1,189 km for these declared integer inputs.
# b²=B_T²+B_R² in that unit, tan(δ/2)=1/b, and conservation gives
# b²=rp²+2rp. Eliminating atan and the unit B vector leaves only rational
# arithmetic and sqrt, all supported by the native Taylor-model layer.
function flyby(v)
    BT, BR = v
    unit = one(BT)
    two = unit + unit
    b2 = BT^2 + BR^2
    denominator = unit + b2
    return [(b2 - unit) / denominator, -two * BT / denominator, -two * BR / denominator, sqrt(unit + b2) - unit]
end
length_unit, radius = 1189, 3390
lower, upper = [3.0, 1.0], [6.0, 4.0]
order, tolerance = 4, 1 // 1000
# The same tolerance applies to each dimensionless output: 6 m/s for velocity
# and 1.189 km for periapsis. Use an atol vector for different physical targets.
methods = (("GuardedTail", GuardedTail()), ("ExtrapolatedTail", ExtrapolatedTail()), ("LastTerms", LastTerms()), ("IntervalBound", IntervalBound()))
fits = [adaptive_map(flyby, lower, upper; estimator, order, atol = tolerance) for (_, estimator) in methods]
oriented = adaptive_map(flyby, lower, upper; estimator = IntervalBound(), order, atol = tolerance, splitter = :oriented, directions = [1 1; -1 1])
@assert all(m -> m.converged, fits) && oriented.converged
for ((label, _), fit) in zip(methods, fits)
    println((label, leaves = length(fit.patches)))
end
println("IntervalBound / diagonal polygons: ", length(oriented.patches), " leaves")

# An independent exact oracle: at (B_T,B_R)=(2,2), b²=8, rp=2 and
# vout/v∞=(7/9,-4/9,-4/9). The speed identity holds symbolically because
# (b²-1)²+4b²=(b²+1)². This point is a formula check, not part of our prior.
oracle = flyby([interval(2), interval(2)])
@assert all(isguaranteed, oracle)
@assert all(in_interval(q, value) for (q, value) in zip((7 // 9, -4 // 9, -4 // 9, 2), oracle))
rng = Xoshiro(2026)
points = [[3 + 3rand(rng), 1 + 3rand(rng)] for _ in 1:1200]
states = flyby.(points)
@assert all(s -> abs(sum(abs2, s[1:3]) - 1) < 1.0e-14, states)

# ## Compare the same B-plane prior and error tolerance
# Only IntervalBound supplies a uniform fit-error certificate. The first three
# methods produce approximate numbers; all four IntervalBound outputs include
# absolute remainders. The velocity projections below illustrate the original
# function, with each input leaf boundary mapped into the outgoing cone.
methods_fig = Figure(size = (1370, 730), fontsize = 13)
for (column, ((label, _), fit)) in enumerate(zip(methods, fits))
    ax = Axis(methods_fig[1, column]; xlabel = "B_T / 1,189 km", ylabel = "B_R / 1,189 km", title = "$label\n$(length(fit.patches)) leaves", aspect = DataAspect())
    image_axis = Axis(methods_fig[2, column]; xlabel = "v_y / v∞", ylabel = "v_z / v∞", title = "Mapped boundary curves (samples)", aspect = DataAspect())
    for patch in fit.patches
        lo, hi = column == 4 ? (inf.(domain(patch)), sup.(domain(patch))) : (patch.lower, patch.upper)
        color = column == 4 ? :darkgreen : :steelblue
        poly!(ax, Rect2f(lo[1], lo[2], hi[1] - lo[1], hi[2] - lo[2]); color = (color, 0.12), strokecolor = color, strokewidth = 0.8)
        t = range(0, 1; length = 13)
        boundary = vcat([[lo[1] + q * (hi[1] - lo[1]), lo[2]] for q in t], [[hi[1], lo[2] + q * (hi[2] - lo[2])] for q in t], [[hi[1] - q * (hi[1] - lo[1]), hi[2]] for q in t], [[lo[1], hi[2] - q * (hi[2] - lo[2])] for q in t])
        values = flyby.(boundary)
        lines!(image_axis, getindex.(values, 2), getindex.(values, 3); color = (color, 0.55), linewidth = 0.8)
    end
end
methods_fig

# ## Whole-polygon altitude classification and the outgoing cone
# Each classification uses the whole leaf enclosure, not a point or midpoint.
# Green is above the surface, red is below, and orange is unresolved. A red
# region represents a collision in this finite-radius planet model: the virtual
# point-mass outgoing asymptote is not an actual post-impact spacecraft path.
# The drawn critical circle is a rounded visualization of bcrit²=(1+R/L)²-1.
classification_fig = Figure(size = (1280, 520), fontsize = 14)
ax = Axis(classification_fig[1, 1]; xlabel = "B_T / 1,189 km", ylabel = "B_R / 1,189 km", title = "Certified altitude classification on polygons", aspect = DataAspect())
counts = [0, 0, 0]
for patch in oriented.patches
    altitude = enclose(oriented, domain(patch))[4] * interval(length_unit) - interval(radius)
    @assert isguaranteed(altitude)
    status = inf(altitude) > 0 ? 1 : sup(altitude) < 0 ? 2 : 3
    counts[status] += 1
    color = (:seagreen, :firebrick, :goldenrod)[status]
    poly!(ax, [Point2f(Float64.(v)) for v in polygon_vertices(domain(patch))]; color = (color, 0.35), strokecolor = color, strokewidth = 0.6)
end
critical = sqrt((1 + radius / length_unit)^2 - 1)
angles = range(0, π / 2; length = 181)
lines!(ax, critical .* cos.(angles), critical .* sin.(angles); color = :black, linestyle = :dash, linewidth = 2, label = "Zero-altitude boundary (rounded)")
xlims!(ax, lower[1], upper[1]); ylims!(ax, lower[2], upper[2])
axislegend(ax, [PolyElement(color = :seagreen), PolyElement(color = :firebrick), PolyElement(color = :goldenrod), LineElement(color = :black, linestyle = :dash)], ["Above surface", "Below surface", "Unresolved cell", "Zero-altitude boundary (rounded)"]; position = :rt, labelsize = 10)
ax = Axis3(classification_fig[1, 2]; xlabel = "v_x / v∞", ylabel = "v_y / v∞", zlabel = "v_z / v∞", title = "Point-mass outgoing asymptotes", aspect = :data, azimuth = 0.7, elevation = 0.35)
altitudes = [(s[4] * length_unit - radius) for s in states]
cloud = scatter!(ax, first.(states), getindex.(states, 2), getindex.(states, 3); color = altitudes, colormap = :balance, colorrange = (-2000, 8000), markersize = 4)
for s in states[1:80:end]
    lines!(ax, [0, s[1]], [0, s[2]], [0, s[3]]; color = (:gray, 0.15))
end
Colorbar(classification_fig[1, 3], cloud; label = "Virtual periapsis altitude (km)")
classification_fig

# Counts distinguish a proven whole-leaf sign from cells whose enclosures still
# cross zero. Splitting them more finely can resolve some, but not all, cells
# at a physical boundary. The count alone is not a probability of impact.
println("Guaranteed above-surface / below-surface / unresolved leaves: ", counts)

# ## Interval altitude bands and the price of validation
# Cell enclosures include physical B_T variation, coefficient rounding and
# approximation remainder. The sampled scalar curve is supplementary. Leaf
# counts are reported separately from construction timing and allocations.
bands_fig = Figure(size = (1180, 440), fontsize = 14)
ax = Axis(bands_fig[1, 1]; xlabel = "B_T / 1,189 km (B_R=2)", ylabel = "periapsis altitude (km)", title = "Original altitude enclosed on complete cells")
edges = collect(range(3, 6; length = 65))
values = [enclose(oriented, [interval(edges[j], edges[j + 1]), interval(2)])[4] * interval(length_unit) - interval(radius) for j in 1:64]
@assert all(isguaranteed, values)
xs = [x for j in 1:64 for x in (edges[j], edges[j + 1])]
band!(ax, xs, [inf(v) for v in values for _ in 1:2], [sup(v) for v in values for _ in 1:2]; color = (:green, 0.4), label = "Whole-cell altitude enclosure")
lines!(ax, edges, [flyby([b, 2.0])[4] * length_unit - radius for b in edges]; color = :black, label = "Original function (samples)")
hlines!(ax, [0]; color = :firebrick, linestyle = :dash, label = "Planet surface")
axislegend(ax; position = :lt, labelsize = 10)
labels = ["Guarded", "Decay", "LastTerms", "Interval\nboxes", "Interval\npolygons"]
ax = Axis(bands_fig[1, 2]; xticks = (1:5, labels), ylabel = "leaf count", title = "Same function, domain, order and tolerance")
barplot!(ax, 1:5, [length(f.patches) for f in vcat(fits, [oriented])]; color = [:steelblue, :steelblue, :steelblue, :darkgreen, :darkgreen])
bands_fig
