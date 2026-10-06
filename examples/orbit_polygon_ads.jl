# # Polygon ADS for an inclined eccentric orbit
# Uncertain eccentric anomaly E and argument of periapsis ω form a curved
# position/velocity ribbon. Compare all four ADS estimators on the same square,
# then use a hexagonal prior to express correlated orbital-phase uncertainty.
# This explicit element-to-state map supports IntervalBound. E is an input;
# no numerical Kepler solve or time integrator is being certified.
using DifferentialAlgebra, IntervalArithmetic, CairoMakie, LinearAlgebra, Random

# The rotation is Rz(Ω) Rx(i) Rz(ω). Length and speed are normalized by a and
# sqrt(μ/a). The fixed e=1/5, i=1/2 and Ω=1/3 are exact rational inputs.
# Lift constants into the input algebra before taking their elementary functions.
# The element-to-state convention is described by [JPL](https://ssd.jpl.nasa.gov/planets/approx_pos.html).
function orbit_state(v)
    E, ω = v
    e, i, Ω = zero(E) + 1 // 5, zero(E) + 1 // 2, zero(E) + 1 // 3
    s, c = sincos(E)
    β = sqrt(1 - e^2)
    radius = 1 - e * c
    x, y, vx, vy = c - e, β * s, -s / radius, β * c / radius
    sw, cw = sincos(ω)
    si, ci = sincos(i)
    sn, cn = sincos(Ω)
    rotate(x, y) = [cn * (cw * x - sw * y) - sn * ci * (sw * x + cw * y), sn * (cw * x - sw * y) + cn * ci * (sw * x + cw * y), si * (sw * x + cw * y)]
    return vcat(rotate(x, y), rotate(vx, vy))
end
lower, upper = [0.0, -0.5], [1.0, 0.5]
order, tolerance = 3, 1 // 1000
methods = (("GuardedTail", GuardedTail()), ("ExtrapolatedTail", ExtrapolatedTail()), ("LastTerms", LastTerms()), ("IntervalBound", IntervalBound()))
boxes = [adaptive_map(orbit_state, lower, upper; estimator, order, atol = tolerance) for (_, estimator) in methods]
polygons = [adaptive_map(orbit_state, lower, upper; estimator, order, atol = tolerance, splitter = :oriented, directions = [1 1; -1 1]) for (_, estimator) in methods]
@assert all(m -> m.converged, vcat(boxes, polygons))
for ((label, _), box, polygon) in zip(methods, boxes, polygons)
    println((label, box_leaves = length(box.patches), polygon_leaves = length(polygon.patches)))
end
@assert all(p -> maximum(sup.(abs.(p.error_bounds))) <= tolerance, boxes[4].patches)
@assert all(p -> maximum(sup.(abs.(p.error_bounds))) <= tolerance, polygons[4].patches)

# ## Estimator versus splitting frame
# Each column changes the error method; each row changes the geometry. The
# diagonal frame follows E+ω, the dominant angular combination near a circular
# orbit. The eccentric orbit still has dependence on the transverse coordinate.
# All eight cases use the same function, order, output scales and tolerance.
# Diagonal frames reduce the heuristic leaf counts here but increase the
# IntervalBound count: a parallelogram cover and conservative composition
# remainders can outweigh alignment. Polygon geometry is not an automatic
# performance improvement; compare the measured partition and construction cost.
partition_fig = Figure(size = (1400, 730), fontsize = 14)
for (row, maps, geometry) in ((1, boxes, "Boxes"), (2, polygons, "Diagonal polygons"))
    for (column, ((label, _), fit)) in enumerate(zip(methods, maps))
        ax = Axis(partition_fig[row, column]; xlabel = "E (rad)", ylabel = "ω (rad)", title = "$geometry / $label\n$(length(fit.patches)) leaves", aspect = DataAspect())
        for patch in fit.patches
            vertices = if row == 2
                polygon_vertices(domain(patch))
            else
                lo, hi = column == 4 ? (inf.(domain(patch)), sup.(domain(patch))) : (patch.lower, patch.upper)
                ([lo[1], lo[2]], [hi[1], lo[2]], [hi[1], hi[2]], [lo[1], hi[2]])
            end
            color = column == 4 ? :darkgreen : :steelblue
            poly!(ax, [Point2f(Float64.(v)) for v in vertices]; color = (color, 0.12), strokecolor = color, strokewidth = 0.8)
        end
    end
end
partition_fig

# ## A correlated prior and its orbital ribbon
# Relative to (E,ω)=(1/2,0), intersect the square with |δE+δω|≤3/4.
# Every vertex is dyadic, so the polygon and its clipped leaf areas are exact.
# A smaller prior is a different physical uncertainty set, not an efficiency
# comparison with the square above. Every new leaf reevaluates the original map.
prior = ConvexPolygon([(0, -1 // 4), (1 // 4, -1 // 2), (1, -1 // 2), (1, 1 // 4), (3 // 4, 1 // 2), (0, 1 // 2)])
correlated = adaptive_map(orbit_state, prior; estimator = IntervalBound(), order, atol = tolerance, directions = [1 1; -1 1])
@assert sum(p -> domain_area(domain(p)), correlated.patches) == domain_area(prior)
rng = Xoshiro(2026)
points = [[rand(rng), rand(rng) - 0.5] for _ in 1:2000]
points = filter(v -> abs(v[1] - 0.5 + v[2]) <= 0.75, points)
states = orbit_state.(points)
@assert all(u -> abs(sum(abs2, u[4:6]) / 2 - 1 / norm(u[1:3]) + 1 / 2) < 1.0e-12, states)
@assert all(u -> abs(sum(abs2, cross(u[1:3], u[4:6])) - (1 - (1 / 5)^2)) < 1.0e-12, states)
println("Hexagonal prior: ", length(correlated.patches), " leaves; full state enclosure: ", enclose(correlated))

# Plot a=10,000 km around a representative Earth of radius 6,378 km. The
# sampled ribbon illustrates the geometry; green rectangles bound the entire
# image of each polygon leaf. Rotation and the six state outputs remain in the
# certified map; the 3D point cloud is not used as an inclusion oracle.
ribbon_fig = Figure(size = (1350, 570), fontsize = 14)
ax = Axis3(ribbon_fig[1, 1]; xlabel = "x (10³ km)", ylabel = "y (10³ km)", zlabel = "z (10³ km)", title = "Inclined orbit and correlated uncertainty", aspect = :data, azimuth = -1.1, elevation = 0.9)
mesh!(ax, Sphere(Point3f(0), 6.378f0); color = :steelblue)
nominal_orbit = orbit_state.([[E, 0.0] for E in range(0, 2π; length = 241)])
lines!(ax, 10first.(nominal_orbit), 10getindex.(nominal_orbit, 2), 10getindex.(nominal_orbit, 3); color = :gray, linewidth = 2)
scatter!(ax, 10first.(states), 10getindex.(states, 2), 10getindex.(states, 3); color = first.(points), colormap = :plasma, colorrange = (0, 1), markersize = 3)
Colorbar(ribbon_fig[2, 1]; colormap = :plasma, limits = (0, 1), vertical = false, label = "Eccentric anomaly E (rad)")
ax = Axis(ribbon_fig[1, 2]; xlabel = "x (10³ km)", ylabel = "y (10³ km)", title = "Whole-leaf position bounds in projection", aspect = DataAspect())
for patch in correlated.patches
    value = enclose(correlated, domain(patch))[1:2] .* interval(10)
    @assert all(isguaranteed, value)
    poly!(ax, Rect2f(inf(value[1]), inf(value[2]), diam(value[1]), diam(value[2])); color = (:green, 0.06), strokecolor = (:darkgreen, 0.35), strokewidth = 0.6)
end
scatter!(ax, 10first.(states), 10getindex.(states, 2); color = (:black, 0.25), markersize = 2)
ribbon_fig

# ## Error evidence and whole-cell enclosures
# Blue indicators are heuristic. Green indicators are uniform original-function
# fit-error bounds, including coefficient widths and absolute remainders.
# Cell bounds include variation within the cell as well as that approximation
# uncertainty. Restricting a query does not invent a smaller remainder.
summary_fig = Figure(size = (1200, 440), fontsize = 14)
ax = Axis(summary_fig[1, 1]; xticks = (1:4, first.(collect(methods))), ylabel = "leaves", title = "Same square: frame and method costs")
barplot!(ax, (1:4) .- 0.18, length.(getproperty.(boxes, :patches)); width = 0.32, color = :steelblue, label = "Boxes")
barplot!(ax, (1:4) .+ 0.18, length.(getproperty.(polygons, :patches)); width = 0.32, color = :darkorange, label = "Diagonal polygons")
axislegend(ax; position = :rt, labelsize = 10)
ax = Axis(summary_fig[1, 2]; xlabel = "E (rad), ω=0", ylabel = "y (10³ km)", title = "Original position on complete cells")
edges = collect(range(0, 1; length = 65))
cells = [enclose(correlated, [interval(edges[j], edges[j + 1]), interval(0)])[2] * interval(10) for j in 1:64]
@assert all(isguaranteed, cells)
xs = [x for j in 1:64 for x in (edges[j], edges[j + 1])]
band!(ax, xs, [inf(v) for v in cells for _ in 1:2], [sup(v) for v in cells for _ in 1:2]; color = (:green, 0.35), label = "Whole-cell interval enclosure")
lines!(ax, edges, [10orbit_state([E, 0.0])[2] for E in edges]; color = :black, label = "Original function (samples)")
axislegend(ax; position = :lt, labelsize = 10)
summary_fig

# The two-body state identities justify the physical interpretation. Inclusion
# comes from Taylor-model arithmetic, not those rounded conservation checks.
# Construction time and allocations can differ even when there are fewer leaves.
