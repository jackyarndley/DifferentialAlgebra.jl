# # Six-variable orbital-element uncertainty
# Map uncertain (a,e,i,Ω,ω,E) to an inertial Cartesian state at one epoch.
# Compare first- and third-order DA propagation, analytic polynomial moments,
# direct interval evaluation and native Taylor-model enclosures. All six
# coordinates are independent uniforms on [-1,1]; E is supplied, not solved.
using DifferentialAlgebra, IntervalArithmetic, CairoMakie, LinearAlgebra, Random

# JPL gives the [orbital-plane coordinates and inertial rotation](https://ssd.jpl.nasa.gov/planets/approx_pos.html).
# Units are kilometres, seconds and radians. Constants below are declared exact
# rationals, including the illustrative Earth μ=398600.4418 km³/s². Wrap them
# explicitly for direct interval evaluation; native polynomial/model operations
# already handle exact constants. Ordinary scalar results are rounded checks.
function cartesian(ξ)
    c(x) = ξ[1] isa Interval ? interval(Float64, x) : zero(ξ[1]) + x
    a = c(10000) + c(100) * ξ[1]
    e = c(1 // 4) + c(1 // 50) * ξ[2]
    i = c(1 // 2) + c(1 // 100) * ξ[3]
    Ω = c(1 // 3) + c(1 // 100) * ξ[4]
    ω = c(1 // 5) + c(1 // 50) * ξ[5]
    E = c(1) + c(1 // 25) * ξ[6]
    s, co = sincos(E)
    β = sqrt(c(1) - e^2)
    speed = sqrt(c(3986004418 // 10000) / a) / (c(1) - e * co)
    sw, cw = sincos(ω)
    si, ci = sincos(i)
    sn, cn = sincos(Ω)
    rotate(x, y) = [cn * (cw * x - sw * y) - sn * ci * (sw * x + cw * y), sn * (cw * x - sw * y) + cn * ci * (sw * x + cw * y), si * (sw * x + cw * y)]
    return vcat(rotate(a * (co - e), a * β * s), rotate(-speed * s, speed * β * co))
end
ξ = variables((:δa, :δe, :δi, :δΩ, :δω, :δE); order = 3)
state = cartesian(ξ)
linear, cubic = CompiledMap(DifferentialAlgebra.trim.(state, 0, 1)), CompiledMap(state)
box = fill(interval(Float64, -1, 1), 6)
stored_ranges = [enclose(p, box) for p in state]

# ## Moments of the retained polynomial
# E[ξ^k]=0 for odd k and 1/(k+1) for even k. Sum complete pairs of monomials
# for covariance: multiplying in the order-three algebra would discard degrees
# four through six and lose part of the second moment. These are moments of
# the stored polynomial, not a certified distribution of the original function.
uniform_moment(α) = any(isodd, α) ? 0.0 : prod(k -> 1 / (k + 1), α)
function uniform_statistics(polys)
    mean = [sum(m.coefficient * uniform_moment(m.exponents) for m in monomials(p); init = 0.0) for p in polys]
    terms = monomials.(polys .- mean)
    covariance = [sum(a.coefficient * b.coefficient * uniform_moment(a.exponents + b.exponents) for a in terms[i], b in terms[j]; init = 0.0) for i in eachindex(polys), j in eachindex(polys)]
    return mean, covariance
end
polynomial_mean, polynomial_covariance = uniform_statistics(state)
@assert minimum(eigvals(Symmetric(polynomial_covariance))) > -1.0e-10
println("Retained-polynomial mean shift from the nominal state: ", polynomial_mean - cartesian(zeros(6)))
println("Retained-polynomial position standard deviations (km): ", sqrt.(diag(polynomial_covariance)[1:3]))

# A four-point Gauss rule in each coordinate integrates degree seven. Thus the
# 4^6 tensor rule integrates the cubic polynomial and its squared components
# exactly in exact arithmetic. This floating implementation is a separate
# numerical check of the monomial moment calculation, not an inclusion proof.
rule = eigen(SymTridiagonal(zeros(4), [k / sqrt(4k^2 - 1) for k in 1:3]))
nodes, weights = rule.values, rule.vectors[1, :] .^ 2
quadrature_points = [collect(nodes[collect(index)]) for index in Iterators.product(ntuple(_ -> 1:4, 6)...)]
quadrature_weights = [prod(weights[collect(index)]) for index in Iterators.product(ntuple(_ -> 1:4, 6)...)]
quadrature_values = stack(cubic(x) for x in vec(quadrature_points))
quadrature_mean = quadrature_values * vec(quadrature_weights)
centered = quadrature_values .- quadrature_mean
quadrature_covariance = (centered .* transpose(vec(quadrature_weights))) * transpose(centered)
@assert norm(quadrature_mean - polynomial_mean, Inf) < 5.0e-9
@assert norm(quadrature_covariance - polynomial_covariance) / norm(polynomial_covariance) < 1.0e-10

# ## Original-function bounds, with the remainder retained
# Physical assumptions hold over the entire declared parameter set: a>0 and
# 0<e<1, so sqrt(1-e²) and 1-e*cos(E) remain valid. All enclosures below use
# outward rounding. Only the Taylor models enclose the original function's
# polynomial approximation error. Stored-polynomial bounds have another target.
direct_ranges = cartesian(box)
models = cartesian(taylor_models(box; order = 3))
snapshots = compile.(models)
model_ranges = enclose.(snapshots)
@assert all(isguaranteed, vcat(direct_ranges, stored_ranges, model_ranges))
subbox = fill(interval(Float64, -1 // 4, 1 // 4), 6)
root_remainders = remainder.(snapshots)
restricted = [enclose(m, subbox) for m in snapshots]
fresh = compile.(cartesian(taylor_models(subbox; order = 3)))
fresh_ranges = enclose.(fresh)
@assert all(isequal_interval(a, b) && decoration(a) == decoration(b) && isguaranteed(a) == isguaranteed(b) for (a, b) in zip(remainder.(snapshots), root_remainders))
@assert all(isguaranteed, vcat(restricted, fresh_ranges))
@assert all(sup(abs(remainder(fresh[j]))) < sup(abs(root_remainders[j])) for j in 1:6)
println("Whole-box original-function state bounds: ", model_ranges)
println("Original model / freshly reevaluated remainder widths: ", (diam.(remainder.(snapshots)), diam.(remainder.(fresh))))

# ## An inclined orbit and its six-dimensional uncertainty cloud
# The cloud and error plots use independent scalar evaluations for illustration.
# A covariance ellipse is a second-moment visualization, not a confidence region
# for the bounded, nonlinear distribution. No normal-distribution assumption
# or sample count turns it into a certificate.
# Local radial/along-track axes use independent scales to make the cloud visible.
rng = Xoshiro(2026)
points = [2rand(rng, 6) .- 1 for _ in 1:1600]
truth = cartesian.(points)
linear_values, cubic_values = linear.(points), cubic.(points)
errors = [[maximum(abs, p[1:3] - q[1:3]) for (p, q) in zip(values, truth)] for values in (linear_values, cubic_values)]
@assert maximum(errors[2]) < maximum(errors[1]) / 50
nominal = cartesian(zeros(6))
radial = normalize(nominal[1:3])
normal = normalize(cross(nominal[1:3], nominal[4:6]))
along = cross(normal, radial)
frame = transpose(hcat(radial, along, normal))
offsets = [frame * (u[1:3] - nominal[1:3]) for u in truth]
frame_covariance = frame * polynomial_covariance[1:3, 1:3] * transpose(frame)
cloud_fig = Figure(size = (1300, 570), fontsize = 14)
ax = Axis3(cloud_fig[1, 1]; xlabel = "x (10³ km)", ylabel = "y (10³ km)", zlabel = "z (10³ km)", title = "Six-variable state uncertainty at one epoch", aspect = :data, azimuth = -0.8, elevation = 0.7)
mesh!(ax, Sphere(Point3f(0), 6.378f0); color = :steelblue)
orbit = [cartesian([0, 0, 0, 0, 0, (E - 1) * 25]) for E in range(0, 2π; length = 241)]
lines!(ax, first.(orbit) ./ 1000, getindex.(orbit, 2) ./ 1000, getindex.(orbit, 3) ./ 1000; color = :gray, linewidth = 2)
scatter!(ax, first.(truth) ./ 1000, getindex.(truth, 2) ./ 1000, getindex.(truth, 3) ./ 1000; color = getindex.(points, 6), colormap = :plasma, colorrange = (-1, 1), markersize = 3)
ax = Axis(cloud_fig[1, 2]; xlabel = "radial offset (km)", ylabel = "along-track offset (km)", title = "Local orbit frame: axes scaled independently")
scatter!(ax, first.(offsets), getindex.(offsets, 2); color = getindex.(points, 6), colormap = :plasma, colorrange = (-1, 1), markersize = 3)
Colorbar(cloud_fig[2, 1:2]; colormap = :plasma, limits = (-1, 1), vertical = false, label = "Normalized eccentric-anomaly uncertainty δE")
factor = eigen(Symmetric(frame_covariance[1:2, 1:2]))
ellipse = factor.vectors * Diagonal(sqrt.(factor.values)) * stack([cos(t), sin(t)] for t in range(0, 2π; length = 181))
shift = frame * (polynomial_mean[1:3] - nominal[1:3])
lines!(ax, ellipse[1, :] .+ shift[1], ellipse[2, :] .+ shift[2]; color = :black, linewidth = 2, label = "Retained polynomial: 1σ ellipse")
axislegend(ax; position = :lt, labelsize = 10)
cloud_fig

# ## Approximation accuracy and propagated coupling
# Numerical sample errors are supplementary checks. In the correlation matrix,
# different units cancel by normalization; the underlying position/velocity
# covariance uses km and km/s. Near-zero log errors use a display floor only.
accuracy_fig = Figure(size = (1180, 460), fontsize = 14)
ax = Axis(accuracy_fig[1, 1]; xlabel = "sampled position error (km)", ylabel = "fraction of points", xscale = log10, title = "Linear versus cubic uncertainty propagation")
for (label, values) in zip(("First order", "Third order"), errors)
    sorted = sort(max.(values, eps(Float64)))
    lines!(ax, sorted, collect(eachindex(sorted)) ./ length(sorted); label, linewidth = 2)
end
axislegend(ax; position = :rb, labelsize = 11)
standard_deviations = sqrt.(diag(polynomial_covariance))
correlation = polynomial_covariance ./ (standard_deviations * transpose(standard_deviations))
labels = ["x", "y", "z", "vx", "vy", "vz"]
ax = Axis(accuracy_fig[1, 2]; xticks = (1:6, labels), yticks = (1:6, labels), title = "Retained-polynomial state correlation", aspect = DataAspect())
heat = heatmap!(ax, 1:6, 1:6, correlation; colormap = :balance, colorrange = (-1, 1))
Colorbar(accuracy_fig[1, 3], heat; label = "Correlation")
accuracy_fig

# ## Bounds and genuine reevaluation on a child domain
# The first panel separates direct original-function ranges from stored-polynomial
# ranges and remainder-aware original-function ranges. The second compares
# restricting an existing model with reevaluating the original expression: only
# reevaluation recomputes a remainder for the smaller physical uncertainty set.
# Physical variation dominates both subbox ranges, so the third panel isolates
# the smaller remainder. Its bound is not the full fit error: retained interval
# coefficients have their own rounding uncertainty.
bounds_fig = Figure(size = (1530, 470), fontsize = 13)
for (column, sets, names, title) in ((1, (direct_ranges, stored_ranges, model_ranges), ("Direct intervals", "Stored polynomial", "Taylor model"), "Whole prior: different enclosure targets"), (2, (restricted, fresh_ranges), ("Restricted original model", "Fresh child models"), "Smaller prior: retain or recompute remainder"))
    ax = Axis(bounds_fig[1, column]; xticks = (1:3, ["x", "y", "z"]), ylabel = "enclosure width (km)", title)
    for (j, (values, label)) in enumerate(zip(sets, names))
        offset = (j - (length(sets) + 1) / 2) * 0.22
        barplot!(ax, (1:3) .+ offset, diam.(values[1:3]); width = 0.2, label)
    end
    axislegend(ax; position = :lt, labelsize = 10)
end
ax = Axis(bounds_fig[1, 3]; xticks = (1:3, ["x", "y", "z"]), ylabel = "bound for remainder only (km)", yscale = log10, title = "Fresh evaluation reduces the remainder")
for (j, (values, label)) in enumerate(((snapshots, "Original models"), (fresh, "Fresh child models")))
    barplot!(ax, (1:3) .+ (j - 1.5) * 0.22, [sup(abs(remainder(m))) for m in values[1:3]]; width = 0.2, fillto = 1.0e-8, label)
end
axislegend(ax; position = :lt, labelsize = 10)
bounds_fig
