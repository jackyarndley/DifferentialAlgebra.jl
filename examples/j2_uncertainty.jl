# # Six-variable injection uncertainty under J2 gravity
# Follow an inclined low-Earth orbit for three revolutions, including the leading
# oblateness term. Propagate six uncertain Cartesian injection coordinates with
# a cubic DA map, then inspect the orbital tube, Earth-fixed ground tracks and
# the growth of radial/along-track/cross-track uncertainty.
using DifferentialAlgebra, OrdinaryDiffEqVerner, CairoMakie, LinearAlgebra, Random, ForwardDiff
using SciMLBase: successful_retcode

# ## Force model and independent expansion orders
# Length is in Earth equatorial radii and μ=1. The acceleration is -∇U for
# U=-1/r + J2*(3z²/r²-1)/(2r³). This is the Cartesian J2 model described in
# [NASA's equations-of-motion check cases](https://nescacademy.nasa.gov/src/flightsim/Reports/NASA-TM-2015-218675-EOM_checkcase_appendices.pdf).
# Drag, higher harmonics, third bodies and maneuver errors are omitted.
function j2_dynamics(u, J2, t)
    x, y, z, vx, vy, vz = u
    r2 = x^2 + y^2 + z^2
    factor = -1 / (r2 * sqrt(r2))
    j = 3J2 / (2r2)
    q = 5z^2 / r2
    return [vx, vy, vz, factor * x * (1 + j * (1 - q)), factor * y * (1 + j * (1 - q)), factor * z * (1 + j * (3 - q))]
end
function potential(r, J2)
    r2 = sum(x -> x^2, r)
    return -1 / sqrt(r2) + J2 * (3r[3]^2 / r2 - 1) / (2r2 * sqrt(r2))
end
energy(u, J2) = sum(abs2, u[4:6]) / 2 + potential(u[1:3], J2)
J2 = 1.08262668e-3
radius, inclination = 1.08, 0.9
nominal = [radius, 0, 0, 0, cos(inclination) / sqrt(radius), sin(inclination) / sqrt(radius)]
widths = [fill(0.001, 3); fill(0.00005, 3)]
initial(ξ) = nominal + widths .* ξ
@assert norm(j2_dynamics(nominal, J2, 0)[4:6] + ForwardDiff.gradient(r -> potential(r, J2), nominal[1:3])) < 1.0e-14
period = 2π * radius^(3 / 2)
times = range(0, 3period; length = 81)
ξ = variables((:δx, :δy, :δz, :δvx, :δvy, :δvz); order = 3)
problem = ODEProblem(j2_dynamics, initial(ξ), (0.0, last(times)), J2)
flow = solve(problem, TaylorMethod(20); abstol = 1.0e-12, reltol = 1.0e-12, saveat = times, dense = false)
@assert successful_retcode(flow) && last(flow.t) == last(times)
@assert max_order() == 3
maps = CompiledMap.(flow.u)
linear_maps = [CompiledMap(DifferentialAlgebra.trim.(u, 0, 1)) for u in flow.u]
nominal_states = [constant_term.(u) for u in flow.u]
numeric = solve(remake(problem; u0 = nominal), Vern9(); abstol = 2.0e-14, reltol = 2.0e-14, saveat = times, dense = false)
two_body = solve(remake(problem; u0 = nominal, p = 0.0), Vern9(); abstol = 2.0e-14, reltol = 2.0e-14, saveat = times, dense = false)
@assert successful_retcode(numeric) && successful_retcode(two_body)
@assert maximum(norm(a - b) for (a, b) in zip(nominal_states, numeric.u)) < 1.0e-8
@assert maximum(abs(energy(u, J2) - energy(nominal, J2)) for u in nominal_states) < 1.0e-9
println("Time order 20; uncertainty order 3; accepted steps: ", flow.stats.naccept)

# This is numerical time integration, not a validated flow. Taylor-model
# remainders for a static map do not certify integration error. Widths are
# independent uniform half-widths (about 6.4 km and 0.40 m/s), not standard
# deviations. The terrestrial units below are illustrative stored Float64s.
earth_radius = 6378.0
time_unit = sqrt(earth_radius^3 / 398600.4418)
earth_rotation = 7.292115e-5 * time_unit
rng = Xoshiro(2026)
samples = [2rand(rng, 6) .- 1 for _ in 1:600]
validation = [[fill(s, 6) for s in (-1.0, 0.0, 1.0)]; [2rand(rng, 6) .- 1 for _ in 1:13]]
function numeric_endpoint(ξ)
    sol = solve(remake(problem; u0 = initial(ξ)), Vern9(); abstol = 2.0e-14, reltol = 2.0e-14, save_everystep = false, dense = false)
    @assert successful_retcode(sol)
    return last(sol.u)
end
references = numeric_endpoint.(validation)
errors = [[norm(m(p)[1:3] - q[1:3]) * earth_radius for (p, q) in zip(validation, references)] for m in (last(linear_maps), last(maps))]
@assert maximum(errors[2]) < maximum(errors[1]) / 20
println("Maximum independent endpoint position discrepancies (km), linear/cubic: ", maximum.(errors))

# ## The orbital tube and the local encounter geometry
# The displayed tube comprises clouds at eleven saved epochs; it is not a
# uniform enclosure between epochs. The cubic maps evaluate all cloud members
# without additional ODE solves. The final frame uses the nominal r and r×v.
function orbit_frame(u)
    radial = normalize(u[1:3])
    normal = normalize(cross(u[1:3], u[4:6]))
    return transpose(hcat(radial, cross(normal, radial), normal))
end
epochs = 1:8:length(times)
clouds = [stack(maps[k](p)[1:3] for p in samples) for k in epochs]
local_clouds = [orbit_frame(nominal_states[k]) * (cloud .- nominal_states[k][1:3]) * earth_radius for (k, cloud) in zip(epochs, clouds)]
orbit_fig = Figure(size = (1280, 520), fontsize = 14)
ax = Axis3(orbit_fig[1, 1]; xlabel = "x (10³ km)", ylabel = "y (10³ km)", zlabel = "z (10³ km)", title = "Inclined LEO: three revolutions under J2", aspect = :data, azimuth = -0.8, elevation = 0.5)
mesh!(ax, Sphere(Point3f(0), 6.378f0); color = :steelblue)
history = stack(nominal_states) * (earth_radius / 1000)
lines!(ax, history[1, :], history[2, :], history[3, :]; color = :black, linewidth = 1.5)
for (k, cloud) in zip(epochs, clouds)
    scatter!(ax, cloud[1, :] * (earth_radius / 1000), cloud[2, :] * (earth_radius / 1000), cloud[3, :] * (earth_radius / 1000); color = fill(times[k] / period, length(samples)), colormap = :plasma, colorrange = (0, 3), markersize = 3)
end
Colorbar(orbit_fig[1, 2]; colormap = :plasma, limits = (0, 3), label = "Time / nominal period")
ax = Axis(orbit_fig[1, 3]; xlabel = "radial displacement (km)", ylabel = "along-track displacement (km)", title = "Final cubic cloud in the local orbit frame")
final_cloud = scatter!(ax, last(local_clouds)[1, :], last(local_clouds)[2, :]; color = last(local_clouds)[3, :], colormap = :balance, markersize = 4)
Colorbar(orbit_fig[1, 4], final_cloud; label = "Cross-track displacement (km)")
validation_offsets = orbit_frame(last(nominal_states)) * (stack(q[1:3] for q in references) .- last(nominal_states)[1:3]) * earth_radius
scatter!(ax, validation_offsets[1, :], validation_offsets[2, :]; color = :black, marker = :cross, markersize = 10, label = "Independent trajectories")
axislegend(ax; position = :lt, labelsize = 10)
orbit_fig

# ## Rotating-Earth ground tracks
# Rotate saved inertial positions by a constant Earth rotation rate. Latitude
# is geocentric, not geodetic. Split polylines at the longitude seam; these
# sampled tracks do not model station visibility or atmospheric refraction.
function ground_track(states)
    longitude = [mod(atan(u[2], u[1]) - earth_rotation * t + π, 2π) - π for (u, t) in zip(states, times)]
    latitude = [atan(u[3], hypot(u[1], u[2])) for u in states]
    for k in reverse(2:length(longitude))
        if abs(longitude[k] - longitude[k - 1]) > π
            insert!(longitude, k, NaN)
            insert!(latitude, k, NaN)
        end
    end
    return rad2deg.(longitude), rad2deg.(latitude)
end
ground_fig = Figure(size = (1220, 430), fontsize = 14)
ax = Axis(ground_fig[1, 1]; xlabel = "Earth-fixed longitude (deg)", ylabel = "geocentric latitude (deg)", title = "Nominal and uncertain ground tracks", limits = (-180, 180, -90, 90), xticks = -180:60:180, yticks = -90:30:90)
for p in samples[1:8]
    lon, lat = ground_track([m(p) for m in maps])
    lines!(ax, lon, lat; color = (:dodgerblue, 0.25), linewidth = 1)
end
lon, lat = ground_track(nominal_states)
lines!(ax, lon, lat; color = :black, linewidth = 2, label = "Nominal J2 orbit")
axislegend(ax; position = :rt, labelsize = 10)
node(u) = atan(cross(u[1:3], u[4:6])[1], -cross(u[1:3], u[4:6])[2])
node_difference = rad2deg.(node.(numeric.u) - node.(two_body.u))
ax = Axis(ground_fig[1, 2]; xlabel = "Time / nominal period", ylabel = "node shift relative to two-body (deg)", title = "Oblateness changes the orbital plane")
lines!(ax, times ./ period, node_difference; color = :darkorange, linewidth = 2)
ground_fig

# ## Dispersion and a separate numerical accuracy check
# Standard deviations here come from the finite uniform cloud, not exact
# moments or interval bounds. Cubic versus linear differences and independent
# endpoint discrepancies help diagnose truncation, without proving inclusion.
diagnostic_fig = Figure(size = (1220, 440), fontsize = 14)
ax = Axis(diagnostic_fig[1, 1]; xlabel = "Time / nominal period", ylabel = "sample standard deviation (km)", title = "Six-coordinate uncertainty grows into along-track spread")
for (j, name) in enumerate(("Radial", "Along track", "Cross track"))
    spread = [sqrt(sum(abs2, c[j, :] .- sum(c[j, :]) / size(c, 2)) / (size(c, 2) - 1)) for c in local_clouds]
    lines!(ax, times[epochs] ./ period, spread; label = name, linewidth = 2)
end
axislegend(ax; position = :lt, labelsize = 10)
ax = Axis(diagnostic_fig[1, 2]; xlabel = "endpoint position discrepancy (km)", ylabel = "fraction of validation points", xscale = log10, title = "Sixteen independent numeric integrations")
for (label, values) in zip(("Linear map", "Cubic map"), errors)
    sorted = sort(max.(values, eps(Float64)))
    stairs!(ax, sorted, collect(eachindex(sorted)) ./ length(sorted); label, linewidth = 2)
end
axislegend(ax; position = :rb, labelsize = 10)
diagnostic_fig
