# SPDX-License-Identifier: MIT #src
# Copyright (c) 2026 Jack Yarndley #src
# # Differential-algebra initial orbit determination
#
# Three right-ascension/declination pairs determine a local orbit branch.
# Differential-algebra initial orbit determination (DAIOD) expands that solution
# in the six measured angles. This tutorial follows the velocity-matching
# formulation of [Armellin2016](@citet): Gauss initialization, two Lambert arcs,
# and inversion of the velocity-defect map. All orbital helpers are included.
#
# The resulting map supports nonlinear uncertainty propagation, prediction of
# another optical observation, and a covariance-based measurement update.
# Distances are kilometres, times are seconds, and angles are radians in one
# Earth-centred inertial frame. The dynamics are two-body; observer positions
# and times are exact. Refraction, light time, and gravitational perturbations
# are omitted. The Lambert solver deliberately covers only the short-way,
# zero-revolution branch needed by this example.
using DifferentialAlgebra
using LinearAlgebra, Random
using OrdinaryDiffEqVerner
using SciMLBase: successful_retcode
using CairoMakie

set_theme!(palette = (color = Makie.to_colormap(:tab10),))

# ## Synthetic data: an inclined Kepler ellipse
#
# This analytic orbit generates observations and provides reference states;
# the IOD functions receive no truth state.
function truth(t; mu = 398600.4418)
    a, e, inclination, M0 = 12000.0, 0.15, deg2rad(50.0), 0.6
    n = sqrt(mu / a^3)
    M = M0 + n * t
    E = M
    for _ in 1:12
        E -= (E - e * sin(E) - M) / (1 - e * cos(E))
    end
    p, q = [1.0, 0.0, 0.0], [0.0, cos(inclination), sin(inclination)]
    r = a * ((cos(E) - e) * p + sqrt(1 - e^2) * sin(E) * q)
    v = a * n / (1 - e * cos(E)) *
        (-sin(E) * p + sqrt(1 - e^2) * cos(E) * q)
    return r, v
end

# The observer moves with a spherical Earth. The chosen arc is above its
# horizon. There are no refraction, light-time, precession, or J2 corrections.
function observer(t)
    latitude = deg2rad(30.0)
    longitude = deg2rad(20.0) + 7.292115e-5 * t
    return 6378.137 * [
        cos(latitude) * cos(longitude),
        cos(latitude) * sin(longitude), sin(latitude),
    ]
end

direction(alpha, delta) = [cos(delta) * cos(alpha), cos(delta) * sin(alpha), sin(delta)]

function observations(times)
    sites = hcat(observer.(times)...)
    angles = zeros(2, length(times))
    for (i, t) in enumerate(times)
        r, _ = truth(t)
        line = r - sites[:, i]
        @assert dot(line, sites[:, i]) > 0 "Synthetic target is below the horizon"
        u = line / norm(line)
        angles[:, i] = [atan(u[2], u[1]), asin(u[3])]
    end
    return angles, sites
end

# ## Solve the sparse eighth-degree radius polynomial
#
# x^8 + a*x^6 + b*x^3 + c = 0. Scale r = length_scale*x before forming
# the companion matrix: dimensional coefficients would span many powers of km.
# Retain every positive real root; a polynomial root is only a candidate.
function positive_radii(a, b, c, length_scale)
    coefficients = [
        c / length_scale^8, 0.0, 0.0, b / length_scale^5,
        0.0, 0.0, a / length_scale^2, 0.0,
    ]
    companion = zeros(8, 8)
    companion[2:8, 1:7] = Matrix{Float64}(I, 7, 7)
    companion[:, 8] = -coefficients
    roots = eigvals(companion)
    return sort(
        [
            real(x) * length_scale for x in roots
                if real(x) > 0 && abs(imag(x)) < 1.0e-9 * max(1, abs(real(x)))
        ]
    )
end

# ## Gauss: eliminate the unknown line-of-sight ranges
#
# At each epoch r_i = R_i + rho_i*u_i. The two-body f,g relations imply
# c1*r1 - r2 + c3*r3 = 0. Expanding f and g about t2 to their leading
# gravity terms gives c1 = a1 + b1*mu/r^3 and c3 = a3 + b3*mu/r^3.
# Solving U*M = R, where columns of U are u_i, gives
# rho2 = A + mu*B/r^3. Squaring r2 = R2 + rho2*u2 produces
#
# r^8 - (A^2 + 2*A*E + R2⋅R2)*r^6
#     - 2*mu*B*(A + E)*r^3 - (mu*B)^2 = 0,  E = R2⋅u2.
#
# This is a short-arc approximation. It does not enforce exact two-body
# motion; the velocity-matching correction below removes that truncation error.
function gauss(angles, sites, times, mu)
    size(angles) == (2, 3) && size(sites) == (3, 3) ||
        throw(ArgumentError("Gauss requires three RA/Dec pairs and observer positions"))
    t1, t3 = times[1] - times[2], times[3] - times[2]
    t1 < 0 < t3 || throw(ArgumentError("Epochs must increase"))
    U = hcat([direction(angles[1, i], angles[2, i]) for i in 1:3]...)
    abs(det(U)) > 1.0e-12 || throw(ArgumentError("Singular line-of-sight geometry"))
    M = U \ sites
    span = t3 - t1
    a1, a3 = t3 / span, -t1 / span
    b1 = t3 * (span^2 - t3^2) / (6span)
    b3 = -t1 * (span^2 - t1^2) / (6span)
    A = M[2, 1] * a1 - M[2, 2] + M[2, 3] * a3
    B = M[2, 1] * b1 + M[2, 3] * b3
    E = dot(sites[:, 2], U[:, 2])
    a = -(A^2 + 2A * E + sum(abs2, sites[:, 2]))
    b, c = -2mu * B * (A + E), -(mu * B)^2
    radii = positive_radii(a, b, c, norm(sites[:, 2]))
    candidates = NamedTuple[]
    for radius in radii
        gravity = mu / radius^3
        weights = [a1 + b1 * gravity, -1.0, a3 + b3 * gravity]
        minimum(abs, weights) > 1.0e-12 || continue
        rho = -(M * weights) ./ weights
        all(x -> isfinite(x) && x > 0, rho) || continue
        positions = sites + U .* reshape(rho, 1, 3)
        f1, f3 = 1 - gravity * t1^2 / 2, 1 - gravity * t3^2 / 2
        g1, g3 = t1 - gravity * t1^3 / 6, t3 - gravity * t3^3 / 6
        denominator = f1 * g3 - f3 * g1
        abs(denominator) > 1.0e-10 * span || continue
        v = (-f3 * positions[:, 1] + f1 * positions[:, 3]) / denominator
        push!(candidates, (; rho, r = positions[:, 2], v, radius))
    end
    isempty(candidates) && error("No positive-range Gauss candidate")
    return candidates
end

# ## Lambert's boundary-value problem, written out
#
# Given r_a, r_b and flight time dt, recover endpoint velocities. This teaching
# solver handles a zero-revolution, short-way branch with |z| <= 4. The finite
# bracket is deliberate: it is sufficient for these short arcs, not a general
# multi-revolution Lambert solver. Parallel/antiparallel endpoints are rejected.
#
# C(z) = sum_k (-z)^k/(2k+2)!, S(z) = sum_k (-z)^k/(2k+3)! are Stumpff
# functions. Horner evaluation avoids cancellation at z=0. Nineteen terms are
# sufficient over this bracket. The same arithmetic later works on Taylor polynomials.
scalar(x::Real) = x
vecnorm(x) = sqrt(sum(xi * xi for xi in x))

const STUMPFF_C = [Float64((-1)^k / factorial(big(2k + 2))) for k in 0:18]
const STUMPFF_S = [Float64((-1)^k / factorial(big(2k + 3))) for k in 0:18]

function stumpff(z)
    C, S, Cp, Sp = zero(z), zero(z), zero(z), zero(z)
    for k in 18:-1:0
        #= Differentiate the Horner recursion along with the polynomial. =#
        Cp, Sp = C + z * Cp, S + z * Sp
        ck = STUMPFF_C[k + 1]
        sk = STUMPFF_S[k + 1]
        C, S = ck + z * C, sk + z * S
    end
    return C, S, Cp, Sp
end

function lambert_equation(z, ra, rb, A, dt, mu)
    C, S, Cp, Sp = stumpff(z)
    y = ra + rb + A * (z * S - 1) / sqrt(C)
    scalar(y) > 0 || throw(DomainError(scalar(y), "Lambert y must be positive"))
    yp = A * ((S + z * Sp) / sqrt(C) - (z * S - 1) * Cp / (2C^(3 / 2)))
    w = y / C
    wp = (yp * C - y * Cp) / C^2
    F = w^(3 / 2) * S + A * sqrt(y) - sqrt(mu) * dt
    Fp = 1.5 * sqrt(w) * wp * S + w^(3 / 2) * Sp + A * yp / (2sqrt(y))
    return F, Fp, y
end

function lambert_root(ra, rb, A, dt, mu)
    #= Scan only to obtain an admissible sign bracket; never cross y <= 0. =#
    previous = nothing
    bracket = nothing
    for z in range(-4.0, 4.0; length = 257)
        C, S, _, _ = stumpff(z)
        y = ra + rb + A * (z * S - 1) / sqrt(C)
        y > 0 || continue
        F = first(lambert_equation(z, ra, rb, A, dt, mu))
        F == 0 && return z
        if previous !== nothing && signbit(F) != signbit(previous[2])
            bracket = (previous[1], z)
            break
        end
        previous = (z, F)
    end
    bracket === nothing && error("No Lambert root in the tutorial's |z| <= 4 branch")
    lo, hi = bracket
    flo = first(lambert_equation(lo, ra, rb, A, dt, mu))
    for _ in 1:60
        mid = (lo + hi) / 2
        fmid = first(lambert_equation(mid, ra, rb, A, dt, mu))
        if signbit(fmid) == signbit(flo)
            lo, flo = mid, fmid
        else
            hi = mid
        end
    end
    return (lo + hi) / 2
end

# The polynomial specialization below lifts this root into the angle variables.
lift_lambert_root(z, ra, rb, A, dt, mu) = z

function lambert(r_a, r_b, dt, mu)
    dt > 0 && mu > 0 || throw(ArgumentError("dt and mu must be positive"))
    ra, rb = vecnorm(r_a), vecnorm(r_b)
    cosine = sum(r_a .* r_b) / (ra * rb)
    abs(scalar(cosine)) < 1 - 1.0e-12 ||
        throw(ArgumentError("Collinear Lambert endpoints"))
    A = sqrt(ra * rb * (1 + cosine)) # Positive A selects the short way.
    z0 = lambert_root(scalar(ra), scalar(rb), scalar(A), dt, mu)
    z = lift_lambert_root(z0, ra, rb, A, dt, mu)
    F, _, y = lambert_equation(z, ra, rb, A, dt, mu)
    abs(scalar(F)) / sqrt(mu) < 1.0e-8 || error("Lambert time equation did not converge")
    f, g, gdot = 1 - y / ra, A * sqrt(y / mu), 1 - y / rb
    v_a = (r_b - f * r_a) / g
    v_b = (gdot * r_b - r_a) / g
    return v_a, v_b
end

# ## Match the incoming and outgoing velocities at the middle epoch
#
# For trial ranges rho, form r_i = R_i + rho_i*u_i, then solve Lambert on
# (t1,t2) and (t2,t3). A single orbit must satisfy the three scalar equations
# F(rho,angles) = v2_left - v2_right = 0.
# Newton solves J*delta_rho = -F. A line search keeps ranges positive and
# reduces the discontinuity. No true range or state is used as a starting guess.
function velocity_defect(rho, angles, sites, times, mu)
    positions = [
        sites[:, i] + rho[i] * direction(angles[1, i], angles[2, i])
            for i in 1:3
    ]
    _, left = lambert(positions[1], positions[2], times[2] - times[1], mu)
    right, _ = lambert(positions[2], positions[3], times[3] - times[2], mu)
    return left - right, vcat(positions[2], (left + right) / 2)
end

function match_ranges(seed, angles, sites, times, mu)
    rho = copy(seed)
    history = Float64[]
    for _ in 1:20
        F, state = velocity_defect(rho, angles, sites, times, mu)
        push!(history, norm(F))
        norm(F) < 2.0e-11 && return (; rho, state, history)
        #= Central differences are visible here for teaching. DAIOD replaces =#
        #= these numerical derivatives with polynomial coefficients. =#
        J = zeros(3, 3)
        for j in 1:3
            step = max(1.0e-3, 1.0e-5 * abs(rho[j]))
            plus, minus = copy(rho), copy(rho)
            plus[j] += step
            minus[j] -= step
            J[:, j] = (
                first(velocity_defect(plus, angles, sites, times, mu)) -
                    first(velocity_defect(minus, angles, sites, times, mu))
            ) / (2step)
        end
        cond(J) < 1.0e10 || error("Ranges are locally unobservable: singular Jacobian")
        correction = -(J \ F)
        accepted = false
        for factor in 0.5 .^ (0:16)
            trial = rho + factor * correction
            all(>(0), trial) || continue
            #= A failed trial can leave the Lambert branch; reduce its step. =#
            trial_F = try
                first(velocity_defect(trial, angles, sites, times, mu))
            catch err
                err isa DomainError || err isa ErrorException || rethrow()
                continue
            end
            if norm(trial_F) < norm(F)
                rho, accepted = trial, true
                break
            end
        end
        accepted || error("Range correction failed to reduce the velocity discontinuity")
    end
    error("Range correction exceeded its iteration limit")
end

# ## Lift the Lambert root and invert the IOD constraints
#
# Once the scalar Lambert root is known, Newton arithmetic in the polynomial
# algebra recovers its derivatives. Each iteration doubles the resolved degree.
# We check every retained coefficient of the time-of-flight equation.
scalar(x::TaylorPolynomial) = constant_term(x)
function lift_lambert_root(z0, ra::TaylorPolynomial, rb, A, dt, mu)
    z = zero(ra) + z0
    for _ in 1:(ceil(Int, log2(max_order() + 1)) + 1)
        F, Fp, _ = lambert_equation(z, ra, rb, A, dt, mu)
        z -= F / Fp
    end
    F = first(lambert_equation(z, ra, rb, A, dt, mu))
    @assert coefficient_norm(F) / sqrt(mu) < 1.0e-8
    return z
end

# Normalize angle offsets by σ and range offsets by 100 km. The first six
# coordinates are ε = (δα₁, δδ₁, δα₂, δδ₂, δα₃, δδ₃)/σ; the last three
# are auxiliary range corrections η. Let F(ε,η) be the velocity defect and
# W = ∂F/∂η at the origin. Invert the square map
#
# ```math
# H(\varepsilon,\eta) = \begin{bmatrix}
# W^{-1}(F(\varepsilon,\eta)-F(0,0)) \\ \varepsilon
# \end{bmatrix}.
# ```
#
# Evaluating H⁻¹ at (-W⁻¹F(0,0), ε) gives the range corrections as functions
# of angle offsets alone. Substitution into the Lambert solution gives the
# six-component state map at the middle observation. The nominal solve above
# makes the constant defect negligible before this local inversion.
# Nonsingular W is essential: a local inverse does not resolve ambiguous
# orbit branches or nearly unobservable observation geometry.
function daiod(nominal_rho, angles, sites, times, mu, sigma; order = 4)
    ξ = variables((:α₁, :δ₁, :α₂, :δ₂, :α₃, :δ₃, :ρ₁, :ρ₂, :ρ₃); order)
    range_scale = 100.0
    rho = nominal_rho + range_scale * ξ[7:9]
    angle_map = angles + reshape(sigma * ξ[1:6], 2, 3)
    defect, _ = velocity_defect(rho, angle_map, sites, times, mu)
    W = linear_part(defect)[:, 7:9]
    @assert cond(W) < 1.0e10
    scaled = inv(W) * defect
    offsets = constant_term.(scaled)
    inverse = invert(vcat(scaled - offsets, ξ[1:6]))
    implicit = evaluate(inverse, vcat(zero.(ξ[7:9]) - offsets, ξ[1:6]))
    ranges = nominal_rho + range_scale * implicit[7:9]
    residual, state = velocity_defect(ranges, angle_map, sites, times, mu)
    residual_norm = maximum(coefficient_norm, residual)
    @assert residual_norm < 1.0e-8
    @assert all(iszero(differentiate(p, j)) for p in vcat(ranges, state) for j in 7:9)
    @assert nvariables(CompiledMap(state)) == 6
    return (; ranges, state, condition = cond(W), residual_norm)
end

# ## Observations and the nominal orbit
#
# An inclined ellipse generates three above-horizon observations. A fixed
# noise realization makes the example reproducible. The IOD solvers receive
# only angles, observer locations, epochs, and μ, never the truth state.
# This geometry has one positive-range Gauss candidate. With ambiguous data,
# refine each admissible candidate and retain a separate map for each branch.
mu = 398600.4418
times = [-120.0, 0.0, 150.0]
sigma = deg2rad(2 / 3600) # Standard deviation of each coordinate angle: 2 arcsec.
clean_angles, sites = observations(times)
angles = clean_angles + sigma * reshape([0.6, -0.4, -0.8, 0.2, 0.3, 0.5], 2, 3)
seeds = gauss(angles, sites, times, mu)
@assert length(seeds) == 1
nominal = match_ranges(only(seeds).rho, angles, sites, times, mu)
orbit = daiod(nominal.rho, angles, sites, times, mu, sigma)
state_map = CompiledMap(orbit.state)
range_map = CompiledMap(orbit.ranges)
@assert norm(state_map(zeros(6)) - nominal.state) < 1.0e-5
exact_state = vcat(truth(0.0)...)
clean_solution = match_ranges(only(gauss(clean_angles, sites, times, mu)).rho, clean_angles, sites, times, mu)
@assert norm(clean_solution.state[1:3] - exact_state[1:3]) < 1.0e-4
@assert norm(clean_solution.state[4:6] - exact_state[4:6]) < 1.0e-7
println("Observation times (s): ", times)
println("Observed RA/Dec (degrees; columns are epochs): ", rad2deg.(angles))
println("Positive-range Gauss candidates: ", length(seeds))
println("Gauss ranges (km): ", only(seeds).rho)
println("Refined ranges (km): ", nominal.rho)
println("Velocity-defect history (km/s): ", nominal.history)
println("Range-Jacobian condition number: ", orbit.condition)
println("Largest map velocity-defect coefficient (km/s): ", orbit.residual_norm)
println("Nominal state [km; km/s]: ", nominal.state)
println("Position error from this noise realization (km): ", norm(nominal.state[1:3] - exact_state[1:3]))
# The map has six inputs even though nine variables were needed to construct
# the inverse. Its sensitivities below are with respect to physical radians.
println("Position sensitivity to the six angles (km/rad): ", linear_part(orbit.state)[1:3, 1:6] / sigma)

# ## Validate finite angle perturbations
#
# Compare degree-one, degree-two, and degree-four maps with fresh nonlinear
# range solves. Test the centre, coordinate axes, all corners of [-4,4]⁶,
# and seeded random points. The scalar reference shares the Lambert model
# but performs no polynomial arithmetic. A separate ODE check below tests
# the recovered orbit against the observations using another algorithm.
orders = [1, 2, 4]
maps = [CompiledMap(DifferentialAlgebra.trim.(orbit.state, 0, q)) for q in orders]
rng = Xoshiro(2026)
corners = stack(collect(point) for point in Iterators.product(ntuple(_ -> (-4.0, 4.0), 6)...))
axis_points = hcat(4Matrix{Float64}(I, 6, 6), -4Matrix{Float64}(I, 6, 6))
validation_points = hcat(zeros(6), axis_points, reshape(corners, 6, :), 8rand(rng, 6, 32) .- 4)
references = [match_ranges(nominal.rho, angles + reshape(sigma * point, 2, 3), sites, times, mu) for point in eachcol(validation_points)]
reference_states = stack(reference.state for reference in references)
predictions = [stack(map(point) for point in eachcol(validation_points)) for map in maps]
position_errors = [vec(sqrt.(sum(abs2, (values - reference_states)[1:3, :]; dims = 1))) for values in predictions]
velocity_errors = [vec(sqrt.(sum(abs2, (values - reference_states)[4:6, :]; dims = 1))) for values in predictions]
@assert maximum(position_errors[3]) < maximum(position_errors[2]) < maximum(position_errors[1])
@assert maximum(position_errors[3]) < 0.01
@assert maximum(velocity_errors[3]) < 1.0e-5
@assert all(all(>(0), range_map(point)) for point in eachcol(validation_points))
for (q, dr, dv) in zip(orders, position_errors, velocity_errors)
    println("Degree ", q, ": maximum position/velocity errors (km, km/s): ", (maximum(dr), maximum(dv)))
end

# ## Turn the map into an uncertainty model
#
# Assign independent standard normals to ε. Thus θ = θobs + σε is an
# uncertainty model in measurement space; the map alone is not a probability
# law. Right ascension here is a coordinate angle, not α cos(δ). Correlated
# errors can instead be introduced with θ = θobs + Lε, LLᵀ = Cθ.
#
# For independent normals E[εᵏ] is zero for odd k and (k-1)!! for even k.
# Sum coefficient lists to obtain the mean and covariance. Products of
# degree-four polynomials require moments through degree eight: ordinary
# multiplication in the order-four algebra would silently omit these terms.
# These are exact moments of the polynomial, not of the full IOD solution.
normal_moment(powers) = any(isodd, powers) ? 0.0 : prod(k -> prod(1:2:(Int(k) - 1); init = 1.0), powers)
function normal_statistics(map)
    mean = [sum(m.coefficient * normal_moment(m.exponents) for m in monomials(p); init = 0.0) for p in map]
    centered = monomials.(map - mean)
    covariance = [
        sum(
            a.coefficient * b.coefficient * normal_moment(a.exponents + b.exponents)
                for a in centered[i], b in centered[j]; init = 0.0
        )
            for i in eachindex(map), j in eachindex(map)
    ]
    return mean, covariance
end
state_mean, state_covariance = normal_statistics(orbit.state)
J = linear_part(orbit.state)[:, 1:6]
linear_covariance = J * J'
println("Nonlinear mean minus nominal state [km; km/s]: ", state_mean - nominal.state)
println("Linear position standard deviations (km): ", sqrt.(diag(linear_covariance)[1:3]))
println("Nonlinear position standard deviations (km): ", sqrt.(diag(state_covariance)[1:3]))
println("Nonlinear state covariance [km, km/s]:")
show(stdout, MIME"text/plain"(), state_covariance)
println()

# ## Propagate the map and predict another observation
#
# Taylor integration in time is independent of the degree-four expansion
# in measurement errors. Propagate that expansion to t = 900 s and project
# onto two orthogonal directions tangent to the nominal line of sight.
# These dimensionless coordinates avoid the right-ascension branch cut;
# near the centre they equal angular offsets in radians to first order.
function two_body(u, mu, t)
    r = u[1:3]
    return vcat(u[4:6], -mu * r / sum(abs2, r)^(3 // 2))
end
future_time = 900.0
flow_problem = ODEProblem(two_body, orbit.state, (0.0, future_time), mu)
flow = solve(flow_problem, TaylorMethod(18); abstol = 1.0e-12, reltol = 1.0e-12, save_everystep = false, dense = false)
@assert successful_retcode(flow) && last(flow.t) == future_time
future_map = last(flow.u)
future_site = observer(future_time)
centre_line = constant_term.(future_map[1:3]) - future_site
alpha0 = atan(centre_line[2], centre_line[1])
delta0 = asin(centre_line[3] / norm(centre_line))
east = [-sin(alpha0), cos(alpha0), 0.0]
north = [-sin(delta0) * cos(alpha0), -sin(delta0) * sin(alpha0), cos(delta0)]
tangent = permutedims(hcat(east, north))
pointing(state) = tangent * (state[1:3] - future_site) / vecnorm(state[1:3] - future_site)
pointing_map = pointing(future_map)

# Independently integrate the recovered state with Vern9. The first checks
# verify that its line of sight reaches all three input observations; the
# later checks compare map propagation with fresh numeric trajectories.
function numeric_flow(state, time)
    time == 0 && return copy(state)
    sol = solve(ODEProblem(two_body, state, (0.0, time), mu), Vern9(); abstol = 2.0e-13, reltol = 2.0e-13, save_everystep = false, dense = false)
    @assert successful_retcode(sol) && last(sol.t) == time
    return last(sol.u)
end
observation_residual = maximum(norm(normalize(numeric_flow(nominal.state, t)[1:3] - sites[:, i]) - direction(angles[1, i], angles[2, i])) for (i, t) in enumerate(times))
@assert observation_residual < 1.0e-9
future_reference = stack(numeric_flow(state, future_time) for state in eachcol(reference_states))
compiled_future = CompiledMap(future_map)
future_values = stack(compiled_future(point) for point in eachcol(validation_points))
future_position_error = maximum(norm(error) for error in eachcol((future_values - future_reference)[1:3, :]))
compiled_pointing = CompiledMap(pointing_map)
pointing_values = stack(compiled_pointing(point) for point in eachcol(validation_points))
pointing_reference = stack(pointing(state) for state in eachcol(future_reference))
arcseconds = rad2deg(1.0) * 3600
pointing_error = arcseconds * maximum(abs, pointing_values - pointing_reference)
@assert future_position_error < 0.02
@assert pointing_error < 0.001
println("Independent ODE observation residual (unit direction): ", observation_residual)
println("Future map maximum position error (km): ", future_position_error)
println("Future pointing maximum error (arcsec): ", pointing_error)

# ## A measurement update and an association statistic
#
# Treat the map-induced state distribution as a prior and add an independent
# fourth optical observation. Its noise is isotropic in the tangent plane,
# with 2 arcsec standard deviation per component. Joint polynomial moments
# provide Cₓᵧ and Cᵧᵧ without finite differences or Monte Carlo sampling.
#
# ```math
# S=C_{yy}+R,\quad K=C_{xy}S^{-1},\quad
# \mu_x^+=\mu_x+K(y_{\mathrm{obs}}-\mu_y),\quad
# P_x^+=P_x-KC_{yx}.
# ```
#
# This is a linear minimum-mean-square update using nonlinear prior moments,
# not the exact posterior of a non-Gaussian distribution. The updated state
# remains at t₂ = 0: cross-covariance transports the later information back.
# The same mean/covariance pair can initialize an EKF or weighted batch fit;
# evaluating the map supplies particles for a particle filter. Do not treat
# the original three observations as independent new data in those methods.
joint_map = vcat(orbit.state, pointing_map)
joint_mean, joint_covariance = normal_statistics(joint_map)
predicted_mean = joint_mean[7:8]
predicted_covariance = joint_covariance[7:8, 7:8]
cross_covariance = joint_covariance[1:6, 7:8]
R = sigma^2 * Matrix{Float64}(I, 2, 2)
S = Symmetric(predicted_covariance + R)
future_truth = vcat(truth(future_time)...)
@assert dot(future_truth[1:3] - future_site, future_site) > 0
measurement = pointing(future_truth) + sigma * [0.25, -0.5]
innovation = measurement - predicted_mean
K = cross_covariance / S
updated_mean = state_mean + K * innovation
# This equivalent form avoids subtracting two almost equal covariance matrices.
_, residual_covariance = normal_statistics(orbit.state - K * pointing_map)
updated_covariance = residual_covariance + K * R * K'
@assert isposdef(Symmetric(updated_covariance))
@assert tr(updated_covariance[1:3, 1:3]) < tr(state_covariance[1:3, 1:3])
@assert isapprox(updated_covariance, state_covariance - K * cross_covariance'; rtol = 1.0e-9)
mahalanobis_squared = dot(innovation, S \ innovation)
println("Predicted pointing standard deviations (arcsec): ", arcseconds * sqrt.(diag(predicted_covariance)))
println("Fourth-observation innovation (arcsec): ", arcseconds * innovation)
println("Squared Mahalanobis association statistic: ", mahalanobis_squared)
println("Updated middle-epoch state [km; km/s]: ", updated_mean)
println("Position standard deviations before/after update (km): ", (sqrt.(diag(state_covariance)[1:3]), sqrt.(diag(updated_covariance)[1:3])))
println("Position errors before/after update (km): ", (norm(state_mean[1:3] - exact_state[1:3]), norm(updated_mean[1:3] - exact_state[1:3])))
# A small Mahalanobis statistic indicates consistency with the predicted
# observation cloud. A χ² gate would require an approximately Gaussian
# innovation; this example does not claim calibrated association probabilities.

# Check all joint moments with tensor Gauss–Hermite quadrature of the map.
# Five nodes per coordinate integrate degrees through nine, covering every
# product of two quartic maps. This checks the coefficient-moment calculation,
# independently of the scalar IOD and ODE checks above.
quadrature = eigen(SymTridiagonal(zeros(5), sqrt.(Float64.(1:4))))
nodes, weights = quadrature.values, quadrature.vectors[1, :] .^ 2
indices = collect(Iterators.product(ntuple(_ -> 1:5, 6)...))
joint_weights = vec([prod(weights[i] for i in index) for index in indices])
compiled_joint = CompiledMap(joint_map)
quadrature_values = stack(compiled_joint(nodes[collect(index)]) for index in vec(indices))
quadrature_mean = quadrature_values * joint_weights
centered_values = quadrature_values .- quadrature_mean
quadrature_covariance = (centered_values .* joint_weights') * centered_values'
deviations = sqrt.(diag(joint_covariance))
mean_check = maximum(abs, (quadrature_mean - joint_mean) ./ deviations)
covariance_check = maximum(abs, (quadrature_covariance - joint_covariance) ./ (deviations * deviations'))
@assert mean_check < 1.0e-8 && covariance_check < 1.0e-8
println("Quadrature checks, standardized mean/covariance errors: ", (mean_check, covariance_check))

# ## Visualize the orbit set and the follow-up prediction
#
# The cloud uses 4,000 Gaussian draws and cheap compiled-map evaluations.
# Gaussian tails are unbounded: sampled tests on [-4,4]⁶ cannot certify the
# whole distribution. Increase order, shorten the domain, or split it and
# revalidate when angular errors grow; keep separate maps for separate roots.
# Ellipses below are covariance contours with squared Mahalanobis radius
# 5.991 (95% for a bivariate Gaussian), not guaranteed 95% nonlinear regions.
samples = randn(rng, 6, 4000)
cloud = stack(compiled_joint(point) for point in eachcol(samples))
radial = normalize(nominal.state[1:3])
normal = normalize(cross(nominal.state[1:3], nominal.state[4:6]))
along = cross(normal, radial)
projection = permutedims(hcat(radial, along))
position_cloud = projection * (cloud[1:3, :] .- nominal.state[1:3])
function ellipse(mean, covariance)
    E = eigen(Symmetric(covariance))
    @assert minimum(E.values) > 0
    phi = range(0, 2π; length = 201)
    return mean .+ sqrt(5.991) * E.vectors * Diagonal(sqrt.(E.values)) * permutedims(hcat(cos.(phi), sin.(phi)))
end
println("Gaussian cloud samples outside the validation box: ", count(point -> maximum(abs, point) > 4, eachcol(samples)), " / ", size(samples, 2))

#-
fig = Figure(size = (1120, 800), fontsize = 15)
colors = Makie.to_colormap(:tab10)
convergence_axis = Axis(fig[1, 1]; xlabel = "Range correction iteration", ylabel = "Velocity defect (km/s)", yscale = log10, title = "Gauss seed → matched Lambert arcs")
scatterlines!(convergence_axis, eachindex(nominal.history), nominal.history; color = colors[1], linewidth = 2)
error_axis = Axis(fig[1, 2]; xlabel = "Position error against scalar IOD (km)", ylabel = "Fraction of validation points", xscale = log10, title = "Finite-perturbation map validation")
for (i, q) in enumerate(orders)
    sorted = sort(position_errors[i])
    stairs!(error_axis, max.(sorted, eps(Float64)), (1:length(sorted)) ./ length(sorted); color = colors[i], label = "Degree $q", linewidth = 2)
end
axislegend(error_axis; position = :lt)
orbit_axis = Axis(fig[2, 1]; xlabel = "Radial offset from nominal (km)", ylabel = "Along-track offset from nominal (km)", title = "Middle-epoch uncertainty and update")
scatter!(orbit_axis, position_cloud[1, :], position_cloud[2, :]; color = (colors[1], 0.15), markersize = 3)
for (mean, covariance, label, color) in ((state_mean, state_covariance, "Prior covariance", colors[1]), (updated_mean, updated_covariance, "Updated covariance", colors[2]))
    contour = ellipse(projection * (mean[1:3] - nominal.state[1:3]), projection * covariance[1:3, 1:3] * projection')
    lines!(orbit_axis, contour[1, :], contour[2, :]; color, label, linewidth = 2)
end
truth_offset = projection * (exact_state[1:3] - nominal.state[1:3])
scatter!(orbit_axis, [truth_offset[1]], [truth_offset[2]]; color = :black, marker = :star5, markersize = 12, label = "Truth")
axislegend(orbit_axis; position = :lt, labelsize = 12)
pointing_axis = Axis(fig[2, 2]; xlabel = "East tangent offset (arcsec)", ylabel = "North tangent offset (arcsec)", title = "Follow-up observation at t = 900 s")
scatter!(pointing_axis, arcseconds * cloud[7, :], arcseconds * cloud[8, :]; color = (colors[1], 0.15), markersize = 3)
contour = arcseconds * ellipse(predicted_mean, S)
lines!(pointing_axis, contour[1, :], contour[2, :]; color = colors[1], linewidth = 2, label = "Prediction + measurement noise")
scatter!(pointing_axis, [arcseconds * measurement[1]], [arcseconds * measurement[2]]; color = colors[2], marker = :cross, markersize = 14, label = "Fourth observation")
axislegend(pointing_axis; position = :lt, labelsize = 11)
fig
