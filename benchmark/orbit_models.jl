# Shared scalar-generic models for the orbit benchmarks. Nondimensional units:
# gravitational parameter μ = 1 and the Kepler semi-major axis a = 1.
primal(x::Real) = x
primal(x::DA.TaylorPolynomial) = DA.constant_term(x)
primal(x::TS.TaylorN) = TS.constant_term(x)
primal(x::ForwardDiff.Dual) = primal(ForwardDiff.value(x))
primal_norm(x::Number, t) = abs(primal(x))
primal_norm(x::AbstractArray, t) = sqrt(sum(v -> abs2(primal(v)), x) / length(x))

function orbit_rhs!(du, u, j2, t)
    x, y, z, vx, vy, vz = u
    r2 = x * x + y * y + z * z
    gravity = -inv(r2 * sqrt(r2))
    du[1], du[2], du[3] = vx, vy, vz
    if iszero(j2)
        du[4], du[5], du[6] = gravity * x, gravity * y, gravity * z
    else
        # The parameter is J₂ R²; the spin axis is the z axis.
        c = -1.5j2 * gravity / r2
        s = 5z * z / r2
        du[4] = gravity * x + c * x * (s - 1)
        du[5] = gravity * y + c * y * (s - 1)
        du[6] = gravity * z + c * z * (s - 3)
    end
    return du
end
orbit_rhs(u, j2, t) = orbit_rhs!(similar(u), u, j2, t)

function propagate(u0, tf, j2, algorithm, tolerance)
    problem = ODEProblem(orbit_rhs!, u0, (0.0, tf), j2)
    solution = solve(
        problem, algorithm; adaptive = true, dt = 0.01,
        abstol = tolerance, reltol = tolerance, internalnorm = primal_norm,
        save_start = false, save_everystep = false, dense = false, maxiters = 1_000_000
    )
    SciMLBase.successful_retcode(solution) || error("Orbit solve failed: $(solution.retcode)")
    return solution
end

# Elliptic Lagrange f/g propagation, used independently of the ODE solver.
# Newton first converges the nominal anomaly, then lifts all Taylor orders.
function kepler(state, tf, order = 0)
    r0, v0 = state[1:3], state[4:6]
    radius0 = sqrt(sum(abs2, r0))
    a = inv(2 / radius0 - sum(abs2, v0))
    primal(a) > 0 || error("This benchmark requires an elliptic orbit")
    scale = sqrt(a)
    sigma = dot(r0, v0)
    m, b, c = tf / (a * scale), sigma / scale, 1 - radius0 / a
    m0, b0, c0 = primal(m), primal(b), primal(c)
    anomaly = m0
    converged = false
    for _ in 1:40
        s, co = sincos(anomaly)
        correction = (anomaly + b0 * (1 - co) - c0 * s - m0) / (1 + b0 * s - c0 * co)
        anomaly -= correction
        if abs(correction) <= 8eps(typeof(m0)) * max(one(m0), abs(anomaly))
            converged = true
            break
        end
    end
    converged || error("Kepler iteration did not converge")
    anomaly = zero(m) + anomaly
    for _ in 1:ceil(Int, log2(order + 1))
        s, co = sincos(anomaly)
        anomaly -= (anomaly + b * (1 - co) - c * s - m) / (1 + b * s - c * co)
    end
    s, co = sincos(anomaly)
    f = 1 - a / radius0 * (1 - co)
    g = a * sigma * (1 - co) + radius0 * scale * s
    position = f * r0 + g * v0
    radius = sqrt(sum(abs2, position))
    fdot = -scale * s / (radius * radius0)
    gdot = 1 - a / radius * (1 - co)
    return vcat(position, fdot * r0 + gdot * v0)
end

# Independently integrate the 42-dimensional state/variational system. Its
# error controller includes the STM, unlike the benchmark's nominal-only norm.
function variational_reference(u0, tf, j2)
    function rhs!(du, u, p, t)
        state = @view u[1:6]
        du[1:6] .= orbit_rhs(state, p, t)
        A = ForwardDiff.jacobian(x -> orbit_rhs(x, p, t), state)
        mul!(reshape(@view(du[7:42]), 6, 6), A, reshape(@view(u[7:42]), 6, 6))
        return nothing
    end
    initial = vcat(u0, vec(Matrix{Float64}(I, 6, 6)))
    problem = ODEProblem(rhs!, initial, (0.0, tf), j2)
    solution = solve(
        problem, Vern9(); abstol = 2.0e-13, reltol = 2.0e-13,
        save_start = false, save_everystep = false, dense = false
    )
    SciMLBase.successful_retcode(solution) || error("Variational reference failed")
    result = solution.u[end]
    return result[1:6], reshape(result[7:42], 6, 6)
end

function scenarios()
    circular = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
    e, inclination = 0.4, 0.5
    speed = sqrt((1 + e) / (1 - e))
    eccentric = [1 - e, 0.0, 0.0, 0.0, speed * cos(inclination), speed * sin(inclination)]
    return [
        (name = "circular", u0 = circular, tf = 2π, j2 = 0.0, algorithm = Tsit5(), tolerance = 1.0e-9),
        (name = "circular", u0 = circular, tf = 2π, j2 = 0.0, algorithm = Vern9(), tolerance = 1.0e-11),
        (name = "eccentric_3periods", u0 = eccentric, tf = 6π, j2 = 0.0, algorithm = Vern9(), tolerance = 1.0e-11),
        (name = "j2_inclined", u0 = eccentric, tf = 2π, j2 = 1.08263e-3 * 0.5^2, algorithm = Vern9(), tolerance = 1.0e-11),
    ]
end
