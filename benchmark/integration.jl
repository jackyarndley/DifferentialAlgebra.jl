"""
Same fixed-step RK4 propagation of a six-variable two-body Taylor map in both
libraries. Run `julia benchmark/setup.jl`, then
`julia --project=benchmark benchmark/integration.jl` from the repository root.
Initialization/compilation are outside timing. GC remains enabled. Every stored
coefficient is checked, degree by degree after removing the perturbation scale.
"""
module IntegrationComparison
import DifferentialAlgebra
import TaylorSeries as TS
using LinearAlgebra, Statistics, Printf
using TOML
using OrdinaryDiffEqVerner: OrdinaryDiffEqVerner, ODEProblem, solve, Vern9

function initial(::Type{T}, order, backend) where {T}
    variables = if backend === :DifferentialAlgebra
        DifferentialAlgebra.init(order, 6)
        [DifferentialAlgebra.variable(i, T) for i in 1:6]
    else
        # OrdinaryDiffEq sometimes calls zero(eltype(u)), which needs a default
        # space. Explicit independent JetSpaces work with our standalone RK4.
        TS.variables!(T, [Symbol("x", i) for i in 1:6]; order, nowarn = true)
    end
    scale = one(T) / 1000
    nominal = T[1, 0, 0, 0, 1, 0]
    return [nominal[i] + scale * variables[i] for i in 1:6]
end

rhs!(du, u, p, t) = rhs!(du, u)
constant(p::DifferentialAlgebra.DA) = DifferentialAlgebra.cons(p)
constant(p::TS.TaylorN) = TS.constant_term(p)
constant(x::Real) = x
constant_norm(u, t) = sqrt(sum(x -> abs2(constant(x)), u))
function vern9(u0, dt, steps)
    problem = ODEProblem(rhs!, u0, (zero(dt), steps * dt))
    solution = solve(
        problem, Vern9(); adaptive = false, dt, save_everystep = false,
        save_start = false, dense = false, internalnorm = constant_norm
    )
    return solution.u[end]
end

function rhs!(du, u)
    r2 = u[1]^2 + u[2]^2 + u[3]^2
    gravity = -inv(r2 * sqrt(r2))
    @inbounds for i in 1:3
        du[i] = u[i + 3]
        du[i + 3] = gravity * u[i]
    end
    return du
end

function rk4(u0, dt, steps)
    u = copy(u0)
    k1, k2, k3, k4, tmp = (similar(u) for _ in 1:5)
    @inbounds for _ in 1:steps
        rhs!(k1, u)
        for i in 1:6
            tmp[i] = u[i] + (dt / 2) * k1[i]
        end
        rhs!(k2, tmp)
        for i in 1:6
            tmp[i] = u[i] + (dt / 2) * k2[i]
        end
        rhs!(k3, tmp)
        for i in 1:6
            tmp[i] = u[i] + dt * k3[i]
        end
        rhs!(k4, tmp)
        for i in 1:6
            u[i] = u[i] + (dt / 6) * (k1[i] + 2k2[i] + 2k3[i] + k4[i])
        end
    end
    return u
end

function measure(f, g; samples = 5)
    f(); g(); f(); g() # compile and warm lookup tables in both libraries
    times, bytes = (Float64[], Float64[]), (Int[], Int[])
    for sample in 1:samples
        for backend in (isodd(sample) ? (1, 2) : (2, 1))
            GC.gc()
            result = backend == 1 ? (@timed f()) : (@timed g())
            push!(times[backend], result.time); push!(bytes[backend], result.bytes)
        end
    end
    return times, Int.(median.(bytes))
end

function compare(::Type{T}, order; steps = 256, samples = 5, integrator = rk4) where {T}
    native = initial(T, order, :DifferentialAlgebra)
    taylor = initial(T, order, :TaylorSeries)
    dt = 2convert(T, π) / steps
    fd = () -> integrator(native, dt, steps)
    ft = () -> integrator(taylor, dt, steps)
    a, b = fd(), ft()
    @assert all(p -> constant(p) isa T, a) && all(p -> constant(p) isa T, b)
    errors = zeros(T, order + 1)
    magnitudes = zeros(T, order + 1)
    for powers in DifferentialAlgebra.getMultiIndices(order, 6)
        degree = sum(powers)
        scale = (one(T) / 1000)^degree
        for i in 1:6
            ca, cb = DifferentialAlgebra.getCoefficient(a[i], powers) / scale, TS.getcoeff(b[i], Int.(powers)) / scale
            errors[degree + 1] = max(errors[degree + 1], abs(ca - cb))
            magnitudes[degree + 1] = max(magnitudes[degree + 1], abs(ca), abs(cb))
        end
    end
    error = maximum(errors ./ max.(one(T), magnitudes))
    @assert error < (T === Float32 ? T(0.005) : T === BigFloat ? big"1e-60" : T(2.0e-10)) error
    nominal = T[1, 0, 0, 0, 1, 0]
    @assert maximum(abs.(DifferentialAlgebra.cons.(a) - nominal)) < T(2.0e-5)
    if integrator === rk4 && T === Float64
        coarse = maximum(abs.(rk4(nominal, dt, steps) - nominal))
        fine = maximum(abs.(rk4(nominal, dt / 2, 2steps) - nominal))
        @assert 10 < coarse / fine < 20 # fourth-order global integration error
    end
    # Independent scalar propagation at a nonzero perturbation, including all
    # six state components. The map's omitted terms shrink with order.
    delta = T[0.02, -0.015, 0.01, -0.02, 0.01, 0.015]
    scalar = integrator(nominal + delta / 1000, dt, steps)
    @assert maximum(abs.(DifferentialAlgebra.evaluate.(a, Ref(delta)) - scalar)) < T(2.0e-5)
    times, (bd, bt) = measure(fd, ft; samples)
    td, tt = median.(times)
    @printf(
        "%s %s order=%d steps=%d DifferentialAlgebra=%.6fs TaylorSeries=%.6fs speedup=%.2fx bytes=%d/%d coefficient_error=%.3g\n",
        integrator, T, order, steps, td, tt, tt / td, bd, bt, error
    )
    flush(stdout)
    return Dict(
        "type" => string(T), "integrator" => string(integrator), "order" => order, "steps" => steps,
        "native_seconds" => td, "taylor_seconds" => tt, "native_bytes" => bd, "taylor_bytes" => bt,
        "native_samples" => times[1], "taylor_samples" => times[2], "error" => Float64(error)
    )
end

function main()
    println("Julia ", VERSION, "; DifferentialAlgebra ", pkgversion(DifferentialAlgebra), "; TaylorSeries ", pkgversion(TS), "; OrdinaryDiffEqVerner ", pkgversion(OrdinaryDiffEqVerner), "; ", Sys.KERNEL, "/", Sys.ARCH, "; CPU ", Sys.CPU_NAME)
    quick = "--quick" in ARGS
    rows = Dict[]
    for order in (2, 3, 5)
        push!(rows, compare(Float64, order; samples = quick ? 3 : 7))
    end
    push!(rows, compare(Float32, 3; samples = quick ? 3 : 7))
    setprecision(256) do
        push!(rows, compare(BigFloat, 3; samples = quick ? 3 : 7))
    end
    push!(rows, compare(Float64, 3; steps = 64, samples = quick ? 3 : 7, integrator = vern9))
    return if !isempty(ARGS) && endswith(last(ARGS), ".toml")
        open(last(ARGS), "w") do io
            TOML.print(
                io, Dict(
                    "julia" => string(VERSION), "differentialalgebra" => string(pkgversion(DifferentialAlgebra)), "taylorseries" => string(pkgversion(TS)),
                    "ordinarydiffeqverner" => string(pkgversion(OrdinaryDiffEqVerner)), "bigfloat_precision" => 256,
                    "kernel" => string(Sys.KERNEL), "cpu" => Sys.CPU_NAME, "samples" => (quick ? 3 : 7), "cases" => rows
                )
            )
        end
    end
end
abspath(PROGRAM_FILE) == (@__FILE__) && main()
end
