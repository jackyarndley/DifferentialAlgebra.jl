"""
Compare dynamic and static arrays without changing the polynomial implementation.
Run `julia --project=benchmark benchmark/staticarrays.jl [output.toml]`.
The coefficient experiment exercises the same multiplication kernel with both
storage layouts; it is a prototype, not a second polynomial implementation.
"""
module StaticArrayComparison
using DifferentialAlgebra, StaticArrays, BenchmarkTools, Statistics, TOML

function measurement(f, evaluations)
    f(); f()
    trial = @benchmark $f() samples = 30 evals = evaluations seconds = 10
    estimate = median(trial)
    return (seconds = estimate.time / 1.0e9, bytes = estimate.memory, allocations = estimate.allocs)
end

function compare(label, dynamic, static, check; evaluations = 1)
    check(dynamic(), static()) || error("Results differ: $label")
    # Reverse the order on a second batch to expose timing noise and cache effects.
    d1, s1 = measurement(dynamic, evaluations), measurement(static, evaluations)
    s2, d2 = measurement(static, evaluations), measurement(dynamic, evaluations)
    d, s = median([d1.seconds, d2.seconds]), median([s1.seconds, s2.seconds])
    println(
        label, ": dynamic/static = ", round(d / s; digits = 3),
        ", bytes = ", d1.bytes, "/", s1.bytes
    )
    return Dict(
        "case" => label, "dynamic_seconds" => d, "static_seconds" => s,
        "speedup" => d / s, "dynamic_bytes" => d1.bytes, "static_bytes" => s1.bytes,
        "dynamic_allocations" => d1.allocations, "static_allocations" => s1.allocations
    )
end

function convolution!(output, a, b, basis, order)
    fill!(output, zero(eltype(output)))
    return DifferentialAlgebra.convolve!(output, a, b, basis, order, Val(true), Val(false))
end

function coefficient_case(order)
    x = variables(6; order)
    basis = x[1].algebra.basis
    n = DifferentialAlgebra.nmonomials()
    ac = [sin(i) / 10 for i in 1:n]
    bc = [cos(i) / 10 for i in 1:n]
    a, b = (coeffs = ac, len = n), (coeffs = bc, len = n)
    sa, sb = (coeffs = SVector{n}(ac), len = n), (coeffs = SVector{n}(bc), len = n)
    out, sout = zeros(n), MVector{n, Float64}(undef)
    return compare(
        "coefficient kernel: order $order, $n coefficients",
        () -> convolution!(out, a, b, basis, order),
        () -> convolution!(sout, sa, sb, basis, order), (a, b) -> a ≈ b; evaluations = 1000
    )
end

# Include conversion and result ownership costs when considering a static kernel
# inside a library whose coefficients have runtime-selected lengths.
function static_product(ac, bc, basis, order, ::Val{N}) where {N}
    a, b = (coeffs = SVector{N}(ac), len = N), (coeffs = SVector{N}(bc), len = N)
    out = MVector{N, Float64}(undef)
    return Vector(convolution!(out, a, b, basis, order))
end
function dynamic_product(ac, bc, basis, order)
    n = length(ac)
    return convolution!(Vector{Float64}(undef, n), (coeffs = ac, len = n), (coeffs = bc, len = n), basis, order)
end
function storage_case(order)
    x = variables(6; order)
    basis = x[1].algebra.basis
    n = DifferentialAlgebra.nmonomials()
    a, b = [sin(i) / 10 for i in 1:n], [cos(i) / 10 for i in 1:n]
    return compare(
        "coefficient storage with conversion: order $order, $n coefficients",
        () -> dynamic_product(a, b, basis, order), () -> static_product(a, b, basis, order, Val(n)),
        (a, b) -> a ≈ b; evaluations = 100
    )
end

function evaluation_case(n, order)
    x = variables(n; order)
    f = CompiledMap([sin(1 + x[i]) + exp(sum(x) / 10) for i in 1:n])
    point = fill(0.01, n)
    out, work = zeros(n), zeros(order + 1)
    spoint, sout, swork = SVector{n}(point), MVector{n, Float64}(undef), MVector{order + 1, Float64}(undef)
    return compare(
        "compiled evaluation: $n variables, order $order",
        () -> evaluate!(out, f, point, work), () -> evaluate!(sout, f, spoint, swork), (a, b) -> a ≈ b; evaluations = 1000
    )
end

function rhs_tuple(u)
    r2 = u[1]^2 + u[2]^2 + u[3]^2
    g = -inv(r2 * sqrt(r2))
    return (u[4], u[5], u[6], g * u[1], g * u[2], g * u[3])
end
rhs(u::Vector) = collect(rhs_tuple(u))
rhs(u::SVector) = SVector(rhs_tuple(u))
function rk4(u, h, steps)
    for _ in 1:steps
        k1 = rhs(u)
        k2 = rhs(u + (h / 2) * k1)
        k3 = rhs(u + (h / 2) * k2)
        k4 = rhs(u + h * k3)
        u = u + (h / 6) * (k1 + 2k2 + 2k3 + k4)
    end
    return u
end
function integration_case(order)
    u = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
    order > 0 && (u = u + variables(6; order) / 1000)
    su = SVector{6}(u)
    check = order == 0 ? ((a, b) -> maximum(abs, a - b) < 1.0e-12) :
        ((a, b) -> maximum(coefficient_norm, a - b) < 1.0e-12)
    return compare(
        "RK4 state: order $order, 64 steps",
        () -> rk4(u, 2π / 64, 64), () -> rk4(su, 2π / 64, 64), check
    )
end

function main()
    rows = Dict[]
    for order in (2, 3, 5)
        push!(rows, coefficient_case(order))
        push!(rows, storage_case(order))
    end
    for (n, order) in ((2, 6), (6, 3))
        push!(rows, evaluation_case(n, order))
    end
    for order in (0, 2, 3, 5)
        push!(rows, integration_case(order))
    end
    return if !isempty(ARGS)
        open(last(ARGS), "w") do io
            TOML.print(
                io, Dict(
                    "julia" => string(VERSION), "staticarrays" => string(pkgversion(StaticArrays)),
                    "cpu" => Sys.CPU_NAME, "kernel" => string(Sys.KERNEL), "cases" => rows
                )
            )
        end
    end
end
abspath(PROGRAM_FILE) == (@__FILE__) && main()
end
